# -*- coding: utf-8 -*-
"""
 * @Date: 2024-12-25 12:06:26
 * @Editors: Jessica_Bryant jessawbryant@gmail.com
* @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
* @LastEditTime: 2025-07-09 19:22:33
* @FilePath: /genome/genome/gene_statistic.py
 * @Description:

This module holds dictionaries that contain codon and amino acid tables, and functions that
 calculate GC, N-ARSC, C-ARSC and Nc.

Code is modified from https://github.com/faylward/pangenomics/blob/7dbc8269f1618492dd4ebf5752c8c4d3deca0e43/get_ARSC.py

 Minor modifications made by Frank Aylward in 5/29/2017 <faylward@vt.edu>
"""

import itertools
import math
import sys
from typing import Callable, Container, Final, Iterable, NamedTuple
from Bio import SeqUtils
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


class CodonTable:
    class CodonFeat(NamedTuple):
        aa: Seq
        n_codon: int
        gc_rank: float
        family: Seq

    def __init__(self, table: int | str = 11):
        """
        Initialize the CodonTable with codon-to-amino acid mappings and codon family groupings for a specified genetic code table.
        
        Constructs mappings from codons to amino acids, identifies stop codons, and organizes codons into families based on amino acid, codon count, and GC content rank. Also builds pseudo codon families for further codon usage analysis.
        
        Parameters:
            table (int | str): NCBI translation table number or name specifying the genetic code (default is 11).
        """
        self.table = table
        _table = {
            i: i.translate(table)  # type: ignore[reportArgumentType]
            for i in (Seq("".join(i)) for i in itertools.product("ACGT", repeat=3))
        }
        aa2codon: dict[Seq, set[Seq]] = {}
        for c, a in _table.items():
            aa2codon.setdefault(a, set()).add(c)
        _table2: dict[Seq, tuple[Seq, int, float]] = {}
        stop_codons: set[Seq] = set()
        sf_aa_codon: dict[int, dict[Seq, list[set[Seq]]]] = {}
        for a, cs in aa2codon.items():
            if a == "*":
                stop_codons.update(cs)
                continue
            n_codon = len(cs)
            cgc = {c: SeqUtils.gc_fraction(c) for c in cs}
            if n_codon == 1:
                gc_rank = {c: math.nan for c in cgc.values()}
            else:
                gc_rank = {c: i for i, c in enumerate(sorted(set(cgc.values())))}
            if n_codon <= 4:
                for c in cs:
                    _table2[c] = a, n_codon, (gc_rank[cgc[c]] / len(cgc))
                sf_aa_codon.setdefault(n_codon, {})[a] = [cs]
            else:
                # split into two more groups
                sub_family: dict[Seq, set[Seq]] = {}
                for c in cs:
                    sub_family.setdefault(c[:2], set()).add(c)
                for cs2 in sub_family.values():
                    for c in cs2:
                        _table2[c] = a, len(cs2), (gc_rank[cgc[c]] / len(cgc))
                    sf_aa_codon.setdefault(len(cs2), {}).setdefault(a, []).append(cs2)
        self.codon_table: Final = {
            k: self.CodonFeat(*v, k[:2]) for k, v in _table2.items()
        }
        self.aa2codon: Final = aa2codon
        self.stop_codons: Final = frozenset(stop_codons)
        # pseudo codon family
        self.pcf: Final = sf_aa_codon

    @classmethod
    def get(cls, table: int | str = 11, _table_cache={}) -> "CodonTable":
        """
        Return a cached CodonTable instance for the specified translation table.
        
        If the CodonTable for the given table number or name does not exist in the cache, it is created and stored for future use.
        
        Parameters:
            table (int | str): NCBI translation table number or name.
        
        Returns:
            CodonTable: Cached CodonTable instance corresponding to the specified translation table.
        """
        if table not in _table_cache:
            _table_cache[table] = cls(table)
        return _table_cache[table]

    def parse(
        self,
        seq_cds: Seq,
        transl_except: Container[int] = frozenset(),
        errorfile_handle=sys.stderr,
        id: int | str | None = "*",
    ):
        """
        Parses a coding DNA sequence (excluding the start codon) to count codon usage and assess GC variability.
        
        Checks for stop codons, incomplete codons, ambiguous bases ('N'), and unknown codons, raising a ValueError with details if an error is encountered. Returns a `_CodonUsage` object containing codon usage counts, translation table, gene codon length, and GC variability statistics.
        
        Parameters:
            seq_cds (Seq): Coding DNA sequence to parse, excluding the start codon.
            transl_except (Container[int], optional): Codon positions (by index) to ignore stop codon errors, typically for translation exceptions.
            errorfile_handle: File-like object for error output.
            id (int | str | None, optional): Identifier for the sequence, used in error messages.
        
        Returns:
            _CodonUsage: Codon usage statistics and GC variability metrics for the parsed sequence.
        
        Raises:
            ValueError: If a stop codon (not in transl_except), incomplete codon, or ambiguous base is encountered.
        """
        codon_used = {c: 0 for c in self.codon_table}
        gene_codon_length = 0
        var_gc_aanum = 0
        var_gc_sumrank = 0.0
        # Parse through codons
        # remove start and stop codon from current DNA sequence
        err_str = ""
        for i in range(3, len(seq_cds), 3):
            # count each codon
            current_codon = seq_cds[i : (i + 3)]

            # identify possible codons that will cause problems
            if err_str:
                pass
            elif current_codon in self.stop_codons:
                if (i // 3) not in transl_except:
                    err_str = f"after stop codon {current_codon}"
                continue
            elif len(current_codon) < 3:
                err_str = f"incomplete codon {current_codon}"
                continue
            elif "N" in current_codon:
                err_str = "N present"
            elif current_codon not in self.codon_table:
                # err_str = "unknown codon"
                codon_used[current_codon] = codon_used.get(current_codon, 0) + 1
                continue
            if err_str:
                print(id, current_codon, i, err_str, sep="\t", file=errorfile_handle)
                raise ValueError(f"Error in {id} at {i}({current_codon}): {err_str}")

            codon_used[current_codon] += 1
            gene_codon_length += 1
            # count the number of codons with GC variability and their rank
            if not math.isnan(self.codon_table[current_codon].gc_rank):
                var_gc_aanum += 1
                var_gc_sumrank += self.codon_table[current_codon].gc_rank
        return self._CodonUsage(
            codon_used, self.table, gene_codon_length, var_gc_aanum, var_gc_sumrank
        )

    @classmethod
    def get_parse(
        cls, *_seq_cds: Seq | SeqRecord, errorfile_handle=sys.stderr, table=None
    ):
        """
        Parses multiple coding DNA sequences or SeqRecords to aggregate codon usage statistics.
        
        Parameters:
            *_seq_cds: One or more coding DNA sequences (Seq or SeqRecord) to be parsed.
            errorfile_handle: File-like object for error reporting (default: sys.stderr).
            table: Optional translation table identifier to use if not specified in SeqRecord annotations.
        
        Returns:
            _CodonUsage: An aggregated codon usage statistics object for all input sequences.
        
        Raises:
            ValueError: If multiple translation tables are detected among the input sequences.
        """
        id2seqs: dict[str | int, Seq] = {}
        id2except: dict[str | int, set[int]] = {}
        code_tables = set()
        id: str | int
        for id, seq_cds in enumerate(_seq_cds):
            if isinstance(seq_cds, SeqRecord):
                id = seq_cds.id or id
                seq = seq_cds.seq[check_frame(seq_cds.annotations) :]
                code_tables.add(seq_cds.annotations.get("transl_table", table))
                id2except[id] = {
                    int(i.split("@")[0])
                    for i in str(seq_cds.annotations.get("transl_except", "")).split(
                        ";"
                    )
                    if i
                }
            else:
                seq = seq_cds
                code_tables.add(table)
                id2except[id] = set()
            id2seqs[id] = seq
        code_tables -= {None}
        if len(code_tables) > 1:
            raise ValueError(f"Multiple translation tables found: {code_tables}")
        code_table = cls.get() if not code_tables else cls.get(list(code_tables)[0])

        return cls._CodonUsage.concat(
            (
                code_table.parse(seq, id2except[id], errorfile_handle, id)
                for id, seq in id2seqs.items()
            ),
            table=code_table.table,
        )

    class _CodonUsage(NamedTuple):
        codon_used: dict[Seq, int]
        table: int | str = 11
        gene_codon_length: int = 0
        var_gc_aanum: int = 0
        var_gc_sumrank: float = 0

        @property
        def gc_variability(self):
            # calculate GC rank sums
            """
            Calculate the average GC rank variability for codons with GC variability.
            
            Returns:
                float: The mean GC rank among codons with variable GC content, or NaN if none are present.
            """
            gc_variability = (
                float(self.var_gc_sumrank) / float(self.var_gc_aanum)
                if self.var_gc_aanum != 0
                else math.nan
            )
            return gc_variability

        def getCFs(self, m: int):
            """
            Yield pseudo codon families of a specified size with usage statistics.
            
            Parameters:
                m (int): The size of the codon family to retrieve.
            
            Yields:
                CodonTable.PseudoCodonFamily: An object containing the total count, F value, codon set, and amino acid for each pseudo codon family of size m.
            """
            code_table = CodonTable.get(self.table)
            for aa in code_table.pcf[m]:
                for cs in code_table.pcf[m][aa]:
                    codon_n = [self.codon_used[codon] for codon in cs]
                    codon_sum = sum(codon_n)
                    yield CodonTable.PseudoCodonFamily(
                        sum(codon_n),
                        sum(((i + 1) / (codon_sum + m)) ** 2 for i in codon_n),
                        frozenset(cs),
                        aa,
                    )

        def wighted_F(self, m: int):
            """
            Calculate the weighted F value for codon families of size `m`.
            
            The weighted F value is used in codon usage analysis to assess codon bias within families of synonymous codons of a given size. If no codons are present in families of size `m`, returns `m`. For singleton families (`m == 1`), returns 1.
            
            Parameters:
                m (int): The size of the codon family.
            
            Returns:
                float: The weighted F value for codon families of size `m`.
            """
            if m == 1:
                return 1
            cfs = list(self.getCFs(m))
            sum_n = sum(cf.n for cf in cfs)
            if sum_n == 0:
                return m
            return sum_n / sum(cf.n * cf.F for cf in cfs)

        @property
        def SCU(self):
            """
            Calculate the effective number of codons (Nc) for the coding sequence using weighted codon family statistics.
            
            Returns:
                float: The effective number of codons (Nc) as defined in the referenced Molecular Biology and Evolution publication.
            """
            pcf = CodonTable.get(self.table).pcf
            return sum(
                sum(len(i) for i in aa_s) * self.wighted_F(i) for i, aa_s in pcf.items()
            )

        @classmethod
        def concat(cls, cus: Iterable["CodonTable._CodonUsage"], table=None):
            """
            Aggregate multiple `_CodonUsage` instances into a single summary object.
            
            Parameters:
            	cus (Iterable[CodonTable._CodonUsage]): An iterable of `_CodonUsage` objects to be combined.
            	table (optional): The translation table to use. If not provided, all input objects must share the same table.
            
            Returns:
            	CodonTable._CodonUsage: A new `_CodonUsage` instance with summed codon usage counts and statistics from all inputs.
            
            Raises:
            	ValueError: If input objects use different translation tables and no table is specified.
            """
            if table is None:
                cus = list(cus)
                code_tables = set(cu.table for cu in cus)
                if len(code_tables) != 1:
                    raise ValueError(
                        f"conflict translation tables found: {code_tables}"
                    )
                table = list(code_tables)[0]
            code_table = CodonTable.get(table)
            codon_used = {c: 0 for c in code_table.codon_table}
            gene_codon_length = 0
            var_gc_aanum = 0
            var_gc_sumrank = 0.0
            for csf in cus:
                for c in codon_used:
                    codon_used[c] += csf.codon_used[c]
                gene_codon_length += csf.gene_codon_length
                var_gc_aanum += csf.var_gc_aanum
                var_gc_sumrank += csf.var_gc_sumrank
            return cls(
                codon_used,
                code_table.table,
                gene_codon_length,
                var_gc_aanum,
                var_gc_sumrank,
            )

    class PseudoCodonFamily(NamedTuple):
        n: int
        F: float
        codons: frozenset[Seq]
        aa: Seq


class ARSC(NamedTuple):
    C: float | int = 0
    N: float | int = 0
    S: float | int = 0
    mw: float = 0  # molecular_weight
    len: int = 0

    @staticmethod
    def AA_DICT() -> dict[str, "ARSC"]:
        """
        Return a dictionary mapping amino acid single-letter codes to their ARSC (carbon, nitrogen, sulfur counts) values.
        
        Returns:
            dict: A dictionary where keys are amino acid codes and values are ARSC instances representing elemental composition.
        """
        return {
            "K": ARSC(C=4, N=1, S=0),
            "R": ARSC(C=4, N=3, S=0),
            "H": ARSC(C=4, N=2, S=0),
            "D": ARSC(C=2, N=0, S=0),
            "E": ARSC(C=3, N=0, S=0),
            "N": ARSC(C=2, N=1, S=0),
            "Q": ARSC(C=3, N=1, S=0),
            "S": ARSC(C=1, N=0, S=0),
            "T": ARSC(C=2, N=0, S=0),
            "Y": ARSC(C=7, N=0, S=0),
            "A": ARSC(C=1, N=0, S=0),
            "V": ARSC(C=3, N=0, S=0),
            "L": ARSC(C=4, N=0, S=0),
            "I": ARSC(C=4, N=0, S=0),
            "P": ARSC(C=3, N=0, S=0),
            "F": ARSC(C=7, N=0, S=0),
            "M": ARSC(C=3, N=0, S=1),
            "W": ARSC(C=9, N=1, S=0),
            "G": ARSC(C=0, N=0, S=0),
            "C": ARSC(C=1, N=0, S=1),
            "U": ARSC(C=1, N=0, S=0),
            "J": ARSC(C=4, N=0, S=0),
            "B": ARSC(C=2, N=0.5, S=0),
            "Z": ARSC(C=3, N=0.5, S=0),
            "O": ARSC(C=10, N=2, S=0),
        }

    @classmethod
    def parse(cls, seq_aa: Seq):
        """
        Parses an amino acid sequence to compute total counts of carbon (C), nitrogen (N), and sulfur (S) atoms, as well as molecular weight and sequence length, excluding stop codons.
        
        Parameters:
        	seq_aa (Seq): Amino acid sequence in single-letter code.
        
        Returns:
        	ARSC: An instance containing total C, N, S counts, molecular weight, and length of the sequence (excluding stop codons).
        """
        aa_dict = cls.AA_DICT()

        # remove invalid char from the string
        stop_pos = seq_aa.rfind("*")
        seq_aa_no_stop = seq_aa[:stop_pos] if stop_pos != -1 else seq_aa
        total_molecular_weight = SeqUtils.molecular_weight(
            seq_aa_no_stop, seq_type="protein"
        )

        # calculate ARSC and Molecular Weight
        arsc_c = sum(aa_dict[x].C for x in seq_aa_no_stop)
        arsc_n = sum(aa_dict[x].N for x in seq_aa_no_stop)
        arsc_s = sum(aa_dict[x].S for x in seq_aa_no_stop)

        valid_len = len(seq_aa_no_stop)

        return cls(
            C=arsc_c, N=arsc_n, S=arsc_s, mw=total_molecular_weight, len=valid_len
        )

    def __add__(self, other):
        """
        Add two ARSC instances by summing their respective fields.
        
        If the other operand is not an ARSC instance, returns the concatenation of the ARSC tuple with the other object.
        """
        if isinstance(other, self.__class__):
            s_d, o_d = self._asdict(), other._asdict()
            return self.__class__(**{k: s_d[k] + o_d[k] for k in s_d})
        return tuple(self) + other

    def scale(self):
        """
        Return a new ARSC instance with all values normalized by sequence length.
        
        If the sequence length is zero, returns the original instance unchanged.
        """
        if self.len == 0:
            return self
        return self.__class__(
            self.C / self.len,
            self.N / self.len,
            self.S / self.len,
            self.mw / self.len,
            1,
        )

    @classmethod
    def concat(cls, arscs: Iterable["ARSC"]):
        """
        Aggregate multiple ARSC instances by summing their elemental counts, molecular weights, and sequence lengths.
        
        Parameters:
            arscs (Iterable[ARSC]): An iterable of ARSC instances to be combined.
        
        Returns:
            ARSC: A new ARSC instance representing the sum of all input ARSCs.
        """
        arsc = cls()
        for arsc_ in arscs:
            arsc += arsc_
        return arsc


def aa_mw(seq_aa: Seq):
    """
    Calculate the average molecular weight of an amino acid sequence, excluding any terminal stop codon.
    
    Parameters:
    	seq_aa (Seq): Amino acid sequence as a Biopython Seq object.
    
    Returns:
    	float: Average molecular weight per residue, excluding the stop codon if present.
    """
    # remove whitespaces and stop codon at end of sequences
    stop_pos = seq_aa.rfind("*")
    seq_aa_no_stop = seq_aa[:stop_pos] if stop_pos != -1 else seq_aa
    total_molecular_weight = SeqUtils.molecular_weight(
        seq_aa_no_stop, seq_type="protein"
    )
    av_molecular_weight = total_molecular_weight / len(seq_aa_no_stop)

    return av_molecular_weight


from .bin_statistic import _BinStatisticContainer
from .gff import extract, _translate
from .translate import check_frame


class GeneStat(NamedTuple):
    codon_usage: CodonTable._CodonUsage
    arsc: ARSC
    n_gene: int

    @classmethod
    def parse(
        cls,
        seq_iter: Iterable[SeqRecord],
        min_contig_len=0,
        min_aa_length=33,
        call_gene_id="infer_gene_id",
    ):
        """
        Parses an iterable of SeqRecord objects to extract gene statistics for coding sequences (CDS).
        
        For each sequence, extracts CDS features, filters by minimum contig and amino acid length, translates CDS to amino acids, and skips pseudogenes or sequences with ambiguous bases. Computes codon usage and amino acid residue side chain (ARSC) statistics for each valid gene.
        
        Parameters:
            seq_iter (Iterable[SeqRecord]): Iterable of sequence records to process.
            min_contig_len (int, optional): Minimum nucleotide length for a contig to be considered.
            min_aa_length (int, optional): Minimum amino acid length for a CDS to be included.
            call_gene_id (str, optional): Method or attribute to use for determining gene IDs.
        
        Returns:
            dict[tuple[str, str], GeneStat]: Dictionary mapping (sequence ID, gene ID) to GeneStat instances containing codon usage and ARSC statistics for each gene.
        """
        seqaa2stat: dict[tuple[str, str], GeneStat] = {}
        for seq in seq_iter:
            if len(seq.seq) < min_contig_len:
                continue
            assert seq.id
            for cds in extract(
                [seq],
                fet_type="CDS",
                translate=False,
                call_gene_id=call_gene_id,
                min_aa_length=min_aa_length,
            ):
                assert cds.id
                aa = _translate(cds).seq
                if cds.features and cds.features[-1].qualifiers.get("pseudo") == [
                    "true"
                ]:
                    try:
                        csf = CodonTable.get_parse(cds)
                        arsc = ARSC.parse(aa)
                    except ValueError:
                        continue
                elif "N" in cds.seq or "X" in aa:
                    continue
                else:
                    csf = CodonTable.get_parse(cds)
                    arsc = ARSC.parse(aa)
                seqaa2stat[seq.id, cds.id] = cls(csf, arsc, 1)
        return seqaa2stat


class GeneStatisticContainer(_BinStatisticContainer):
    def __init__(
        self,
        seqiter: Iterable[SeqRecord],
        source_file,
        min_contig_len=0,
        loader: Callable[
            [Iterable[SeqRecord], int, int], dict[tuple[str, str], GeneStat]
        ] = GeneStat.parse,
        min_aa_len=33,
    ):
        """
        Initializes the GeneStatisticContainer by loading gene statistics from a sequence iterable.
        
        Parameters:
            seqiter (Iterable[SeqRecord]): An iterable of sequence records containing gene data.
            source_file: Reference to the source file from which sequences are loaded.
            min_contig_len (int, optional): Minimum contig length required for a sequence to be included. Defaults to 0.
            loader (Callable, optional): Function to parse sequence records and extract gene statistics. Defaults to GeneStat.parse.
            min_aa_len (int, optional): Minimum amino acid length required for a gene to be included. Defaults to 33.
        """
        self._aa_stats = loader(seqiter, min_contig_len, min_aa_len)
        self.source_file = source_file

    def statistic(self):
        """
        Aggregate codon usage and amino acid statistics across all genes in the container.
        
        Returns:
            GeneStatistic: A named tuple containing the effective number of codons (SCU), GC variability, translation table, scaled ARSC values (C, N, S), average protein molecular weight, average protein length, and total gene count.
        """
        csf_all = CodonTable._CodonUsage.concat(
            (seq_cds.codon_usage for seq_cds in self._aa_stats.values())
        )
        arsc = ARSC.concat((seq_cds.arsc for seq_cds in self._aa_stats.values()))
        arsc_scale = arsc.scale()
        genes_num = len(self._aa_stats)
        return self.GeneStatistic(
            scu=csf_all.SCU,
            gc_variability=csf_all.gc_variability,
            table=csf_all.table,
            C_ARSC=arsc_scale.C,
            N_ARSC=arsc_scale.N,
            S_ARSC=arsc_scale.S,
            avg_protein_mw=arsc.mw / genes_num,
            avg_protein_len=arsc.len / genes_num,
            genes_num=genes_num,
        )

    class GeneStatistic(NamedTuple):
        scu: float
        gc_variability: float
        table: str | int
        C_ARSC: float
        N_ARSC: float
        S_ARSC: float
        avg_protein_mw: float
        avg_protein_len: float
        genes_num: int
