# -*- coding: utf-8 -*-
"""
 * @Date: 2024-12-25 12:06:26
 * @Editors: Jessica_Bryant jessawbryant@gmail.com
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2025-02-16 11:21:38
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
from typing import Callable, Final, Iterable, NamedTuple
from Bio import SeqUtils
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


class CodonTable:
    class CodonFeat(NamedTuple):
        aa: Seq
        n_codon: int
        gc_rank: float
        family: Seq

    def __init__(self, table=11):
        self.table = table
        _table = {
            i: i.translate(table)
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
    def get(cls, table=11, _table_cache={}) -> "CodonTable":
        if table not in _table_cache:
            _table_cache[table] = cls(table)
        return _table_cache[table]

    def parse(
        self,
        seq_cds: Seq,
        transl_except: set[int] = set(),
        errorfile_handle=sys.stderr,
        id="*",
    ):
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
        if not code_tables:
            code_table = cls.get()
        else:
            code_table = cls.get(list(code_tables)[0])
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
            gc_variability = (
                float(self.var_gc_sumrank) / float(self.var_gc_aanum)
                if self.var_gc_aanum != 0
                else math.nan
            )
            return gc_variability

        def getCFs(self, m: int):
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
            Reference:

              An Improved Implementation of Effective Number of Codons (Nc)

              Molecular Biology and Evolution, Volume 30, Issue 1, January 2013, Pages 191–196, https://doi.org/10.1093/molbev/mss201
            """
            pcf = CodonTable.get(self.table).pcf
            return sum(
                sum(len(i) for i in aa_s) * self.wighted_F(i) for i, aa_s in pcf.items()
            )

        @classmethod
        def concat(cls, cus: Iterable["CodonTable._CodonUsage"], table=None):
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
            "U": ARSC(C=1, N=0, S=0),
            "O": ARSC(C=10, N=2, S=0),
        }

    @classmethod
    def parse(cls, seq_aa: Seq):
        """
        This functions takes an amino acid sequence coded in single letters and returns N/C ARSC and Molecular Weight
        N and S counts come from page 30 of 'Understanding Bioinformatics' by Zvelbil and Baum
        molecular weights from http://www.webqc.org/aminoacids.php
        """
        aa_dict = cls.AA_DICT()

        # remove invalid char from the string
        seq_aa_no_stop = seq_aa[: seq_aa.rfind("*")]
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
        if isinstance(other, self.__class__):
            s_d, o_d = self._asdict(), other._asdict()
            return self.__class__(**{k: s_d[k] + o_d[k] for k in s_d.keys()})
        return tuple(self) + other

    def scale(self):
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
        arsc = cls()
        for arsc_ in arscs:
            arsc += arsc_
        return arsc


def aa_mw(seq_aa: Seq):
    """
    This functions takes a protein sequence coded in nucleotides and returns N-ARSC.
    Internal stop codons will cause this script problems
    """
    # remove whitespaces and stop codon at end of sequences
    seq_aa_no_stop = seq_aa[: seq_aa.rfind("*")]
    total_molecular_weight = SeqUtils.molecular_weight(
        seq_aa_no_stop, seq_type="protein"
    )
    av_molecular_weight = total_molecular_weight / len(seq_aa_no_stop)

    return av_molecular_weight


from .bin_statistic import _BinStatisticContainer
from .gff import extract, _translate, check_frame


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
        self._aa_stats = loader(seqiter, min_contig_len, min_aa_len)
        self.source_file = source_file

    def statistic(self):
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
            avg_protien_len=arsc.len / genes_num,
            genes_num=genes_num,
        )

    class GeneStatistic(NamedTuple):
        scu: float
        gc_variability: float
        table: float
        C_ARSC: float
        N_ARSC: float
        S_ARSC: float
        avg_protein_mw: float
        avg_protien_len: float
        genes_num: int
