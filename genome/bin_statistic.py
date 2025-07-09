# -*- coding: utf-8 -*-
"""
* @Date: 2022-10-15 17:05:11
* @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
* @LastEditTime: 2025-07-09 19:27:35
* @FilePath: /genome/genome/bin_statistic.py
* @Description:
"""


import math
import pickle
import shutil
from functools import cached_property
from pathlib import Path
from typing import Callable, Final, Iterable, NamedTuple, Union

import numpy as np
import pandas as pd
from Bio import SeqIO, SeqRecord
from numpy import mean

from .gff import Parse

PathLike = Union[str, Path]


class Contig2Bin:
    """
    Usage: Contig2Bin(contig2bin_tsv, contigs).output(outdir)

    will read a table such without head such as:
        (format: contig\tbin)
        ```
        k141_128353\tconcoct_0
        k141_15265\tconcoct_0
        k141_164694\tconcoct_0
        k141_172306\tconcoct_0
        k141_18944\tconcoct_0
        ```
    """

    columns: Final = ["contig", "bin"]

    def __init__(self, contig2bin_tsv, contigs, outdir=None):
        """
        Initialize a Contig2Bin instance by parsing contig-to-bin mappings and contig sequences.
        
        Parameters:
            contig2bin_tsv: Path, pandas DataFrame, or Series specifying the mapping of contig IDs to bin names.
            contigs: Path to a FASTA file or an iterable of SeqRecord objects representing contig sequences.
            outdir: Optional; output directory for writing binned sequence files.
        """
        self.contig2bin_tsv = self.parse_contig2bin(contig2bin_tsv)
        self._contigs = self.parse_contigs(contigs)
        self._outdir = None
        if outdir:
            self.outdir = outdir

    @classmethod
    def parse_contig2bin(cls, contig2bin_tsv: PathLike | pd.DataFrame | pd.Series):
        """
        Parses a contig-to-bin mapping table from a TSV file, pandas DataFrame, or Series into a standardized DataFrame.
        
        Supports both standard and VAMB-formatted tables, normalizing contig IDs and bin names for downstream processing.
        
        Parameters:
        	contig2bin_tsv: Input mapping as a file path, DataFrame, or Series. Accepts VAMB-format tables with contig and bin columns.
        
        Returns:
        	DataFrame indexed by normalized contig IDs with a single 'bin' column.
        """
        if isinstance(contig2bin_tsv, pd.Series):
            pd_raw = contig2bin_tsv.reset_index()
            pd_raw.columns = cls.columns
        elif isinstance(contig2bin_tsv, pd.DataFrame):
            if set(contig2bin_tsv.columns) >= set(cls.columns):
                pd_raw = contig2bin_tsv[cls.columns]
            elif len(contig2bin_tsv.columns) == 1:
                pd_raw = contig2bin_tsv.reset_index()
                pd_raw.columns = cls.columns
            else:
                raise ValueError(f"columns of table should be {cls.columns}")
        else:
            pd_raw = pd.read_csv(contig2bin_tsv, sep="\t", names=cls.columns)
        return pd_raw.assign(
            contig_id=lambda df: df["contig"].apply(lambda x: x.strip().split()[0]),
            bin=lambda df: df["bin"].apply(lambda x: str(x).strip()),
        ).set_index("contig_id")[["bin"]]

    @cached_property
    def bin2seqs(self):
        """
        Return a dictionary mapping each bin name to its contig sequences.
        
        Returns:
            bin2seqs (dict): A dictionary where each key is a bin name and each value is a dictionary mapping contig names to their corresponding SeqRecord objects.
        """
        bin2seqs: dict[str, dict[str, SeqRecord.SeqRecord]] = {
            b: {} for b in self.contig2bin_tsv["bin"].unique()
        }
        for i in self.contigs:
            if i.name in self.contig2bin_tsv.index:
                bin2seqs[str(self.contig2bin_tsv.loc[i.name, "bin"])][i.name] = i
        return bin2seqs

    def parse_contigs(self, contigs: PathLike | Iterable[SeqRecord.SeqRecord]):
        """
        Parses contig sequences from a file path or an iterable of SeqRecord objects.
        
        If a file path is provided, attempts to parse it as a FASTA file; if unsuccessful, falls back to a GFF parser. If an iterable is provided, returns a callable that yields the contained SeqRecords.
        
        Parameters:
            contigs: Path to a contig file or an iterable of SeqRecord objects.
        
        Returns:
            A callable that yields SeqRecord objects when invoked.
        
        Raises:
            ValueError: If the input is neither a valid path nor an iterable of SeqRecord objects.
        """
        if isinstance(contigs, (str, Path)):
            try:
                next(SeqIO.parse(contigs, "fasta"))
                return lambda: SeqIO.parse(contigs, "fasta")
            except ValueError:
                return Parse(contigs)
        elif isinstance(contigs, Iterable):
            return lambda: list(contigs)
        else:
            raise ValueError("contigs should be a path or a iterable of SeqRecord")

    @property
    def contigs(self) -> Iterable[SeqRecord.SeqRecord]:
        """
        Returns an iterable of contig sequences associated with the current instance.
        """
        return self._contigs()

    @property
    def outdir(self):
        """
        Returns the output directory path for writing binned sequence files.
        
        Raises:
            ValueError: If the output directory has not been set.
        """
        if self._outdir is None:
            raise ValueError(".outdir not set")
        return self._outdir

    @outdir.setter
    def outdir(self, outdir: Path):
        """
        Sets the output directory for bin files and creates the directory if it does not exist.
        
        Parameters:
            outdir (Path): Path to the desired output directory.
        
        Returns:
            Path: The created or existing output directory path.
        """
        td = Path(outdir)
        td.mkdir(parents=True, exist_ok=True)
        self._outdir = outdir
        return td

    def __call__(self, outdir: PathLike):
        """
        Writes all binned contig sequences to separate FASTA files in the specified output directory.
        
        This method efficiently writes each bin's sequences to disk in one operation and returns a `Binput` instance representing the output directory and bins.
        
        Parameters:
            outdir: Path to the directory where bin FASTA files will be written.
        
        Returns:
            Binput: An object representing the output directory, bin IDs, and file suffix.
        """
        self.outdir = Path(outdir)
        bin2seqs = self.bin2seqs
        for b, seqs in bin2seqs.items():
            SeqIO.write(seqs.values(), self.outdir / f"{b}.fa", "fasta-2line")
        return self.output

    @property
    def output(self):
        """
        Returns a `Binput` instance representing the expected output directory and file structure for binned sequences, regardless of whether the files currently exist on disk.
        """
        return Binput(self.outdir, list(self.bin2seqs), ".fa")

    def extract1(self, bin_name: str | int):
        """
        Writes the sequences of a specified bin to a FASTA file in the output directory.
        
        Parameters:
        	bin_name (str | int): The name of the bin, or its integer index in the unique bin list.
        
        Returns:
        	Path: The path to the generated FASTA file containing the bin's sequences.
        """
        if isinstance(bin_name, int):
            b = self.contig2bin_tsv["bin"].unique()[bin_name]
        else:
            b = bin_name
        bout = self.outdir / f"{b}.fa"
        SeqIO.write(self.bin2seqs[b].values(), bout, "fasta-2line")
        return bout

    def extract(self, *bin_names: str | int):
        """
        Generator yielding FASTA file paths for the specified bin names.
        
        Parameters:
        	bin_names: One or more bin names or indices to extract.
        
        Returns:
        	Generator yielding the paths to the extracted FASTA files for each specified bin.
        """
        return (self.extract(b) for b in bin_names)


class Binput(NamedTuple):
    bindir: Path
    binids: list[str]
    suffix: str

    @classmethod
    def parse(
        cls,
        bin_output: PathLike,
        bin_input: PathLike,
        support: Union[PathLike, str],
        keep_if_avail=True,
    ):
        """
        Creates a directory of binned sequence files from an input directory or single file.
        
        If the input is a directory containing files with the specified suffix, it either keeps the directory as-is (if `keep_if_avail` is True and the suffix ends with ".fa") or copies files to the output directory with a ".fa" suffix. If the input is a single file, it generates binned FASTA files using the provided contig-to-bin mapping.
        
        Parameters:
            bin_output: Path to the output directory for binned sequences.
            bin_input: Path to the input directory or file containing bin sequences.
            support: Suffix string or path to a contig-to-bin mapping file.
            keep_if_avail (bool): If True and input directory files end with ".fa", the input directory is used directly.
        
        Returns:
            An instance of the class representing the binned sequence directory.
        """
        bin_input = Path(bin_input)
        (bin_output_ := Path(bin_output)).mkdir(parents=True, exist_ok=True)
        if bin_input.is_dir():
            assert list(bin_input.glob(f"*{support}")), "input is not a valid bin path"
            suffix = str(support)
            support_str_len = len(suffix)
            if suffix.endswith(".fa") and keep_if_avail:
                binids = [i.stem for i in bin_input.glob(f"*{support}")]
                bindir = bin_input
            else:
                suffix = ".fa"
                bindir = bin_output_
                binids = []
                for bin_file in bin_input.glob(f"*{support}"):
                    binids.append(
                        bin_name := bin_file.name[:-support_str_len].rstrip(".")
                    )
                    shutil.copy(bin_file, bindir / f"{bin_name}.fa")
            return cls(bindir, binids, suffix)
        else:
            assert bin_input.is_file() and Path(support).is_file()
            return Contig2Bin(bin_input, support)(bin_output_)

    def fas(self):
        """
        Yield paths to FASTA files for each bin in the binned sequence directory.
        
        Returns:
            Generator of file paths corresponding to each bin's FASTA file.
        """
        return (self.bindir / f"{i}{self.suffix}" for i in self.binids)

    def fas_with(self, suffix: str):
        """
        Return a generator of FASTA file paths for each bin using the specified file suffix.
        
        Parameters:
        	suffix (str): The file extension or suffix to use for the FASTA files.
        
        Returns:
        	Generator[Path]: Paths to FASTA files for each bin with the given suffix.
        """
        return Binput(self.bindir, self.binids, suffix).fas()


format_bin_input = Binput.parse


def contig2bin(outdir: PathLike, contig2bin_tsv: PathLike, contigs: PathLike):
    """
    Deprecated wrapper for extracting binned sequences from contigs using a contig-to-bin mapping.
    
    Writes sequences for each bin to separate FASTA files in the specified output directory. Use the `Contig2Bin` class directly for new code.
    """
    import warnings

    warnings.warn(
        "contig2bin() is deprecated; use the Contig2Bin class instead",
        DeprecationWarning,
        stacklevel=2,
    )
    return Contig2Bin(contig2bin_tsv, contigs)(outdir)


def calculateN50(seqLens: list[int]):
    """
    Calculate the N50 statistic for a list of sequence lengths.
    
    The N50 is the length of the shortest sequence at which the sum of lengths of all longer or equal sequences covers at least half of the total length.
    
    Parameters:
        seqLens (list[int]): List of sequence lengths.
    
    Returns:
        int: The N50 value.
    """
    seqLens_ = sorted(seqLens, reverse=True)
    thresholdN50 = sum(seqLens_) / 2.0

    testSum = 0
    for seqLen in seqLens_:
        testSum += seqLen
        if testSum >= thresholdN50:
            N50 = seqLen
            break

    return N50


class SeqStat(NamedTuple):
    len: int = 0
    gc_pct: float = 0.0
    gc: int = 0
    at: int = 0
    n: int = 0
    n_cds: int = 0
    len_cds: int = 0

    @classmethod
    def parse(
        cls, seq_iter: Iterable[SeqRecord.SeqRecord], min_contig_len=0, min_aa_len=33
    ):
        """
        Parse sequence records and compute statistics for each, including length, GC content, ambiguous bases, and coding sequence features.
        
        Parameters:
            seq_iter (Iterable[SeqRecord.SeqRecord]): Iterable of sequence records to analyze.
            min_contig_len (int, optional): Minimum contig length to include in the statistics. Defaults to 0.
            min_aa_len (int, optional): Minimum amino acid length for coding sequences (CDS) to be counted. Defaults to 33.
        
        Returns:
            dict[str, SeqStat]: Dictionary mapping sequence IDs to their computed statistics.
        """
        _seq_stats: dict[str, SeqStat] = {}
        min_gene_len = int(min_aa_len) * 3
        for rec_enum, rec in enumerate(seq_iter):
            if len(rec.seq) < min_contig_len:
                continue
            sequence_arr = np.frombuffer(bytes(rec.seq).upper(), dtype="S1")
            base_count = {
                b.decode(): int(n)
                for b, n in zip(*np.unique(sequence_arr, return_counts=True))
            }
            at = sum(base_count.get(i, 0) for i in "ATWU")
            gc = sum(base_count.get(i, 0) for i in "CGS")

            gc_content = 0.0 if (gcat := gc + at) <= 0 else float(gc) / gcat

            cds_mask = np.zeros(len(rec.seq))
            n_cds = 0
            for fet in rec.features:
                if fet.type == "CDS" and fet.location is not None:
                    if len(fet) < min_gene_len:
                        continue
                    cds_mask[fet.location.start : fet.location.end] = 1
                    n_cds += 1
            seq_n = sequence_arr == b"N"
            cds_mask[seq_n] = 0
            _seq_stats[rec.id or str(rec_enum)] = cls(
                len(rec),
                gc_content,
                gc,
                at,
                int(sum(seq_n)),
                n_cds,
                int(np.sum(cds_mask)),
            )
        return _seq_stats

    @classmethod
    def quick_parse(cls, seq_iter: Iterable[SeqRecord.SeqRecord], min_contig_len=0):
        """
        Quickly parses sequence records to obtain their lengths, ignoring all other statistics.
        
        Parameters:
        	seq_iter (Iterable[SeqRecord.SeqRecord]): An iterable of sequence records to process.
        	min_contig_len (int, optional): Minimum sequence length to include; shorter sequences are skipped.
        
        Returns:
        	_dict (dict[str, SeqStat]): Dictionary mapping sequence IDs to SeqStat instances containing only sequence length.
        """
        _seq_stats: dict[str, SeqStat] = {}
        for seq in seq_iter:
            if len(seq.seq) < min_contig_len:
                continue
            assert seq.id
            _seq_stats[seq.id] = cls(len(seq), 0, 0, 0, 0, 0, 0)

        return _seq_stats


class _BinStatisticContainer:
    @classmethod
    def read_gff(
        cls,
        filename: PathLike,
        refernce_file: PathLike | None = None,
        min_contig_len=0,
        min_aa_len=33,
    ):
        """
        Reads a GFF file and returns an instance initialized with parsed sequence statistics.
        
        If a reference file is provided, resets the parser's reference before parsing. Filters sequences and coding regions by minimum contig and coding sequence length.
        
        Parameters:
            filename: Path to the GFF file.
            refernce_file: Optional path to a reference file for resetting the parser.
            min_contig_len: Minimum contig length to include.
            min_aa_len: Minimum coding sequence length to include.
        
        Returns:
            An instance of the class initialized with parsed sequence statistics.
        """
        parser = Parse(filename)
        if refernce_file:
            parser = parser.reset_reference(refernce_file)
        return cls(parser(), filename, min_contig_len, min_aa_len=min_aa_len)

    @classmethod
    def read_gff_parser(cls, parser: Parse, min_contig_len=0, min_aa_len=33):
        """
        Create an instance from a GFF parser object, loading sequence statistics with optional minimum contig and CDS length filters.
        
        Parameters:
        	parser (Parse): A parser object yielding sequence records and providing a GFF file reference.
        	min_contig_len (int, optional): Minimum contig length to include. Defaults to 0.
        	min_aa_len (int, optional): Minimum coding sequence (CDS) length to include. Defaults to 33.
        
        Returns:
        	Instance of the class initialized with parsed sequence statistics.
        """
        return cls(parser(), parser.gff_file, min_contig_len, min_aa_len=min_aa_len)

    @classmethod
    def read_contig(cls, filename, format="fasta", min_contig_len=0):
        """
        Reads contig sequences from a file and initializes the container with sequence statistics.
        
        Parameters:
            filename (str or Path): Path to the contig file.
            format (str): File format for parsing sequences (default is "fasta").
            min_contig_len (int): Minimum contig length to include.
        
        Returns:
            An instance of the class initialized with parsed sequence statistics.
        """
        return cls(SeqIO.parse(filename, format), filename, min_contig_len)

    def __init__(
        self,
        seqiter: Iterable[SeqRecord.SeqRecord],
        source_file,
        min_contig_len=0,
        loader: Callable[
            [Iterable[SeqRecord.SeqRecord], int], dict[str, SeqStat]
        ] = SeqStat.parse,
        **parse_kwargs,
    ):
        """
        Initialize the container with sequence statistics parsed from an iterable of sequence records.
        
        Parameters:
            seqiter: An iterable of SeqRecord objects representing sequences to be analyzed.
            source_file: The source file or identifier associated with the sequence data.
            min_contig_len: Minimum contig length to include in statistics (default: 0).
            loader: Callable used to parse sequence statistics from the iterable (default: SeqStat.parse).
            **parse_kwargs: Additional keyword arguments passed to the loader function.
        """
        self._seq_stats = loader(seqiter, min_contig_len, **parse_kwargs)
        self.source_file = source_file
        self.min_contig_len = min_contig_len

    @classmethod
    def to_data_frame(cls, states: dict):
        """
        Convert a dictionary of namedtuple statistics to a pandas DataFrame with consistent column ordering.
        
        Parameters:
            states (dict): Dictionary mapping keys to namedtuple instances containing statistical fields.
        
        Returns:
            pd.DataFrame: DataFrame where each row corresponds to a key in `states` and columns represent namedtuple fields, ordered consistently.
        """
        all_fields = sorted({i._fields for i in states.values()})
        field_order = [i for j in all_fields for i in j]
        fields = sorted({i for j in all_fields for i in j}, key=field_order.index)

        bs = pd.DataFrame(
            {i: v._asdict() for i, v in states.items()},
        ).T
        return bs[fields]


class BinStatisticContainer(_BinStatisticContainer):
    @classmethod
    def quick_read_contig(cls, filename, format="fasta", min_contig_len=0):
        """
        Quickly loads sequence statistics from a contig file using minimal parsing.
        
        Parameters:
            filename (str or Path): Path to the contig file.
            format (str): File format for parsing sequences (default is "fasta").
            min_contig_len (int): Minimum contig length to include in statistics (default is 0).
        
        Returns:
            BinStatisticContainer: An instance containing basic sequence statistics for the contigs.
        """
        return cls(
            SeqIO.parse(filename, format),
            filename,
            min_contig_len,
            SeqStat.quick_parse,
        )

    def seq_stats(self, min_contig_len=0):
        """
        Yield sequence statistics for contigs with length greater than or equal to the specified minimum.
        
        Parameters:
            min_contig_len (int, optional): Minimum contig length to include. Defaults to 0.
        
        Returns:
            Generator[tuple[str, SeqStat]]: Pairs of sequence ID and corresponding SeqStat for qualifying contigs.
        """
        _min_contig_len = max(min_contig_len, self.min_contig_len)
        return (
            (seq_id, seq_stat)
            for seq_id, seq_stat in self._seq_stats.items()
            if seq_stat.len >= _min_contig_len
        )

    class BinStatistic(NamedTuple):
        gc: float
        gc_std: float
        bp_size: int
        max_contig_len: int
        contigs_num: int
        contig_n50: int
        ambiguous_bases_num: int
        contig_cutoff: int
        coding_density: float
        genes_num: int

    def statistic(self, min_contig_len: int | None = None):
        """
        Compute summary statistics for the genome bin, including GC content, sequence length metrics, ambiguous base count, coding density, and gene count.
        
        Parameters:
            min_contig_len (int, optional): Minimum contig length to include in the statistics. Defaults to the container's minimum if not specified.
        
        Returns:
            BinStatistic: Named tuple containing GC content, GC standard deviation, total base pairs, maximum contig length, number of contigs, N50, ambiguous base count, contig cutoff, coding density (ratio of coding length to non-ambiguous bases), and gene count.
        """
        if min_contig_len is None:
            min_contig_len = 0
        _min_contig_len = max(min_contig_len, self.min_contig_len)
        gc, gc_std = self.calculate_gc_std(_min_contig_len)
        gss = self.calculate_seq_stats(_min_contig_len)
        coding_len, num_orfs = self.calculate_prot_coding_length(_min_contig_len)

        return self.BinStatistic(
            gc=gc,
            gc_std=gc_std,
            bp_size=gss.sum,
            max_contig_len=gss.max,
            contigs_num=gss.num,
            contig_n50=gss.n50,
            ambiguous_bases_num=gss.numN,
            contig_cutoff=min_contig_len,
            coding_density=float(coding_len) / (gss.sum - gss.numN),
            genes_num=num_orfs,
        )

    def calculate_gc_std(
        self,
        min_contig_len=0,
        min_seq_len_gc_std=1000,
    ):
        """
        Calculate fraction of nucleotides that are G or C.
        """
        total_gc = 0
        total_at = 0
        gc_pct_per_seq = []
        for _, seq_stat in self.seq_stats(min_contig_len):
            total_gc += seq_stat.gc
            total_at += seq_stat.at

            if seq_stat.len > min_seq_len_gc_std:
                gc_pct_per_seq.append(seq_stat.gc_pct)

        if (total_gc + total_at) > 0:
            gc_pct = float(total_gc) / (total_gc + total_at)
        else:
            gc_pct = 0.0

        var_gc = 0
        if len(gc_pct_per_seq) > 1:
            var_gc = mean(list(map(lambda x: (x - gc_pct) ** 2, gc_pct_per_seq)))

        return gc_pct, math.sqrt(var_gc)

    class _GenomeSeqsStatistics(NamedTuple):
        sum: int
        max: int
        num: int
        n50: int
        numN: int

    def calculate_seq_stats(self, min_contig_len=0):
        """Calculate scaffold length statistics (min length, max length, total length, N50, # contigs)."""
        contig_lens = []
        numAmbiguousBases = 0
        for _, seq_stat in self.seq_stats(min_contig_len):
            contig_lens.append(seq_stat.len)

            numAmbiguousBases += seq_stat.n

        contig_N50 = calculateN50(contig_lens)

        return self._GenomeSeqsStatistics(
            sum(contig_lens),
            max(contig_lens),
            len(contig_lens),
            contig_N50,
            numAmbiguousBases,
        )

    def calculate_prot_coding_length(self, min_contig_len=0):
        """
        Calculate the total coding sequence length and number of coding sequences for the genome bin.
        
        Parameters:
            min_contig_len (int, optional): Minimum contig length to include in the calculation. Defaults to 0.
        
        Returns:
            tuple: A tuple (len_cds, n_cds) where len_cds is the total length of coding sequences and n_cds is the total number of coding sequences across all contigs meeting the length threshold.
        """
        len_aa = 0
        n_aa = 0
        for _, seq_stat in self.seq_stats(min_contig_len):
            len_aa += seq_stat.len_cds
            n_aa += seq_stat.n_cds

        return len_aa, n_aa

    def dump(self, filename: PathLike | None = None):
        """
        Serialize the current instance to a pickle file.
        
        If no filename is provided, the output file is named using the source file with a '-stat.pkl' suffix. The method ensures the output directory exists before writing.
        """
        self.parse = lambda: ()
        if filename is None:
            pickle_filename = Path(f"{self.source_file}-stat.pkl")
        else:
            pickle_filename = Path(filename)
        pickle_filename.parent.mkdir(parents=True, exist_ok=True)
        with open(pickle_filename, "wb") as po:
            pickle.dump(self, po)

    @classmethod
    def load(cls, filename: PathLike) -> "BinStatisticContainer":
        """
        Load a BinStatisticContainer instance from a pickle file.
        
        Parameters:
            filename: Path to the pickle file containing a serialized BinStatisticContainer.
        
        Returns:
            BinStatisticContainer: The deserialized instance loaded from the file.
        """
        with open(filename, "rb") as pi:
            return pickle.load(pi)
