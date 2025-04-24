# -*- coding: utf-8 -*-
"""
* @Date: 2022-10-15 17:05:11
* @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
* @LastEditTime: 2025-04-24 13:24:39
* @FilePath: /genome/genome/bin_statistic.py
* @Description:
"""


from functools import cached_property
import math
import pickle
from pathlib import Path
import shutil
from typing import Callable, Iterable, NamedTuple, Union, Final, overload

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
        self.contig2bin_tsv = self.parse_contig2bin(contig2bin_tsv)
        self._contigs = self.parse_contigs(contigs)
        self._outdir = None
        if outdir:
            self.outdir = outdir

    @classmethod
    def parse_contig2bin(cls, contig2bin_tsv: PathLike | pd.DataFrame | pd.Series):
        """
        also support vamb-format table (format: contig\tbin):
        ```
        k141_135186 flag=1 multi=23.0000 len=8785\tvamb-1012                                                                                                         metabat2_90_90  metadecoder  vamb
        k141_420558 flag=0 multi=35.2019 len=10561\tvamb-1012
        k141_42780 flag=0 multi=26.0231 len=10767\tvamb-1012
        k141_746716 flag=1 multi=21.0002 len=9798\tvamb-1012
        k141_91422 flag=1 multi=50.0000 len=2703\tvamb-1012
        ```
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
        bin2seqs: dict[str, dict[str, SeqRecord.SeqRecord]] = {
            b: {} for b in self.contig2bin_tsv["bin"].unique()
        }
        for i in self.contigs:
            if i.name in self.contig2bin_tsv.index:
                bin2seqs[self.contig2bin_tsv.loc[i.name, "bin"]][i.name] = i
        return bin2seqs

    def parse_contigs(self, contigs: PathLike | Iterable[SeqRecord.SeqRecord]):
        if isinstance(contigs, str) or isinstance(contigs, Path):
            try:
                return lambda: SeqIO.parse(contigs, "fasta")
            except ValueError as e:
                return Parse(contigs)
        elif isinstance(contigs, Iterable):
            return lambda: list(contigs)
        else:
            raise ValueError("contigs should be a path or a iterable of SeqRecord")

    @property
    def contigs(self) -> Iterable[SeqRecord.SeqRecord]:
        return self._contigs()

    @property
    def outdir(self):
        if self._outdir is None:
            raise ValueError(".outdir not set")
        return self._outdir

    @outdir.setter
    def outdir(self, outdir: Path):
        td = Path(outdir)
        td.mkdir(parents=True, exist_ok=True)
        self._outdir = outdir
        return td

    def __call__(self, outdir: PathLike):
        self.outdir = Path(outdir)
        bin2seqs = self.bin2seqs
        for b, seqs in bin2seqs.items():
            SeqIO.write(seqs.values(), self.outdir / f"{b}.fa", "fasta-2line")
        return self.output

    @property
    def output(self):
        return Binput(self.outdir, list(self.bin2seqs), ".fa")

    @overload
    def extract(self, bin_name: str | int) -> Path: ...

    @overload
    def extract(self, bin_name: str | int, *bin_names: str | int) -> Iterable[Path]: ...

    def extract(self, bin_name: str | int, *bin_names: str | int):
        if isinstance(bin_name, int):
            b = self.contig2bin_tsv["bin"].unique()[bin_name]
        else:
            b = bin_name
        bout = self.outdir / f"{b}.fa"
        SeqIO.write(self.bin2seqs[b].values(), bout, "fasta-2line")
        if not bin_names:
            return bout
        else:
            return (*(bout,), *(self.extract(b) for b in bin_names))


class Binput(NamedTuple):
    bin_input: Path
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
        if {param kept_if_avail}:
            Only if bin_input is a dir and required genomes endswith "fa",
            bin_output will be kept as bin_input.

        This is designed to support snakemake to handle and keep intermediate files
        """
        bin_input = Path(bin_input)
        (bin_output_ := Path(bin_output)).mkdir(parents=True, exist_ok=True)
        if bin_input.is_dir():
            assert list(bin_input.glob(f"*{support}")), "input is not a valid bin path"
            if str(support).endswith(".fa") and keep_if_avail:
                binids: list[str] = [str(i)[:-3] for i in bin_input.glob(f"*{support}")]
                suffix = str(support)
                bin_input_dir = bin_input
            else:
                suffix = ".fa"
                bin_input_dir = bin_output_
                support_str_len = len(str(support))
                binids = []
                for bin_file in bin_input.glob(f"*{support}"):
                    binids.append(
                        bin_name := bin_file.name[:-support_str_len].rstrip(".")
                    )
                    shutil.copy(bin_file, bin_input_dir / f"{bin_name}.fa")
            return cls(bin_input_dir, binids, suffix)
        else:
            assert bin_input.is_file() and Path(support).is_file()
            return Contig2Bin(bin_input, support)(bin_output_)

    def fas(self):
        return (self.bin_input / f"{i}{self.suffix}" for i in self.binids)

    def fas_with(self, suffix: str):
        return Binput(self.bin_input, self.binids, suffix).fas()


def contig2bin(outdir: PathLike, contig2bin_tsv: PathLike, contigs: PathLike):
    return Contig2Bin(contig2bin_tsv, contigs)(outdir)


def calculateN50(seqLens: list[int]):
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
        _seq_stats: dict[str, SeqStat] = {}
        min_gene_len = int(min_aa_len) * 3
        for seq in seq_iter:
            if len(seq.seq) < min_contig_len:
                continue
            upper_seq = seq.seq.upper()
            a, c, g, t, u, n = (upper_seq.count(base) for base in "ACGTUN")

            at = a + u + t
            gc = g + c

            gcContent = 0.0 if (gcat := gc + at) <= 0 else float(gc) / gcat

            cds_mask = np.zeros(len(seq.seq))
            n_cds = 0
            for fet in seq.features:
                if fet.type == "CDS" and fet.location is not None:
                    if len(fet) < min_gene_len:
                        continue
                    cds_mask[fet.location.start : fet.location.end] = 1
                    n_cds += 1
            assert seq.id
            _seq_stats[seq.id] = cls(
                len(seq), gcContent, gc, at, n, n_cds, int(np.sum(cds_mask))
            )
        return _seq_stats

    @classmethod
    def quick_parse(cls, seq_iter: Iterable[SeqRecord.SeqRecord], min_contig_len=0):
        """
        warning: Every statement except contig length are disabled
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
        parser = Parse(filename)
        if refernce_file:
            parser = parser.reset_reference(refernce_file)
        return cls(parser(), filename, min_contig_len, min_aa_len=min_aa_len)

    @classmethod
    def read_gff_parser(cls, parser: Parse, min_contig_len=0, min_aa_len=33):
        return cls(parser(), parser.gff_file, min_contig_len, min_aa_len=min_aa_len)

    @classmethod
    def read_contig(cls, filename, format="fasta", min_contig_len=0):
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
        self._seq_stats = loader(seqiter, min_contig_len, **parse_kwargs)
        self.source_file = source_file
        self.min_contig_len = min_contig_len

    @classmethod
    def to_data_frame(cls, states: dict):
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
        return cls(
            SeqIO.parse(filename, format),
            filename,
            min_contig_len,
            SeqStat.quick_parse,
        )

    def seq_stats(self, min_contig_len=0):
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
            coding_density=float(coding_len) / gss.sum,
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
        """Calculate coding density of putative genome bin."""
        len_aa = 0
        n_aa = 0
        for _, seq_stat in self.seq_stats(min_contig_len):
            len_aa += seq_stat.len_cds
            n_aa += seq_stat.n_cds

        return len_aa, n_aa

    def dump(self, filename: PathLike | None = None):
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
        with open(filename, "rb") as pi:
            return pickle.load(pi)
