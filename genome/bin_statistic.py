# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-15 17:05:11
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-05-30 10:31:02
 * @FilePath: /genome/genome/bin_statistic.py
 * @Description:
"""


import math
from pathlib import Path
from pickle import dump, load
from typing import Iterable, NamedTuple, Union

import pandas as pd
from Bio import SeqIO, SeqRecord
from numpy import mean

from .gff import Parse

PathLike = Union[str, Path]


def contig2bin(outdir: PathLike, contig2bin_tsv: PathLike, contigs: PathLike):
    """
    will read a table such without head such as:
        (format: contig\tbin)
        ```
        k141_128353\tconcoct_0
        k141_15265\tconcoct_0
        k141_164694\tconcoct_0
        k141_172306\tconcoct_0
        k141_18944\tconcoct_0
        ```

    fix: will also support vamb-format table:
        (format: contig\tbin)
        ```
        k141_135186 flag=1 multi=23.0000 len=8785\tvamb-1012                                                                                                         metabat2_90_90  metadecoder  vamb
        k141_420558 flag=0 multi=35.2019 len=10561\tvamb-1012
        k141_42780 flag=0 multi=26.0231 len=10767\tvamb-1012
        k141_746716 flag=1 multi=21.0002 len=9798\tvamb-1012
        k141_91422 flag=1 multi=50.0000 len=2703\tvamb-1012
        ```
    """
    contig2bin_ = (
        pd.read_csv(contig2bin_tsv, sep="\t", names=["contig", "bin"])
        .assign(contig_id=lambda df: df["contig"].apply(lambda x: x.strip().split()[0]))
        .set_index("contig_id")[["bin"]]
    )

    td = Path(outdir)
    td.mkdir(parents=True, exist_ok=True)

    bin2seqs: dict[str, set[str]] = {b: set() for b in contig2bin_["bin"].unique()}
    for i in Parse(contigs)():
        if i.name in contig2bin_.index:
            bin2seqs[contig2bin_.loc[i.name, "bin"]].add(i.format("fasta-2line"))
    for b, seqs in bin2seqs.items():
        with open(td / f"{b}.fa", "w") as po:
            for seq in seqs:
                print(seq, file=po)

    return td, list(bin2seqs)


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
    n_aa: int = 0
    len_aa: int = 0

    @classmethod
    def load_seq_stats(cls, seq_iter: Iterable[SeqRecord.SeqRecord], min_contig_len=0):
        _seq_stats: dict[str, SeqStat] = {}
        for seq in seq_iter:
            if len(seq.seq) < min_contig_len:
                continue
            upper_seq = seq.seq.upper()
            a, c, g, t, u, n = (upper_seq.count(base) for base in "ACGTUN")

            at = a + u + t
            gc = g + c

            gcContent = 0.0 if (gcat := gc + at) <= 0 else float(gc) / gcat

            cc: list[int] = [
                # fet: SeqFeature.SeqFeature
                int(fet.location.end - fet.location.start)
                for fet in seq.features
                if fet.type == "CDS" and fet.location
            ]

            assert seq.id
            _seq_stats[seq.id] = cls(len(seq), gcContent, gc, at, n, len(cc), sum(cc))
        return _seq_stats

    @classmethod
    def quick_load_seq_stats(
        cls, seq_iter: Iterable[SeqRecord.SeqRecord], min_contig_len=0
    ):
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


class BinStatisticContainer:
    """
    modified from checkm
    """

    @classmethod
    def read_gff(
        cls,
        filename: PathLike,
        refernce_file: PathLike | None = None,
        min_contig_len=0,
    ):
        parser = Parse(filename)
        if refernce_file:
            parser = parser.reset_reference(refernce_file)

        return cls(SeqStat.load_seq_stats(parser()), filename, min_contig_len)

    @classmethod
    def read_gff_parser(
        cls,
        parser: Parse,
        min_contig_len=0,
    ):
        return cls(SeqStat.load_seq_stats(parser()), parser.gff_file, min_contig_len)

    @classmethod
    def read_contig(cls, filename, format="fasta", min_contig_len=0):
        return cls(
            SeqStat.load_seq_stats(SeqIO.parse(filename, format)),
            filename,
            min_contig_len,
        )

    @classmethod
    def quick_read_contig(cls, filename, format="fasta", min_contig_len=0):
        return cls(
            SeqStat.quick_load_seq_stats(SeqIO.parse(filename, format)),
            filename,
            min_contig_len,
        )

    @classmethod
    def read_seqiter(
        cls, seqiter: Iterable[SeqRecord.SeqRecord], filename, min_contig_len=0
    ):
        return cls(SeqStat.load_seq_stats(seqiter), filename, min_contig_len)

    def seq_stats(self, min_contig_len=0):
        _min_contig_len = max(min_contig_len, self.min_contig_len)
        return (
            (seq_id, seq_stat)
            for seq_id, seq_stat in self._seq_stats.items()
            if seq_stat.len >= _min_contig_len
        )

    def __init__(
        self,
        seq_stats: dict[str, SeqStat],
        source_file,
        min_contig_len=0,
    ):
        self._seq_stats = seq_stats
        self.source_file = source_file
        self.min_contig_len = min_contig_len

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
            len_aa += seq_stat.len_aa
            n_aa += seq_stat.n_aa

        return len_aa, n_aa

    def dump(self, filename: PathLike | None = None):
        self.parse = lambda: ()
        if filename is None:
            pickle_filename = Path(f"{self.source_file}-stat.pkl")
        else:
            pickle_filename = Path(filename)
        pickle_filename.parent.mkdir(parents=True, exist_ok=True)
        with open(pickle_filename, "wb") as po:
            dump(self, po)

    @classmethod
    def load(cls, filename: PathLike) -> "BinStatisticContainer":
        with open(filename, "rb") as pi:
            return load(pi)
