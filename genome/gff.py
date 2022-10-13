# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-12 19:32:50
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-10-13 14:14:27
 * @FilePath: /genome/genome/gff.py
 * @Description:
"""

import gzip
import math
from pathlib import Path
from pickle import dump, load
from typing import Any, Final, Generator, Literal, NamedTuple, TextIO, Union, Iterable

from BCBio import GFF
from Bio import SeqFeature, SeqRecord
from numpy import mean

PathLike = Union[str, Path]
gff_out_format = Literal["faa", "fna"]


def parse(
    gff_files: Union[PathLike, TextIO],
    base_dict: dict[str, SeqRecord.SeqRecord] = None,
    limit_info: dict[str, Any] = None,
    target_lines=None,
) -> Generator[SeqRecord.SeqRecord, None, None]:
    """
    High level interface to parse GFF files into SeqRecords and SeqFeatures.
    Add type hints to this function.
    Can automatically read gzip files.
    """
    if str(gff_files).endswith(".gz"):
        gff_files = gzip.open(gff_files, "r")  # type: ignore  # I'm sure this will return a TextIO
    return GFF.parse(gff_files, base_dict, limit_info, target_lines)


def write(
    recs: Iterable[SeqRecord.SeqRecord], out_handle: Union[PathLike, TextIO], include_fasta=False
):
    """
    High level interface to write GFF files into SeqRecords and SeqFeatures.
    Add type hints to this function.
    Can automatically write gzip files.
    """
    if str(out_handle).endswith(".gz"):
        out_handle = gzip.open(out_handle, "r")  # type: ignore  # I'm sure this will return a TextIO
    return GFF.write(recs, out_handle, include_fasta)


def gff_extract_protein_fa(
    gff_file: PathLike, out_format: gff_out_format = "faa", min_aa_length=33
):
    """
    min_aa_length acturally refer to at least 32 aa complete protein,
    for a complete terminal codon.
    """
    min_gene_length = min_aa_length * 3

    rec: SeqRecord.SeqRecord = None
    for rec in parse(gff_file):
        fet: SeqFeature.SeqFeature = None
        for fet in rec.features:
            if fet.type != "CDS":
                continue
            seq: SeqRecord.SeqRecord = fet.extract(rec)
            if len(seq) < min_gene_length:
                continue
            if out_format == "faa":
                seq = seq.translate()
            seq.id = rec.id + "_" + str(int(fet.id.rsplit("_", 1)[1]))
            seq.description = " # ".join(
                (
                    str(i)
                    for i in (
                        "",
                        fet.location.start,
                        fet.location.end,
                        fet.location.strand,
                        ";".join(
                            f"{k}={','.join(v)}" for k, v in fet.qualifiers.items()
                        ),
                    )
                )
            ).strip()
            yield seq


def calculateN50(seqLens: list[int]):
    thresholdN50 = sum(seqLens) / 2.0

    seqLens.sort(reverse=True)

    testSum = 0
    for seqLen in seqLens:
        testSum += seqLen
        if testSum >= thresholdN50:
            N50 = seqLen
            break

    return N50


class BinStatistic:
    """
    modified from checkm
    """

    class _SeqStat(NamedTuple):
        len: int = 0
        gc_pct: float = 0.0
        gc: int = 0
        at: int = 0
        n: int = 0

    def __init__(self, gff_file: PathLike):
        # def prokka_gff_bin_statistic(gff_file: PathLike):

        self.gff_file = Path(gff_file)
        self._seq_stat: dict[str, "BinStatistic._SeqStat"] = {}

        gc, gc_std = self.calculate_gc_std()
        gss = self.calculate_seq_stats()
        coding_len, num_orfs = self.calculate_prot_coding_length()

        self.gc: Final[float] = gc
        self.gc_std: Final[float] = gc_std
        self.bp_size: Final[int] = gss.sum
        self.max_contig_len: Final[int] = gss.max
        self.contigs_num: Final[int] = gss.num
        self.contig_n50: Final[int] = gss.n50
        self.ambiguous_bases_num: Final[int] = gss.numN
        self.coding_density: Final[float] = (
            float(coding_len) / gss.sum if coding_len > 0 and gss.sum else -1
        )
        self.genes_size: Final[int] = num_orfs

    def __repr__(self):
        return (
            f"Statistics of bin <{self.gff_file}>: \n"
            f"GC:              {self.gc:4f}\n"
            f"GC std:          {self.gc_std:4f}\n"
            f"Genome size:     {self.bp_size}\n"
            f"Longest contig:  {self.max_contig_len}\n"
            f"Contig number:   {self.contigs_num}\n"
            f"Contig N50:      {self.contig_n50}\n"
            f"Ambiguous bases: {self.ambiguous_bases_num}\n"
            f"Coding density:  {self.coding_density:4f}\n"
            f"predicted genes: {self.genes_size}"
        )

    def contigs(self) -> Generator[SeqRecord.SeqRecord, None, None]:
        # read contigs
        return (rec for rec in parse(self.gff_file))

    @property
    def seq_stats(self):
        if not self._seq_stat:
            for seq in self.contigs():
                a, c, g, t, u, n = (seq.seq.lower().count(base) for base in "acgtun")

                at = a + u + t
                gc = g + c

                if (gc + at) > 0:
                    gcContent = float(gc) / (gc + at)
                else:
                    gcContent = 0.0

            self._seq_stat[seq.id] = type(self)._SeqStat(len(seq), gcContent, gc, at, n)
        return self._seq_stat.items()

    def calculate_gc_std(
        self,
        min_seq_len_gc_std=1000,
    ):
        """
        Calculate fraction of nucleotides that are G or C.
        """
        total_gc = 0
        total_at = 0
        gc_pct_per_seq = []
        for _, seq_stat in self.seq_stats:
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

    def calculate_seq_stats(self):
        """Calculate scaffold length statistics (min length, max length, total length, N50, # contigs)."""
        contig_lens = []
        numAmbiguousBases = 0
        for _, seq_stat in self.seq_stats:
            contig_lens.append(seq_stat.len)

            numAmbiguousBases += seq_stat.n

        contig_N50 = calculateN50(contig_lens)

        return type(self)._GenomeSeqsStatistics(
            sum(contig_lens),
            max(contig_lens),
            len(contig_lens),
            contig_N50,
            numAmbiguousBases,
        )

    class _ContigCodings(NamedTuple):
        contig_len: int
        coding_lens: list[int]

    def calculate_prot_coding_length(self):
        """Calculate coding density of putative genome bin."""
        contigs_codings: dict[str, "BinStatistic._ContigCodings"] = {}
        rec: SeqRecord.SeqRecord
        for rec in self.contigs():
            cc = contigs_codings.setdefault(
                rec.id, type(self)._ContigCodings(len(rec.seq), [])
            )

            fet: SeqFeature.SeqFeature
            for fet in rec.features:
                if fet.type != "CDS":
                    continue
                cc.coding_lens.append(int(fet.location.end - fet.location.start))

        return (
            sum(
                coding_len
                for cc in contigs_codings.values()
                for coding_len in cc.coding_lens
            ),
            sum(len(cc.coding_lens) for cc in contigs_codings.values()),
        )

    def dump(self, filename: PathLike = None):
        if filename is None:
            pickle_filename = Path(f"{self.gff_file}-stat.pkl")
        else:
            pickle_filename = Path(filename)
        pickle_filename.parent.mkdir(parents=True, exist_ok=True)
        with open(pickle_filename, "wb") as po:
            dump(self, po)

    @classmethod
    def load(cls, filename: PathLike):
        with open(filename, "rb") as pi:
            return load(pi)
