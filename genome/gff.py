# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-12 19:32:50
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-10-12 20:32:45
 * @FilePath: /genome/genome/gff.py
 * @Description:
"""

import math
from pathlib import Path
from typing import Any, Final, Generator, Literal, NamedTuple, TextIO, Union

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
    """
    return GFF.parse(gff_files, base_dict, limit_info, target_lines)


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

    def __init__(self, gff_file: PathLike):
        # def prokka_gff_bin_statistic(gff_file: PathLike):

        self.gff_file = Path(gff_file)

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

    def calculate_gc_std(
        self,
        seqStats: dict[str, dict[str, float]] = None,
        min_seq_len_gc_std=1000,
    ):
        """
        Calculate fraction of nucleotides that are G or C.
        """
        totalGC = 0
        totalAT = 0
        gcPerSeq = []
        for seq in self.contigs():
            a, c, g, t, u = (seq.seq.lower().count(base) for base in "acgtu")

            at = a + u + t
            gc = g + c

            totalGC += gc
            totalAT += at

            if (gc + at) > 0:
                gcContent = float(gc) / (gc + at)
            else:
                gcContent = 0.0

            if seqStats:
                seqStats[seq.id]["GC"] = gcContent

            if len(seq) > min_seq_len_gc_std:
                gcPerSeq.append(gcContent)

        if (totalGC + totalAT) > 0:
            GC = float(totalGC) / (totalGC + totalAT)
        else:
            GC = 0.0

        varGC = 0
        if len(gcPerSeq) > 1:
            varGC = mean(list(map(lambda x: (x - GC) ** 2, gcPerSeq)))

        return GC, math.sqrt(varGC)

    class _GenomeSeqsStatistics(NamedTuple):
        sum: int
        max: int
        num: int
        n50: int
        numN: int

    def calculate_seq_stats(self, seqStats=None):
        """Calculate scaffold length statistics (min length, max length, total length, N50, # contigs)."""
        contig_lens = []
        numAmbiguousBases = 0
        for contig in self.contigs():
            contig_len = len(contig)
            contig_lens.append(contig_len)

            if seqStats:
                seqStats[contig.id]["Length"] = contig_len

            numAmbiguousBases += contig.seq.lower().count("n")

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

        return sum(
            coding_len
            for cc in contigs_codings.values()
            for coding_len in cc.coding_lens
        ), sum(len(cc.coding_lens) for cc in contigs_codings.values())
