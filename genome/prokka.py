# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-11 13:49:35
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-10-12 16:09:38
 * @FilePath: /genome/genome/prokka.py
 * @Description:
"""

import os
import shutil
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import NamedTuple, Iterable, Literal, Union
import math
from numpy import mean

from BCBio import GFF
from Bio import SeqIO, SeqRecord, SeqFeature
from snakemake import main as smk


prokka_kindom = Literal["Archaea", "Bacteria", "Mitochondria", "Viruses"]
PathLike = Union[str, Path]


def prokka_gff_onethread(
    genome: Union[PathLike, Iterable[SeqRecord.SeqRecord]],
    kingdom: prokka_kindom = "Bacteria",
    gff_out: PathLike = "",
) -> Path:
    # infer gff_out automatically if not given in some cases
    if not gff_out:
        if not isinstance(genome, str) and not isinstance(genome, Path):
            raise ValueError("inital filename must provided")
        if not str(genome).endswith(".fa"):
            raise ValueError("genome file must endswith '.fa'")
        gff_out_ = Path(str(genome)[:-3] + f"-prokka.{kingdom}.gff")
    else:
        gff_out_ = Path(gff_out)

    with NamedTemporaryFile("w", suffix=".fa", delete=True) as tmpf:
        tpmf_out = Path(f"{tmpf.name[:-3]}-prokka.{kingdom}.gff")

        if not isinstance(genome, str) and not isinstance(genome, Path):
            SeqIO.write(genome, tmpf, "fasta")
        else:
            if not Path(genome).is_file():
                raise FileNotFoundError(f"file {genome} does not exist, please check.")
            with open(genome) as fi:
                while True:
                    # read 16 Kib words one time
                    block = fi.read(65536)
                    if not block:
                        break
                    tmpf.write(block)
                tmpf.flush()

        smk_workflow = Path(__file__).parent.parent / "workflow"
        target_smk_file = smk_workflow / "genome.smk"
        smk_params = f"-s {target_smk_file} " f"{tpmf_out} " f"--use-conda " f"-c1 -rp"

        try:
            os.system(f"ls {tmpf.name}")
            print("params:", "snakemake", smk_params)
            smk(smk_params)
        except SystemExit as se:
            if se.code:
                print(se.code, se.with_traceback(None))
                raise RuntimeError("snakemake seems not run successfully.")
            else:
                gff_out_.parent.mkdir(parents=True, exist_ok=True)
                shutil.move(tpmf_out, gff_out_)
                return gff_out_

    raise NotImplementedError("")


def prokka_gff_multithread(
    genomes: Iterable[PathLike],
    kingdom: prokka_kindom = "Bacteria",
    gff_out_dir: PathLike = "",
    threads: int = 8,
) -> Iterable[Path]:
    """WARNING: untested"""
    # if many genomes are provided, the file must exist
    _genome_files = [Path(file).expanduser().absolute() for file in genomes]
    for genome in _genome_files:
        if not str(genome).endswith(".fa"):
            raise ValueError("genome file must endswith '.fa'")

    genome_files = []
    if gff_out_dir:
        gff_out_dir_ = Path(gff_out_dir).expanduser().absolute()
        gff_out_dir_.mkdir(parents=True, exist_ok=True)

        _genome_files_dict = {file.name: file for file in _genome_files}
        if len(_genome_files_dict) != len(_genome_files):
            raise ValueError("cannot collect genome_files with same name")

        for genome_name, genome_path in _genome_files_dict.items():
            new_file = gff_out_dir_ / genome_name
            if new_file.exists():
                if new_file != genome_path:
                    raise FileExistsError(genome_path)
            else:
                shutil.copy(genome_path, new_file)
            genome_files.append(new_file)
    else:
        genome_files.extend(_genome_files)

    smk_workflow = Path(__file__).parent.parent / "workflow"
    target_smk_file = smk_workflow / "genome.smk"
    tpmf_outs = [f"{str(genome)[:-3]}-prokka.{kingdom}.gff" for genome in genomes]
    tpmf_outs_str = " ".join(tpmf_outs)
    smk_params = (
        f"-s {target_smk_file} " f"{tpmf_outs_str} " f"--use-conda " f"-c{threads} -rp"
    )
    try:
        print("params:", "snakemake", smk_params)
        smk(smk_params)
    except SystemExit as se:
        if se.code:
            print(se.code, se.with_traceback(None))
            raise RuntimeError("snakemake seems not run successfully.")

    # clean up temp added genomes
    if gff_out_dir:
        for new_file in genome_files:
            if new_file not in _genome_files:
                shutil.rmtree(new_file)

    return [Path(tpmf_out) for tpmf_out in tpmf_outs]


gff_out_format = Literal["faa", "fna"]


def prokka_gff_extract_protein_fa(
    gff_file: PathLike, out_format: gff_out_format = "faa"
):
    with open(gff_file) as in_handle:
        rec: SeqRecord.SeqRecord = None
        for rec in GFF.parse(in_handle):
            fet: SeqFeature.SeqFeature = None
            for fet in rec.features:
                if fet.type != "CDS":
                    continue
                seq: SeqRecord.SeqRecord = fet.extract(rec)
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


def calculate_gc(
    seqs: Iterable[SeqRecord.SeqRecord],
    seqStats: dict[str, dict[str, float]] = None,
    min_seq_len_gc_std=1000,
):
    """
    Calculate fraction of nucleotides that are G or C.
    modified from checkm
    """
    totalGC = 0
    totalAT = 0
    gcPerSeq = []
    for seq in seqs:
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


def calculateN50(seqLens):
    thresholdN50 = sum(seqLens) / 2.0

    seqLens.sort(reverse=True)

    testSum = 0
    for seqLen in seqLens:
        testSum += seqLen
        if testSum >= thresholdN50:
            N50 = seqLen
            break

    return N50


class GenomeSeqsStatistics(NamedTuple):
    sum: int
    max: int
    num: int
    n50: int
    numN: int


def calculate_seq_stats(contigs: Iterable[SeqRecord.SeqRecord], seqStats=None):
    """Calculate scaffold length statistics (min length, max length, total length, N50, # contigs)."""
    contig_lens = []
    numAmbiguousBases = 0
    for contig in contigs:
        contig_len = len(contig)
        contig_lens.append(contig_len)

        if seqStats:
            seqStats[contig.id]["Length"] = contig_len

        numAmbiguousBases += contig.seq.lower().count("n")

    contig_N50 = calculateN50(contig_lens)

    return GenomeSeqsStatistics(
        sum(contig_lens),
        max(contig_lens),
        len(contig_lens),
        contig_N50,
        numAmbiguousBases,
    )


class _ContigCodings(NamedTuple):
    contig_len: int
    coding_lens: list[int]


def calculate_prot_coding_length(contigs: Iterable[SeqRecord.SeqRecord]):
    """Calculate coding density of putative genome bin."""
    contigs_codings: dict[str, _ContigCodings] = {}
    rec: SeqRecord.SeqRecord
    for rec in contigs:
        cc = contigs_codings.setdefault(rec.id, _ContigCodings(len(rec.seq), []))

        fet: SeqFeature.SeqFeature
        for fet in rec.features:
            if fet.type != "CDS":
                continue
            cc.coding_lens.append(int(fet.location.end - fet.location.start))

    return sum(
        coding_len for cc in contigs_codings.values() for coding_len in cc.coding_lens
    ), sum(len(cc.coding_lens) for cc in contigs_codings.values())


class BinStatistic(NamedTuple):
    gc: float
    gc_std: float
    bp_size: int
    max_contig_len: int
    contigs_num: int
    contig_n50: int
    ambiguous_bases_num: int
    coding_density: float
    genes_size: int


def prokka_gff_bin_statistic(gff_file: PathLike):

    # read scaffolds
    contigs = [rec for rec in GFF.parse(gff_file)]

    # calculate GC statistics
    gc, gc_std = calculate_gc(contigs)

    # calculate statistics related to contigs and scaffolds
    gss = calculate_seq_stats(contigs)

    # calculate coding density statistics
    coding_len, num_orfs = calculate_prot_coding_length(contigs)

    bin_stats = BinStatistic(
        gc,
        gc_std,
        gss.sum,
        gss.max,
        gss.num,
        gss.n50,
        gss.numN,
        (float(coding_len) / gss.sum if coding_len > 0 and gss.sum else -1),
        num_orfs,
    )

    return bin_stats
