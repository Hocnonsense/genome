# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-11 20:11:35
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-10-12 19:00:46
 * @FilePath: /genome/test/genome/test_prokka.py
 * @Description:
__file__ = "/home/hwrn/software/genome/test/genome/test_prokka.py"
"""

import os
from pathlib import Path
from genome.prokka import (
    prokka_gff_onethread,
    prokka_gff_multithread,
    prokka_gff_extract_protein_fa,
    prokka_gff_bin_statistic,
)
from Bio import SeqIO

test_temp = Path(__file__).parent.parent / "temp"
test_files = Path(__file__).parent.parent / "file"


def test_prokka_gff_onethread():
    genome = test_files / "metadecoder.1.fa"
    expect = test_files / "metadecoder.1-prokka.Bacteria.gff"

    test_out = test_temp / "metadecoder.1.gff"
    gff = prokka_gff_onethread(genome, "Bacteria", test_out)
    assert gff == test_out.expanduser().absolute()

    print(f"expect file: {expect}", flush=True)
    os.system(f"md5sum {expect}")
    print(f"test output file: {test_out}", flush=True)
    os.system(f"md5sum {test_out}")


def test_prokka_gff_multithread():
    genomes = [test_files / f"metadecoder.{i}.fa" for i in (2, 3)]
    expects = [test_files / f"metadecoder.{i}-prokka.Bacteria.gff" for i in (2, 3)]

    test_outs = [test_temp / f"metadecoder.{i}-prokka.Bacteria.gff" for i in (2, 3)]
    gffs = prokka_gff_multithread(genomes, "Bacteria", test_files)
    assert all(
        gff == test_out.expanduser().absolute()
        for gff, test_out in zip(gffs, test_outs)
    )


def test_prokka_gff_to_faa():
    gff = test_files / "metadecoder.1-prokka.Bacteria.gff"
    expect = test_files / "metadecoder.1.prokka.faa"

    test_out = test_temp / "metadecoder.1.prokka.faa"
    SeqIO.write(prokka_gff_extract_protein_fa(gff), test_out, "fasta-2line")

    print(f"expect file: {expect}", flush=True)
    os.system(f"md5sum {expect}")
    print(f"test output file: {test_out}", flush=True)
    os.system(f"md5sum {test_out}")


def test_prokka_gff_bin_statistic():
    gff = test_files / "metadecoder.1-prokka.Bacteria.gff"
    bs = prokka_gff_bin_statistic(gff)
    assert int(bs.gc * 100) == 56
    assert int(bs.gc_std * 10000) == 206
    assert bs.bp_size == 1040626
    assert bs.max_contig_len == 9951
    assert bs.contigs_num == 293
    assert bs.contig_n50 == 3411
    assert bs.ambiguous_bases_num == 0
    assert int(bs.coding_density * 100) == 0.79
    assert bs.genes_size == 1130
