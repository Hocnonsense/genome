# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-11 20:11:35
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-06-18 11:14:52
 * @FilePath: /genome/tests/genome/test_prokka.py
 * @Description:
__file__ = "test/genome/test_prokka.py"
"""

import os

from Bio import SeqIO
from genome.bin_statistic_ext import format_bin_input
from genome.prokka import prokka_gff_onethread, prokka_gff_multithread

from tests import Path, temp_output, test_files, test_temp
from tests.genome._decorator import pytest_mark_resource


@pytest_mark_resource
def test_prokka_gff_onethread():
    genome = test_files / "binny_contigs_4bins.fa"
    SeqIO.write(
        (i for i, _ in zip(SeqIO.parse(genome, "fasta"), range(10))),
        test_out_fa := test_temp / "binny_contigs_4bins-top10.fa",
        "fasta",
    )
    expect = test_files / "binny_contigs_4bins-prokka_Bacteria.gff"

    test_out = test_temp / "binny_contigs_4bins.gff"
    gff = prokka_gff_onethread(genome, "Bacteria", test_out)
    assert gff == test_out.expanduser().absolute()

    print(f"expect file: {expect}", flush=True)
    os.system(f"md5sum {expect}")
    print(f"test output file: {test_out}", flush=True)
    os.system(f"md5sum {test_out}")


@pytest_mark_resource
def test_prokka_gff_multithread():
    ctg2mag = test_files / "binny_unitem_unanimous.tsv"
    support = test_files / "binny_contigs_4bins.fa"
    bin_input_dir, binids, suffix = format_bin_input(
        bin_output=test_temp / "prodigal",
        bin_input=ctg2mag,
        support=support,
        keep_if_avail=False,
    )
    test_outs = [bin_input_dir / (binid + "-prokka_Bacteria.gff") for binid in binids]
    gffs = prokka_gff_multithread(
        (bin_input_dir / (binid + suffix) for binid in binids), "Bacteria", test_files
    )
    assert all(
        gff == test_out.expanduser().absolute()
        for gff, test_out in zip(gffs, test_outs)
    )
