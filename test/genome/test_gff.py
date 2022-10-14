# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-12 19:53:55
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-10-14 20:20:39
 * @FilePath: /genome/test/genome/test_gff.py
 * @Description:
__file__ = "/home/hwrn/software/genome/test/genome/test_prokka.py"
"""

import os
from pathlib import Path
from genome.gff import Parse, BinStatistic
from Bio import SeqIO

test_temp = Path(__file__).parent.parent / "temp"
test_files = Path(__file__).parent.parent / "file"


def test_gff_extract_protein_fa():
    gff = test_files / "metadecoder.1-prokka.Bacteria.gff"
    expect = test_files / "metadecoder.1.prokka.faa"

    test_out = test_temp / "metadecoder.1.prokka.faa"
    SeqIO.write(
        sorted(Parse(gff).extract(min_aa_length=33), key=lambda x: x.id),
        test_out,
        "fasta-2line",
    )

    print(f"expect file: {expect}", flush=True)
    os.system(f"md5sum {expect}")
    print(f"test output file: {test_out}", flush=True)
    os.system(f"md5sum {test_out}")


def test_prokka_gff_bin_statistic():
    gff = test_files / "metadecoder.1-prokka.Bacteria.gff"
    bs = BinStatistic(gff)
    assert int(bs.gc * 100) == 56
    assert int(bs.gc_std * 10000) == 206
    assert bs.bp_size == 1040626
    assert bs.max_contig_len == 9951
    assert bs.contigs_num == 293
    assert bs.contig_n50 == 3411
    assert bs.ambiguous_bases_num == 0
    assert int(bs.coding_density * 100) == 78
    assert bs.genes_size == 1130
