# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-11 20:11:35
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-10-11 22:56:08
 * @FilePath: /genome/test/genome/test_prokka.py
 * @Description:
__file__ = "/home/hwrn/software/genome/test/genome/test_prokka.py"
"""

import os
from pathlib import Path
from genome.prokka import prokka_gff_onethread, prokka_gff_extract_protein_fa
from Bio import SeqIO

test_temp = Path(__file__).parent.parent/"temp"
test_files = Path(__file__).parent.parent/"file"


def prokka_gff_onethread():
    genome = test_files / "metadecoder.1.fa"
    expect = test_files / "metadecoder.1-prokka.Bacteria.gff"

    test_out = test_temp / "metadecoder.1.gff"
    gff = prokka_gff_onethread(
        genome, "Bacteria", test_out
    )
    assert gff == test_out.expanduser().absolute()

    print(f"expect file: {expect}", flush=True)
    os.system(f"md5sum {expect}")
    print(f"test output file: {test_out}", flush=True)
    os.system(f"md5sum {test_out}")


def test_prokka_gff_to_faa():
    gff = test_files / "metadecoder.1-prokka.Bacteria.gff"
    expect = test_files / "metadecoder.1.prokka.faa"

    test_out = test_temp / "metadecoder.1.prokka.faa"
    SeqIO.write(prokka_gff_extract_protein_fa(gff), test_out, "fasta-2line")

    print(f"expect file: {expect}", flush=True)
    os.system(f"md5sum {expect}")
    print(f"test output file: {test_out}", flush=True)
    os.system(f"md5sum {test_out}")
