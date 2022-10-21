# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-12 19:53:55
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-10-21 16:34:54
 * @FilePath: /genome/test/genome/test_gff.py
 * @Description:
__file__ = "/home/hwrn/software/genome/test/genome/test_prokka.py"
"""

import os
from pathlib import Path
from genome.gff import Parse
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


def test_gff_reset_reference():
    gff = test_files / "metadecoder.1-prokka.Bacteria.gff"
    genome = test_files / "metadecoder.1.fa"
    expect = test_files / "metadecoder.1.prokka.faa"

    test_out = test_temp / "metadecoder.1.prokka.faa"
    SeqIO.write(
        sorted(Parse(gff).reset_reference(genome).extract(min_aa_length=33), key=lambda x: x.id),
        test_out,
        "fasta-2line",
    )

    print(f"expect file: {expect}", flush=True)
    os.system(f"md5sum {expect}")
    print(f"test output file: {test_out}", flush=True)
    os.system(f"md5sum {test_out}")
