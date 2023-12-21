# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-11 20:11:35
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-10-22 21:34:55
 * @FilePath: /genome/test/genome/test_prokka.py
 * @Description:
__file__ = "test/genome/test_prokka.py"
"""

import os
from pathlib import Path
from genome.prokka import prokka_gff_onethread, prokka_gff_multithread

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
