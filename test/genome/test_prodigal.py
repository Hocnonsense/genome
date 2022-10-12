# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-11 20:11:35
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-10-12 19:26:14
 * @FilePath: /genome/test/genome/test_prodigal.py
 * @Description:
__file__ = "/home/hwrn/software/genome/test/genome/test_prodigal.py"
"""

import os
from pathlib import Path
from genome.prodigal import prodigal_gff_onethread

test_temp = Path(__file__).parent.parent / "temp"
test_files = Path(__file__).parent.parent / "file"


def test_prokka_gff_onethread():
    genome = test_files / "metadecoder.4.fa"
    expect = test_files / "metadecoder.4-prodigal.single.gff"

    test_out = test_temp / "metadecoder.4.gff"
    gff = prodigal_gff_onethread(genome, "single", test_out)
    assert gff == test_out.expanduser().absolute()

    print(f"expect file: {expect}", flush=True)
    os.system(f"md5sum {expect}")
    print(f"test output file: {test_out}", flush=True)
    os.system(f"md5sum {test_out}")
