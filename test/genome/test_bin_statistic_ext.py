# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-12 19:53:55
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-11-24 17:07:41
 * @FilePath: /genome/test/genome/test_bin_statistic_ext.py
 * @Description:
__file__ = "/home/hwrn/software/genome/test/genome/test_bin_statistic.py"
"""

from pathlib import Path

from genome.bin_statistic_ext import checkm

test_temp = Path(__file__).parent.parent / "temp"
test_files = Path(__file__).parent.parent / "file"


def test_checkm():
    test_fa = test_files / "02_assem..TY.041_cut.fa"
    test_c2b = test_temp / "union" / "dastool-all.tsv"
    test_out = test_temp / "union" / "dastool-all-checkm.tsv"
    (df := checkm(test_c2b, test_fa)).to_csv(test_out, sep="\t", index=False)
