# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-12 19:53:55
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-12-21 22:07:56
 * @FilePath: /genome/tests/genome/test_bin_statistic_ext.py
 * @Description:
__file__ = "test/genome/test_bin_statistic.py"
"""

from pathlib import Path

from genome.bin_statistic_ext import checkm

try:
    from _decorator import temp_output, test_temp, test_files
except (ModuleNotFoundError, ImportError):
    from tests.genome._decorator import temp_output, test_temp, test_files


@temp_output
def test_checkm(test_temp: Path):
    bin_input = test_files / "binny_unitem_unanimous.tsv"
    support = test_files / "binny_contigs_4bins.fa"
    test_out = test_temp / "union" / "binny_unitem_unanimous-checkm.tsv"
    (df := checkm(bin_input, support, test_temp / "checkm")).to_csv(
        test_out, sep="\t", index=False
    )
