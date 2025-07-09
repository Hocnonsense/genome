# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-12 19:53:55
* @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
* @LastEditTime: 2025-07-09 19:31:54
* @FilePath: /genome/tests/genome/test_bin_statistic_ext.py
 * @Description:
__file__ = "test/genome/test_bin_statistic.py"
"""

import pytest

from genome.bin_statistic_ext import checkm, format_bin_input

from tests import Path, temp_output, test_files, test_temp
from tests.genome._decorator import MARK_LIMIT_RESOURCE


@temp_output
def test_format_bin_input(test_temp: Path):
    bin_input = test_files / "binny_unitem_unanimous.tsv"
    support = test_files / "binny_contigs_4bins.fa"
    assert format_bin_input(
        bin_output=f"{test_temp}/bin_input",
        bin_input=bin_input,
        support=support,
    )[1:] == (
        [
            "unitem_unanimous-1",
            "unitem_unanimous-2",
            "unitem_unanimous-3",
            "unitem_unanimous-4",
            "unitem_unanimous-5",
            "unitem_unanimous-6",
        ],
        ".fa",
    )


@pytest.mark.skipif(MARK_LIMIT_RESOURCE, reason="no enough memory; improve required")
@temp_output
def test_checkm(test_temp: Path):
    bin_input = test_files / "binny_unitem_unanimous.tsv"
    support = test_files / "binny_contigs_4bins.fa"
    test_out = test_temp / "union" / "binny_unitem_unanimous-checkm.tsv"
    # warning: now we don't support checkm to call prodigal itself
    (df := checkm(bin_input, support, test_temp / "checkm")).to_csv(
        test_out, sep="\t", index=False
    )
