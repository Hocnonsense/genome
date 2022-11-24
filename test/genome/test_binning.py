# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-25 20:53:06
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-11-24 17:05:24
 * @FilePath: /genome/test/genome/test_binning.py
 * @Description:
 * @Description:
__file__ = "/home/hwrn/software/genome/test/genome/test_binning.py"
"""

from pathlib import Path

import pandas as pd

from genome.binning import BinningConfig, bin_union

test_temp = Path(__file__).parent.parent / "temp"
test_files = Path(__file__).parent.parent / "file"


def test_binning_config_to_config():
    bc = BinningConfig(
        MIN_BIN_CONTIG_LEN=2500,
        contig=str(test_files / "02_assem..TY.041_cut.fa"),
        bams=[],
        bams_ls="",
        jgi=str(test_files / "02_assem..TY.041_cut-jgi.depth"),
        bin_single=str(test_files / "binsingle"),
        bin_union_dir=str(test_temp / "union"),
    )
    bc.to_config(test_temp / "test_binning.yaml")
    bc1 = BinningConfig()
    bc1.to_config(test_temp / "test_binning1.yaml")


def test_unitem_profile():
    binunion_tsv = bin_union(
        "unitem_unanimous",
        "",
        test_temp / "union",
        test_files / "02_assem..TY.041_cut.fa",
        test_files / "02_assem..TY.041_cut-jgi.depth",
        bin_single=test_files / "binsingle",
        threads=12,
    )
    pd.read_csv(binunion_tsv.ctg2mag, sep="\t", names=["Contig", "Bin"])


def test_dastool():
    binunion_tsv = bin_union(
        "dastool",
        "all",
        test_temp / "union",
        test_files / "02_assem..TY.041_cut.fa",
        test_files / "02_assem..TY.041_cut-jgi.depth",
        bin_single=test_files / "binsingle",
        threads=12,
    )
