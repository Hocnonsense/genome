# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-25 20:53:06
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-11-15 10:14:01
 * @FilePath: /genome/test/genome/test_binning.py
 * @Description:
 * @Description:
__file__ = "/home/hwrn/software/genome/test/genome/test_binning.py"
"""

from pathlib import Path

import pandas as pd

from genome.binning import BinningConfig, bin_union, checkm, contig2bin

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


def test_unitem_profile():
    binunion_tsv = bin_union(
        "unitem_consensus",
        "",
        test_temp / "union",
        test_files / "02_assem..TY.041_cut.fa",
        test_files / "02_assem..TY.041_cut-jgi.depth",
        bin_single=test_files / "binsingle",
        threads=12,
    )


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


def test_contig2bin():
    test_fa = test_files / "02_assem..TY.041_cut.fa"
    test_c2b = test_temp / "union" / "dastool-all.tsv"
    temp_bin = test_temp / "dastool-all"

    contig2bin(temp_bin, test_c2b, test_fa)


def test_checkm():
    test_fa = test_files / "02_assem..TY.041_cut.fa"
    test_c2b = test_temp / "union" / "dastool-all.tsv"
    test_out = test_temp / "union" / "dastool-all-checkm.tsv"
    (df := checkm(test_c2b, test_fa)).to_csv(test_out, sep="\t", index=False)
