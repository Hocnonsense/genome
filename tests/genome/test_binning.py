# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-25 20:53:06
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-06-18 11:14:01
 * @FilePath: /genome/tests/genome/test_binning.py
 * @Description:
__file__ = "test/genome/test_binning.py"
"""
# """

import pandas as pd
import yaml

from genome.binning import (
    BinningConfig,
    BinningInput,
    bin_union,
    check_bams,
    default_bin_methods,
)

from tests import Path, temp_output, test_files, test_temp
from tests.genome._decorator import pytest_mark_resource


@temp_output
def test_binning_input_to_config(test_temp: Path):
    f = BinningInput(
        contig=str(test_files / "binny_contigs_4bins.fa"),
        lsbams=check_bams(test_temp / "binny_test_binning", []),
        jgi=str(test_files / "binny_reads_4bins-jgi.tsv"),
    ).dump_prefix(test_temp / "binny_test_binning", relpath=False)

    assert f == Path(test_temp / "binny_test_binning-bins.yaml")
    with open(f) as yi:
        assert yaml.safe_load(yi) == {
            "contig": str(test_files / "binny_contigs_4bins.fa"),
            "jgi": str(test_files / "binny_reads_4bins-jgi.tsv"),
            "lsbams": str(test_temp / "binny_test_binning-bams.list"),
        }


@temp_output
def test_binning_config_to_config(test_temp: Path):
    f1 = BinningConfig(
        MIN_BIN_CONTIG_LEN=2500, bin_config=test_temp / "binny_test_binning-bins.yaml"
    ).to_config(test_temp / "test_binning_config.yaml")
    with open(f1) as yi:
        assert yaml.safe_load(yi) == {
            "MIN_BIN_CONTIG_LEN": 2500,
            "bin_methods": default_bin_methods,
        }


@pytest_mark_resource
@temp_output
def test_unitem_profile(test_temp: Path):
    binunion_tsv = bin_union(
        method="unitem_unanimous",
        marker="",
        bin_prefix=test_temp / "binny_test_binning",
        contig=test_files / "binny_contigs_4bins.fa",
        jgi=test_files / "binny_reads_4bins-jgi.tsv",
        bams=[
            test_files / "binny_reads_4bins.bam",
        ],
        bin_methods=["metabat2_60_60", "metabat2_75_75", "concoct"],
        threads=12,
    )
    len(
        pd.read_csv(binunion_tsv.ctg2mag, sep="\t", names=["Contig", "Bin"])
        .groupby("Bin")["Contig"]
        .apply(set)
        .to_dict()
    )
    test_files / "binny_unitem_unanimous.tsv"
