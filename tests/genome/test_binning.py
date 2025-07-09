# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-25 20:53:06
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-07-02 23:27:05
 * @FilePath: /genome/tests/genome/test_binning.py
 * @Description:
__file__ = "test/genome/test_binning.py"
"""
# """

import os
import pandas as pd
import yaml

from tests import Path, temp_output, test_files, test_temp
from tests.genome._decorator import pytest_mark_resource

from genome.binning import (
    BinningConfig,
    BinningInput,
    bin_union,
    check_bams,
    default_bin_methods,
    rules_dir,
)


@temp_output
def test_binning_input_to_config(test_temp: Path):
    """
    Test that a BinningInput instance correctly generates a YAML configuration file with expected paths for contig, JGI, and BAM list files.
    """
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
    """
    Runs the bin union process using the 'unitem_unanimous' method and loads the resulting contig-to-bin assignments.
    
    This test executes the bin union workflow with specified binning methods and input files, then reads and groups the resulting contig-to-bin mapping to verify output structure.
    """
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


@temp_output
def test_smk_filter__rename_filtered_ls_tsv(test_temp: Path):
    """
    Test that the Snakemake filter workflow correctly renames filtered bin files and updates associated metadata.
    
    This test creates mock CheckM output and a list of bin FASTA files, runs the `filter.smk` Snakemake rule to rename bins, and verifies that both the TSV metadata and the list file reflect the new naming convention.
    """
    prefix = test_temp / "{any}-bins/filter/{method}{marker}-{check}_bins".format(
        any="any", method="method", marker="marker", check="check"
    )
    checkm_raw_str = (
        "Bin Id\tMarker lineage\tCompleteness\tContamination\tStrain heterogeneity\n"
        "bin_1\tBacteria\t99.99\t0.11\t0.33\n"
        "bin_2\tBacteria\t76.42\t5.23\t63.24\n"
        "bin_3\tBacteria\t99.99\t0.22\t0.44\n"
    )
    prefix = Path(prefix)
    prefix.mkdir(parents=True, exist_ok=True)
    with prefix.with_suffix(".tsv").open("w") as f:
        f.write(checkm_raw_str)
    with prefix.with_suffix(".ls").open("w") as f:
        for i in range(1, 4):
            bini = prefix / f"bin_{i}.fa"
            bini.touch()
            f.write(f"{bini}\n")
    o_prefix = Path(f"{prefix}-rename/Hoho_bins")
    assert not os.system(
        f"cd {test_temp} && "
        "python -m snakemake "
        f"-s {Path(rules_dir).resolve() / "binning" / "filter.smk"} "
        f"{o_prefix}.ls "
        "-p -c 1 "
        "> dryrun.log"
    )
    with o_prefix.with_suffix(".tsv").open() as f:
        assert f.read() == checkm_raw_str.replace("bin_", "Hoho-000")
    with o_prefix.with_suffix(".ls").open() as f:
        assert f.read() == "".join(
            f"{o_prefix}/Hoho-000{i}/Hoho-000{i}.fa\n" for i in range(1, 4)
        )
