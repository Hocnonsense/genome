# -*- coding: utf-8 -*-
"""
 * @Date: 2022-11-24 16:23:50
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-11-24 16:28:59
 * @FilePath: /genome/genome/bin_statistic_ext.py
 * @Description:
"""
import os
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Optional, Union

import pandas as pd

from .bin_statistic import contig2bin

PathLike = Union[str, Path]


@dataclass
class CheckMFakeOptions:
    from checkm.defaultValues import DefaultValues

    subparser_name = "lineage_wf"
    output_dir: str = "./"
    bin_input: str = ""
    extension: str = "fa"
    file: str = "stdout"
    threads: int = 8

    bReducedTree: bool = False
    bKeepAlignment: bool = False
    bNucORFs: bool = False
    bCalledGenes: bool = False
    unique: int = 10
    multi: int = 10
    bForceDomain: int = False
    bNoLineageSpecificRefinement: int = False
    bIndividualMarkers: int = False
    bSkipAdjCorrection: int = False
    bSkipPseudoGeneCorrection: int = False
    aai_strain: float = 0.9
    alignment_file: Optional[str] = None
    bIgnoreThresholds: bool = False
    e_value: float = DefaultValues.E_VAL
    length: float = DefaultValues.LENGTH
    bQuiet: bool = False
    bTabTable = True
    pplacer_threads = 1

    def run(self):
        import checkm

        versionFile = open(os.path.join(checkm.__path__[0], "VERSION"))

        from checkm.logger import logger_setup
        from checkm.main import OptionsParser

        logger_setup(
            self.output_dir,
            "checkm.log",
            "CheckM",
            versionFile.readline().strip(),
            self.bQuiet,
        )

        checkmParser = OptionsParser()
        checkmParser.parseOptions(self)
        return self.output_dir


def checkm(
    bin_input: PathLike,
    support: Union[PathLike, str],
    output_dir=None,
    threads=10,
    **checkm_options,
):
    bin_input = Path(bin_input)
    with TemporaryDirectory() as _td:
        file = f"{_td}/checkm.tsv"
        if bin_input.is_dir():
            assert list(bin_input.glob(f"*{support}")), "input is not a valid bin path"
            support = str(support)
            bin_input_ = bin_input
        else:
            assert bin_input.is_file() and Path(support).is_file()
            bin_input_ = contig2bin(f"{_td}/bin_input", bin_input, support)
            support = "fa"
        if output_dir is None:
            output_dir = f"{_td}/checkm"
        CheckMFakeOptions(
            file=file,
            bin_input=str(bin_input_),
            output_dir=output_dir,
            extension=support,
            threads=threads,
            **checkm_options,  # type: ignore  # confirmed by users
        ).run()
        return pd.read_csv(file, sep="\t")


def gunc(
    bin_input: PathLike,
    support: Union[PathLike, str],
    output_dir=None,
    threads=10,
):
    pass
