# -*- coding: utf-8 -*-
"""
 * @Date: 2022-11-24 16:23:50
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-11-26 14:23:02
 * @FilePath: /genome/genome/bin_statistic_ext.py
 * @Description:
"""
import os
import shutil
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Optional, Union

from snakemake import main as smk
import pandas as pd

from .bin_statistic import contig2bin
from .prodigal import prodigal_multithread

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


def format_bin_input(
    bin_output: PathLike,
    bin_input: PathLike,
    support: Union[PathLike, str],
    keep_if_avail=True,
):
    """
    if {param kept_if_avail}:
        Only if bin_input is a dir and required genomes endswith "fa",
        bin_output will be kept as bin_input.

    This is designed to support snakemake to handle and keep intermediate files
    """
    bin_input = Path(bin_input)
    (bin_output_ := Path(bin_output)).mkdir(parents=True, exist_ok=True)
    if bin_input.is_dir():
        assert list(bin_input.glob(f"*{support}")), "input is not a valid bin path"
        if str(support).endswith("fa") and keep_if_avail:
            suffix = str(support)
            bin_input_dir = bin_input
        else:
            suffix = "fa"
            bin_input_dir = bin_output_
            support_str_len = len(str(support))
            for bin_file in bin_input.glob(f"*{support}"):
                bin_name = bin_file.name[:-support_str_len].strip(".")
                shutil.copy(bin_file, bin_input_dir / f"{bin_name}.fa")
    else:
        assert bin_input.is_file() and Path(support).is_file()
        bin_input_dir = contig2bin(bin_output_, bin_input, support)
        suffix = "fa"
    return bin_input_dir, suffix


def checkm(
    bin_input: PathLike,
    support: Union[PathLike, str],
    output_dir=None,
    threads=10,
    **checkm_options,
):
    with TemporaryDirectory() as _td:
        file = f"{_td}/checkm.tsv"
        bin_input_, support_ = format_bin_input(
            bin_output=f"{_td}/bin_input",
            bin_input=bin_input,
            support=support,
            keep_if_avail=True,
        )
        if output_dir is None:
            output_dir = f"{_td}/checkm"
        CheckMFakeOptions(
            file=file,
            bin_input=str(bin_input_),
            output_dir=output_dir,
            extension=support_,
            threads=threads,
            **checkm_options,  # type: ignore  # confirmed by users
        ).run()
        return pd.read_csv(file, sep="\t")


def gunc(
    bin_input: PathLike,
    support: Union[PathLike, str],
    GUNC_DB: PathLike,
    output_dir: PathLike = None,
    threads=10,
):
    """
    @param support:
        if endswith "faa", will just use bin_input as directory of faa files.
    @param output_dir:
        if not None: will just detect output files in the output_dir and return
        should not be a existing directory

    include 3 steps:
    1.  extrat genome with "contig2bin"
    2.  annot gene with "prodigal_gff_multithread"
    3.  gunc

    """
    if output_dir:
        for gunc_tsv_file in Path(output_dir).glob("GUNC.*maxCSS_level.tsv"):
            return pd.read_csv(gunc_tsv_file, sep="\t")
    with TemporaryDirectory() as _td:
        (bin_faa_dir := Path(f"{_td}/out-bins_faa")).mkdir(parents=True, exist_ok=True)
        if str(support).endswith("faa"):
            bin_input_dir = Path(bin_input)
            for bin_file in bin_input_dir.glob(f"*{support}"):
                shutil.copy(bin_file, bin_faa_dir)
            suffix = str(support)
        else:
            bin_input_dir, suffix = format_bin_input(
                bin_output=f"{_td}/bin_fa_input",
                bin_input=bin_input,
                support=support,
                keep_if_avail=True,
            )
            prodigal_multithread(
                bin_input_dir.glob(f"*{suffix}"),
                mode="meta",
                out_dir=bin_faa_dir,
                suffix="faa",
                threads=threads,
            )

        # gunc_out_tsv = f"{_td}/out-gunc.tsv"
        # gunc_out_dir = f"{_td}/out-gunc-dir"

        smk_workflow = Path(__file__).parent.parent / "workflow"
        smk_conda_env = Path(__file__).parent.parent / ".snakemake" / "conda"
        target_smk_file = smk_workflow / "genome.smk"
        tpmf_outs = f"{_td}/out-gunc.tsv"
        smk_params = (
            f"-s {target_smk_file} "
            f"{tpmf_outs} "
            f"--config GUNC_DB='{GUNC_DB}' "
            f"--use-conda "
            f"--conda-prefix {smk_conda_env} "
            f"-c{threads} -rp "
        )
        try:
            print("params:", "snakemake", smk_params)
            smk(smk_params)
        except SystemExit as se:
            if se.code:
                print(se.code, se.with_traceback(None))
                raise RuntimeError("snakemake seems not run successfully.")

        if output_dir:
            Path(output_dir).mkdir(parents=True, exist_ok=True)
            for i in Path(f"{_td}/out-gunc-dir").glob("*"):
                shutil.move(i, output_dir)

        return pd.read_csv(tpmf_outs, sep="\t")


def bin_filter(
    bin_out_dir: PathLike,
    bin_input: PathLike,
    support: Union[PathLike, str],
    GUNC_DB: PathLike,
    checkm_output_dir: Optional[Union[PathLike, pd.DataFrame]] = None,
    gunc_output_dir: Optional[Union[PathLike, pd.DataFrame]] = None,
    threads=10,
):
    """
    filter input genomes carefully, applying two checks
    1.  checkm: 50, 10
    2.  gunc: passed (in single mode)

    other genomes will be put into bin_out_dir / "discard"

    intermediate checkm and gunc table can be provided as DataFrame;
    or will be used as output dir
    """
    with TemporaryDirectory() as _td:
        # region format input to prodigal single faa
        (bin_faa_dir := Path(f"{_td}/out-bins_faa")).mkdir(parents=True, exist_ok=True)
        if str(support).endswith("faa"):
            raise ValueError("genome nucleotide sequences required only")
        bin_input_dir, suffix = format_bin_input(
            bin_output=f"{_td}/bin_fa_input",
            bin_input=bin_input,
            support=support,
            keep_if_avail=False,
        )
        for bin_faa in prodigal_multithread(
            bin_input_dir.glob(f"*{suffix}"),
            mode="single",
            out_dir=bin_faa_dir,
            suffix="faa",
            threads=threads,
        ):
            bin_faa.rename(str(bin_faa)[:-20] + ".faa")
        # endregion format input to prodigal single faa

        if isinstance(checkm_output_dir, pd.DataFrame):
            checkm_tsv = checkm_output_dir
        else:
            checkm_output_dir_ = Path(checkm_output_dir or f"{_td}/checkm")
            checkm_output_dir_.mkdir(parents=True, exist_ok=True)
            checkm_output_file = checkm_output_dir_ / f"checkm.tsv"
            CheckMFakeOptions(
                file=str(bin_faa_dir),
                bin_input=str(bin_input_dir),
                output_dir=str(checkm_output_dir_),
                extension="faa",
                threads=threads,
            ).run()
            checkm_tsv = pd.read_csv(checkm_output_file, sep="\t")

        if isinstance(gunc_output_dir, pd.DataFrame):
            gunc_tsv = gunc_output_dir
        else:
            gunc_tsv = gunc(
                bin_input=str(bin_faa_dir),
                support="faa",
                GUNC_DB=GUNC_DB,
                output_dir=gunc_output_dir,
                threads=threads,
            )

        checkm_gunc = checkm_tsv.merge(
            gunc_tsv[["genome", "taxonomic_level", "pass.GUNC"]].rename(
                {"genome": "Bin Id"}, axis=1
            )
        )
        checkm_gunc_filter = checkm_gunc[
            (checkm_gunc["Completeness"] >= 50) & (checkm_gunc["Contamination"] <= 10)
        ]

        (bin_out_dir_ := Path(bin_out_dir)).mkdir(parents=True, exist_ok=True)
        for bin_fa in checkm_gunc_filter["Bin Id"]:
            shutil.move(bin_input_dir / f"{bin_fa}.fa", bin_out_dir_ / f"{bin_fa}.fa")

        shutil.move(bin_input_dir, bin_out_dir_ / "discard")
        checkm_gunc.to_csv(
            bin_out_dir_ / "discard" / "checkm_gunc.tsv", sep="\t", index=False
        )

        return checkm_gunc_filter
