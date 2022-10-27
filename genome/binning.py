# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-25 16:45:32
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-10-27 19:11:54
 * @FilePath: /genome/genome/binning.py
 * @Description:
"""


import os
from dataclasses import dataclass
from pathlib import Path
from tempfile import NamedTemporaryFile, TemporaryDirectory
from typing import Final, NamedTuple, Optional, Union, overload

import pandas as pd
import yaml
from snakemake import main as smk

from genome.gff import Parse

PathLike = Union[str, Path]
AVAIL_MIN_BIN_CONTIG_LEN: Final = 1000
default_bin_methods = [
    "metabat2_60_60",
    "metabat2_60_75",
    "metabat2_60_90",
    "metabat2_75_60",
    "metabat2_75_75",
    "metabat2_75_90",
    "metabat2_90_60",
    "metabat2_90_75",
    "metabat2_90_90",
    "maxbin2_40",
    "maxbin2_107",
    "concoct",
    "metadecoder",
    "vamb",
]


class BinningConfig(NamedTuple):
    MIN_BIN_CONTIG_LEN: int = 1500
    contig: str = "data/{site}/02_assem..{site}_cut.fa"
    bams: list[str] = ["data/{site}/02_assem..{site}..{layer}_cut.bam"]
    bams_ls: str = "data/{site}/02_assem..{site}_cut.bam.list"
    jgi: str = "data/{site}/02_assem..{site}_cut-jgi.depth"
    bin_single: str = "pipe/{site}/04_bin/single"
    bin_union_dir: str = "pipe/{site}/04_bin/"
    bin_methods: list[str] = default_bin_methods

    def to_config(self, config_file: PathLike):
        Path(config_file).parent.mkdir(parents=True, exist_ok=True)
        with open(config_file, "w") as c:
            yaml.dump(self._asdict(), stream=c, allow_unicode=True)


def DAS_Tool(
    marker: str = "",
    bin_union_dir: PathLike = None,
    contig: PathLike = "",
    jgi: PathLike = "",
    bams: Union[list[PathLike], PathLike] = None,
    bin_single: PathLike = None,
    bin_methods: list[str] = None,
    min_bin_contig_len: int = 1500,
    threads: int = 10,
):
    dastool_out = "dastool"
    if marker:
        dastool_out += f"-{marker}"
    # infer bin methods automatically if not given in some cases
    if not bin_methods:
        bin_methods = default_bin_methods
    min_bin_contig_len = max(int(min_bin_contig_len), AVAIL_MIN_BIN_CONTIG_LEN)

    if not Path(contig).is_file() and not bin_union_dir:
        bin_union_dir = Path(".")
    if not bin_union_dir:
        bin_union_dir = Path(str(contig) + "-bins")
    if not bin_single:
        bin_single = Path(bin_union_dir) / "single"
    if not bams:
        bams = []
        bams_ls = Path(str(bin_union_dir) + "-bams.list")
    elif isinstance(bams, list):

        bams_ls = Path(str(bin_union_dir) + "-bams.list")
    else:
        bams_ls, bams = Path(bams), []
    bc = BinningConfig(
        min_bin_contig_len,
        str(contig),
        [str(i) for i in bams],
        str(bams_ls),
        str(jgi),
        str(bin_single),
        str(bin_union_dir),
        bin_methods,
    )
    with NamedTemporaryFile("w", suffix=".yaml", delete=True) as tmpf:
        # setattr(tmpf, "__str__", lambda self: self.name)
        # tmpf = "test/temp/test_binning.yaml"
        tpmf_out = Path(bin_union_dir) / f"{dastool_out}.tsv"
        bc.to_config(tmpf.name)

        smk_workflow = Path(__file__).parent.parent / "workflow"
        smk_conda_env = Path(__file__).parent.parent / ".snakemake" / "conda"
        target_smk_file = smk_workflow / "binning" / "__init__.smk"
        smk_params1 = (
            f"-s {target_smk_file} "
            f"filtered_contig "
            f"--use-conda "
            f"--conda-prefix {smk_conda_env} "
            f"-c{threads} -rp "
            f"--configfile {tmpf.name} "
        )

        try:
            print("params:", "snakemake", smk_params1)
            smk(smk_params1)
        except SystemExit as se:
            if se.code:
                print(se.code, se.with_traceback(None))
                raise RuntimeError("snakemake seems not run successfully.")
            else:
                smk_params2 = (
                    f"-s {target_smk_file} "
                    f"{tpmf_out} "
                    f"--use-conda "
                    f"--conda-prefix {smk_conda_env} "
                    f"-c{threads} -rp "
                    f"--configfile {tmpf.name} "
                )

        try:
            os.system(f"ls {tmpf.name}")
            print("params:", "snakemake", smk_params2)
            smk(smk_params2)
        except SystemExit as se:
            if se.code:
                print(se.code, se.with_traceback(None))
                raise RuntimeError("snakemake seems not run successfully.")
            else:
                return tpmf_out

    raise NotImplementedError("")


def contig2bin(outdir: PathLike, contig2bin_tsv: PathLike, contigs: PathLike):
    contig2bin_ = pd.read_csv(
        contig2bin_tsv, sep="\t", names=["contig", "bin"], index_col=0
    )

    td = Path(outdir)
    td.mkdir(parents=True, exist_ok=True)

    try:
        binfiles = {b: open(td / f"{b}.fa", "w") for b in contig2bin_["bin"].unique()}
        for i in Parse(contigs)():
            if i.name in contig2bin_.index:
                binfiles[contig2bin_.loc[i.name, "bin"]].write(i.format("fasta-2line"))
    finally:
        for bf in binfiles.values():
            bf.close()

    return td


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

        from checkm.main import OptionsParser
        from checkm.logger import logger_setup

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
