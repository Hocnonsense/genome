# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-25 16:45:32
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-10-23 13:21:22
 * @FilePath: /genome/genome/binning.py
 * @Description:
"""


import os
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Final, Literal, NamedTuple

import yaml
from snakemake import main as smk

PathLike = str | Path
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


class BinningOutput(NamedTuple):
    ctg2mag: Path
    out_dir: Path | None = None

    @classmethod
    def from_prefix(cls, prefix: PathLike):
        return cls(Path(f"{prefix}.tsv"), Path(f"{prefix}-dir"))


class BinningInput(NamedTuple):
    contig: str = "{any}.fa"
    lsbams: str = "{any}.bams.ls"
    jgi: str = "{any}-jgi.tsv"

    @classmethod
    def file_from_prefix(cls, prefix: PathLike):
        return Path(str(prefix) + "-bins.yaml")

    def dump_prefix(self, prefix: PathLike, relpath=True):
        config = self.file_from_prefix(prefix)
        config.parent.mkdir(parents=True, exist_ok=True)
        if relpath:
            _asdict = {k: os.path.relpath(v, config.parent) for k, v in self._asdict().items()}
        else:
            _asdict = self._asdict()

        with open(config, "w") as fo:
            yaml.safe_dump(_asdict, fo)
        return config


@dataclass
class BinningConfig:
    MIN_BIN_CONTIG_LEN: int = 1500
    bin_methods: list[str] | None = None
    bin_config: PathLike = "{any}"

    def __post_init__(self):
        bin_config = str(self.bin_config)
        assert bin_config.endswith("-bins.yaml")
        self._bin_perfix_dir = bin_config[:-5]

    def to_config(self, config_file: PathLike):
        Path(config_file).parent.mkdir(parents=True, exist_ok=True)
        _asdict = dict(
            MIN_BIN_CONTIG_LEN=self.MIN_BIN_CONTIG_LEN,
            bin_methods=self.bin_methods or default_bin_methods,
        )
        with open(config_file, "w") as c:
            yaml.dump(_asdict, stream=c, allow_unicode=True)
        return Path(config_file)

    def output(self, basename: str = "{method}[-{marker}]"):
        return BinningOutput.from_prefix(
            Path(self._bin_perfix_dir) / "union" / basename
        )

    def run(
        self,
        method: Literal[
            "dastool", "unitem_greedy", "unitem_consensus", "unitem_unanimous"
        ],
        marker: str = "",
        threads=10,
    ):
        """
        @param method: only given methods are supported

        @param marker: name for identify, will add after "dastool"
            - [if not given], name will be "dastool" only
            - else, name will be like "dastool-{marker}"
        """
        out_basename: str = method
        if marker:
            out_basename += f"-{marker}"

        with TemporaryDirectory() as _td:
            assert Path(self.bin_config).is_file()

            tpmf_out = self.output(out_basename)
            tmp_config = self.to_config(f"{_td}/config")

            smk_workflow = Path(__file__).parent.parent / "workflow"
            smk_conda_env = Path(__file__).parent.parent / ".snakemake" / "conda"
            target_smk_file = smk_workflow / "binning" / "genomecall.smk"

            smk_params2 = (
                f"-s {target_smk_file} "
                f"{tpmf_out.ctg2mag} "
                f"--nolock "
                # f"--drop-metadata "  # add this if necessary
                f"--use-conda "
                f"--conda-prefix {smk_conda_env} "
                f"-c{threads} -rp "
                f"--configfile {tmp_config} "
            )

            try:
                os.system(f"ls {tmp_config}")
                print("params:", "snakemake", smk_params2)
                smk(smk_params2)
            except SystemExit as se:
                if se.code:
                    print(se.code, se.with_traceback(None))
                    raise RuntimeError("snakemake seems not run successfully.")
                else:
                    return tpmf_out

        raise NotImplementedError("")


def check_bams(bin_prefix: PathLike, bams: list[PathLike] | PathLike | None = None):
    if not bams:
        bams = []
        lsbams = Path(str(bin_prefix) + "-bams.list")
    elif isinstance(bams, list):
        lsbams = Path(str(bin_prefix) + "-bams.list")
    else:
        lsbams, bams = Path(bams), []
    if not lsbams.is_file():
        lsbams.parent.mkdir(parents=True, exist_ok=True)
        with open(lsbams, "w") as fo:
            for i in bams:
                print(i, file=fo)
    return str(lsbams)


def bin_union(
    method: Literal[
        "dastool", "unitem_greedy", "unitem_consensus", "unitem_unanimous"
    ] = "dastool",
    marker: str = "",
    bin_prefix: PathLike = "",
    contig: PathLike = "",
    jgi: PathLike = "",
    bams: list[PathLike] | PathLike | None = None,
    bin_methods: list[str] | None = None,
    min_bin_contig_len: int = 1500,
    threads: int = 10,
):
    """
    @param method: only given methods are supported

    @param marker: name for identify, will add after "dastool"
        - [if not given], name will be "dastool" only
        - else, name will be like "dastool-{marker}"

    @param bin_prefix: prefix of output, so file will be like "{bin_prefix}-bins/dastool{-marker}.tsv"
        - [if not given], name will follow @param contig (if @param contig is skipped, will use "./bins")
        - else, will use given path

    @param contig: used as reference by following methods
        - in fact, it should not be skipped

    @param jgi, bams: used in single binning
        - can be skipped, if all single binning result is generated

    @param bin_methods: bin methods to use
        - [if not given], will use all 9 metabat2, 2 maxbin2, concoct, vamb, and metadecoder

    @param min_bin_contig_len: min length of contig

    @param threads: threads
    """
    # infer bin methods automatically if not given in some cases
    if not bin_methods:
        bin_methods = default_bin_methods
    min_bin_contig_len = max(int(min_bin_contig_len), AVAIL_MIN_BIN_CONTIG_LEN)

    if not Path(contig).is_file() and not bin_prefix:
        bin_prefix = "bins"
    if not bin_prefix:
        bin_prefix = str(contig)

    lsbams = check_bams(bin_prefix, bams)

    bc = BinningConfig(
        MIN_BIN_CONTIG_LEN=min_bin_contig_len,
        bin_methods=bin_methods,
        bin_config=BinningInput(
            contig=str(contig),
            lsbams=str(lsbams),
            jgi=str(jgi),
        ).dump_prefix(prefix=bin_prefix),
    )
    return bc.run(method, marker, threads)
