# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-25 16:45:32
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-11-24 16:27:09
 * @FilePath: /genome/genome/binning.py
 * @Description:
"""


import os
import shutil
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Final, Literal, NamedTuple, Optional, Union

import yaml
from snakemake import main as smk

from .bin_statistic import contig2bin
from .bin_statistic_ext import checkm

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


class BinningOutput(NamedTuple):
    ctg2mag: Path
    out_dir: Optional[Path] = None

    @classmethod
    def from_prefix(cls, prefix: PathLike):
        return cls(Path(f"{prefix}.tsv"), Path(f"{prefix}-dir"))


class BinningConfig(NamedTuple):
    MIN_BIN_CONTIG_LEN: int = 1500
    contig: str = "{any}.fa"
    bams: Optional[list[str]] = None
    lsbams: str = "{any}.bams.ls"
    jgi: str = "{any}-jgi.tsv"
    bin_single: str = "{any}-bins/single"
    bin_union_dir: str = "{any}-bins/union"
    bin_methods: Optional[list[str]] = None

    def to_config(self, config_file: PathLike):
        Path(config_file).parent.mkdir(parents=True, exist_ok=True)
        _asdict = self._asdict()
        _asdict["bams"] = _asdict["bams"] or []
        _asdict["bin_methods"] = _asdict["bin_methods"] or default_bin_methods
        with open(config_file, "w") as c:
            yaml.dump(_asdict, stream=c, allow_unicode=True)

    def touch_contig_jgi_bams(self):
        from Bio import SeqIO

        bin_single = Path(self.bin_single)
        bin_single.mkdir(parents=True, exist_ok=True)
        contig = bin_single / f"contig.{self.MIN_BIN_CONTIG_LEN}.fa"
        if not contig.is_file():
            SeqIO.write(
                (
                    i
                    for i in SeqIO.parse(self.contig, "fasta")
                    if self.MIN_BIN_CONTIG_LEN <= len(i.seq)
                ),
                contig,
                format="fasta",
            )
            os.system(f"touch -amcr {self.contig} {contig}")

        jgi = bin_single / f"vamb-jgi.tsv"
        if not jgi.is_file():
            with open(self.jgi) as fi:
                header = next(fi)
                jgi_lines: dict[str, str] = {i.split()[0]: i for i in fi}
            with open(jgi, "w") as fo, open(contig) as fa:
                fo.write(header)
                for line in fa:
                    if line.startswith(">"):
                        fo.write(jgi_lines[line[1:].split()[0]])
                fo.flush()
            os.system(f"touch -amcr {self.jgi} {jgi}")

        lsbams = bin_single / f"bams.ls"
        if Path(self.lsbams).is_file():
            shutil.copy(self.lsbams, lsbams)
            os.system(f"touch -amcr {self.lsbams} {lsbams}")

    def output(self, basename):
        bin_union_dir = Path(self.bin_union_dir)
        return BinningOutput.from_prefix(bin_union_dir / basename)

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
            self.touch_contig_jgi_bams()

            tmpc = self._replace(
                jgi=str(Path(self.bin_single) / f"vamb-jgi.tsv"),
                lsbams=str(Path(self.bin_single) / f"bams.ls"),
            )
            tmp_config = f"{_td}/config"
            tmpc.to_config(tmp_config)

            tpmf_out = self.output(out_basename)

            smk_workflow = Path(__file__).parent.parent / "workflow"
            smk_conda_env = Path(__file__).parent.parent / ".snakemake" / "conda"
            target_smk_file = smk_workflow / "binning" / "__init__.smk"

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


def check_bams(
    bin_union_dir: PathLike, bams: Union[list[PathLike], PathLike] = None
) -> tuple[PathLike, list[PathLike]]:
    if not bams:
        bams = []
        lsbams = Path(str(bin_union_dir) + "-bams.list")
    elif isinstance(bams, list):
        lsbams = Path(str(bin_union_dir) + "-bams.list")
    else:
        lsbams, bams = Path(bams), []
    return lsbams, bams


def bin_union(
    method: Literal[
        "dastool", "unitem_greedy", "unitem_consensus", "unitem_unanimous"
    ] = "dastool",
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
    """
    @param method: only given methods are supported

    @param marker: name for identify, will add after "dastool"
        - [if not given], name will be "dastool" only
        - else, name will be like "dastool-{marker}"

    @param bin_union_dir: dirname of output, so file will be like "{bin_union_dir}/dastool{-marker}.tsv"
        - [if not given], name will follow @param contig (if @param contig is skipped, will use ".")
        - else, will use given path

    @param contig: used as reference by following methods
        - in fact, it should not be skipped

    @param jgi, bams: used in single binning
        - can be skipped, if all single binning result is generated

    @param bin_single: to keep each single bin method result ctg2mag.tsv, should be maintained by user
        - [if not given], will be store in "{bin_union_dir}/single"

    @param bin_methods: bin methods to use
        - [if not given], will use all 9 metabat2, 2 maxbin2, concoct, vamb, and metadecoder

    @param min_bin_contig_len: min length of contig

    @param threads: threads
    """
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

    lsbams, bams = check_bams(bin_union_dir, bams)

    bc = BinningConfig(
        min_bin_contig_len,
        str(contig),
        [str(i) for i in bams],
        str(lsbams),
        str(jgi),
        str(bin_single),
        str(bin_union_dir),
        bin_methods,
    )
    return bc.run(method, marker, threads)
