# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-25 16:45:32
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-11-16 22:08:47
 * @FilePath: /genome/genome/binning.py
 * @Description:
"""


import os
import shutil
from dataclasses import dataclass
from pathlib import Path
from tempfile import NamedTemporaryFile, TemporaryDirectory
from typing import Final, NamedTuple, Optional, Union, Literal

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

    def touch_contig_jgi(self):
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

    def get_bams(self):
        if Path(self.lsbams).is_file():
            with open(self.lsbams) as bams:
                return [Path(bam) for bam in bams.read().strip().split()]
        else:
            return [Path(bam) for bam in self.bams]

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
            # with NamedTemporaryFile("w", suffix=".yaml", delete=True) as tmpf:
            tmp_lsbams = f"{_td}/fake-bams.ls"
            tmp_bams = []
            if Path(self.lsbams).is_file():
                shutil.copy(self.lsbams, tmp_lsbams)
            else:
                if not self.bams:
                    raise FileNotFoundError("bams must be given if lsbams is not file")
                for i, bam in enumerate(self.bams):
                    tmp_bams.append(tmp_bam := f"{_td}/{i}.bam")
                    os.system(f"ln -s {Path(bam).expanduser().absolute()} {tmp_bam}")

            self.touch_contig_jgi()

            tmpc = self._replace(
                jgi=str(Path(self.bin_single) / f"vamb-jgi.tsv"),
                bams=tmp_bams,
                lsbams=tmp_lsbams,
            )
            tmp_cofig = f"{_td}/config"
            tmpc.to_config(tmp_cofig)

            tpmf_out = self.output(out_basename)

            smk_workflow = Path(__file__).parent.parent / "workflow"
            smk_conda_env = Path(__file__).parent.parent / ".snakemake" / "conda"
            target_smk_file = smk_workflow / "binning" / "__init__.smk"

            smk_params2 = (
                f"-s {target_smk_file} "
                f"{tpmf_out.ctg2mag} "
                f"--use-conda "
                f"--conda-prefix {smk_conda_env} "
                f"-c{threads} -rp "
                f"--configfile {tmp_cofig} "
            )

            try:
                os.system(f"ls {tmp_cofig}")
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
