# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-12 16:35:45
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-11-25 23:03:19
 * @FilePath: /genome/genome/prodigal.py
 * @Description:
"""


import os
import shutil
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Iterable, Literal, Union, Optional

from Bio import SeqIO, SeqRecord
from snakemake import main as smk


PathLike = Union[str, Path]
prodigal_mode = Literal["single", "meta"]


def check_genome_length_prodigal(
    genome: Union[PathLike, Iterable[SeqRecord.SeqRecord]]
):
    """check if givem genome size is long enough to use prodigal single mode"""
    if not isinstance(genome, str) and not isinstance(genome, Path):
        genome_iter = genome
    else:
        if not Path(genome).is_file():
            raise FileNotFoundError(f"file {genome} does not exist, please check.")
        genome_iter = SeqIO.parse(genome, "fasta")
    return sum(len(i) for i in genome_iter) >= 20000


def prodigal_gff_onethread(
    genome: Union[PathLike, Iterable[SeqRecord.SeqRecord]],
    mode: prodigal_mode = "single",
    gff_out: PathLike = "",
) -> Optional[Path]:
    # infer gff_out automatically if not given in some cases
    if not gff_out:
        if not isinstance(genome, str) and not isinstance(genome, Path):
            raise ValueError("inital filename must provided")
        if not str(genome).endswith(".fa"):
            raise ValueError("genome file must endswith '.fa'")
        gff_out_ = Path(str(genome)[:-3] + f"-prodigal.{mode}.gff")
    else:
        gff_out_ = Path(gff_out)

    with NamedTemporaryFile("w", suffix=".fa", delete=True) as tmpf:
        tpmf_out = Path(f"{tmpf.name[:-3]}-prodigal.{mode}.gff")

        if not isinstance(genome, str) and not isinstance(genome, Path):
            SeqIO.write(genome, tmpf, "fasta")
        else:
            if not Path(genome).is_file():
                raise FileNotFoundError(f"file {genome} does not exist, please check.")
            with open(genome) as fi:
                while True:
                    # read 16 Kib words one time
                    block = fi.read(65536)
                    if not block:
                        break
                    tmpf.write(block)
                tmpf.flush()
        if mode == "single" and not check_genome_length_prodigal(genome):
            return None

        smk_workflow = Path(__file__).parent.parent / "workflow"
        smk_conda_env = Path(__file__).parent.parent / ".snakemake" / "conda"
        target_smk_file = smk_workflow / "genome.smk"
        smk_params = (
            f"-s {target_smk_file} "
            f"{tpmf_out} "
            f"--use-conda "
            f"--conda-prefix {smk_conda_env} "
            f"-c1 -rp "
        )

        try:
            os.system(f"ls {tmpf.name}")
            print("params:", "snakemake", smk_params)
            smk(smk_params)
        except SystemExit as se:
            if se.code:
                print(se.code, se.with_traceback(None))
                raise RuntimeError("snakemake seems not run successfully.")
            else:
                gff_out_.parent.mkdir(parents=True, exist_ok=True)
                shutil.move(tpmf_out, gff_out_)
                return gff_out_

    raise NotImplementedError("")


def prodigal_multithread(
    genomes: Iterable[PathLike],
    mode: prodigal_mode = "single",
    out_dir: PathLike = "",
    suffix: Literal["gff", "faa", "fna"] = "gff",
    threads: int = 8,
) -> Iterable[Path]:
    """
    If in single mode, length should not shorter than 20000 bp.
    """
    # if many genomes are provided, the file must exist
    if mode == "single":
        _genome_files = [
            Path(file).expanduser().absolute()
            for file in genomes
            if check_genome_length_prodigal(file)
        ]
    else:
        _genome_files = [
            Path(file).expanduser().absolute()
            for file in genomes
            if check_genome_length_prodigal(file)
        ]
    if not _genome_files:
        return []
    for genome in _genome_files:
        if not str(genome).endswith(".fa"):
            raise ValueError("genome file must endswith '.fa'")

    genome_files: list[Path] = []
    if out_dir:
        gff_out_dir_ = Path(out_dir).expanduser().absolute()
        gff_out_dir_.mkdir(parents=True, exist_ok=True)

        _genome_files_dict = {file.name: file for file in _genome_files}
        if len(_genome_files_dict) != len(_genome_files):
            raise ValueError("cannot collect genome_files with same name")

        for genome_name, genome_path in _genome_files_dict.items():
            new_file = gff_out_dir_ / genome_name
            if new_file.exists():
                if new_file != genome_path:
                    raise FileExistsError(new_file)
            else:
                shutil.copy(genome_path, new_file)
            genome_files.append(new_file)
    else:
        genome_files.extend(_genome_files)

    smk_workflow = Path(__file__).parent.parent / "workflow"
    smk_conda_env = Path(__file__).parent.parent / ".snakemake" / "conda"
    target_smk_file = smk_workflow / "genome.smk"
    tpmf_outs = [
        f"{str(genome)[:-3]}-prodigal.{mode}.{suffix}" for genome in genome_files
    ]
    tpmf_outs_str = " ".join(tpmf_outs)
    smk_params = (
        f"-s {target_smk_file} "
        f"{tpmf_outs_str} "
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

    # clean up temp added genomes
    s_genome_files = [str(i) for i in _genome_files]
    if out_dir:
        for new_file in genome_files:
            if str(new_file) not in s_genome_files:
                new_file.unlink()

    return [Path(tpmf_out) for tpmf_out in tpmf_outs]
