# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-12 16:35:45
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-10-15 21:30:18
 * @FilePath: /genome/genome/prodigal.py
 * @Description:
"""


import os
import shutil
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Iterable, Literal, Union

from Bio import SeqIO, SeqRecord
from snakemake import main as smk


PathLike = Union[str, Path]
prodigal_mode = Literal["single", "meta"]


def prodigal_gff_onethread(
    genome: Union[PathLike, Iterable[SeqRecord.SeqRecord]],
    mode: prodigal_mode = "single",
    gff_out: PathLike = "",
) -> Path:
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

        smk_workflow = Path(__file__).parent.parent / "workflow"
        target_smk_file = smk_workflow / "genome.smk"
        smk_params = f"-s {target_smk_file} " f"{tpmf_out} " f"--use-conda " f"-c1 -rp"

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
