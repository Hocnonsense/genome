# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-11 13:49:35
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-10-11 16:25:53
 * @FilePath: /genome/genome/prokka.py
 * @Description:
"""

import os
import shutil
from pathlib import Path
from typing import Iterable, Literal, Union
from Bio import SeqRecord, SeqIO
from snakemake import main as smk
from tempfile import NamedTemporaryFile


prokka_kindom = Literal["Archaea", "Bacteria", "Mitochondria", "Viruses"]


def prokka_gff_onethread(
    genome: Union[str, Path, Iterable[SeqRecord.SeqRecord]],
    kingdom: prokka_kindom = "Bacteria",
    gff_out: Union[str, Path] = "",
) -> Path:
    # infer gff_out automatically if not given in some cases
    if not gff_out:
        if not isinstance(genome, str) and not isinstance(genome, Path):
            raise ValueError("inital filename must provided")
        if not str(genome).endswith(".fa"):
            raise ValueError("genome file must endswith '.fa'")
        gff_out_ = Path(str(genome)[:-3] + f"-prokka.{kingdom}.gff")
    else:
        gff_out_ = Path(gff_out)

    with NamedTemporaryFile("w", suffix=".fa", delete=True) as tmpf:
        tpmf_out = Path(f"{tmpf.name[:-3]}-prokka.{kingdom}.gff")

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
