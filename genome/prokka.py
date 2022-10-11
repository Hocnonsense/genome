# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-11 13:49:35
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-10-11 22:57:55
 * @FilePath: /genome/genome/prokka.py
 * @Description:
"""

import os
import shutil
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Callable, Generator, Iterable, Literal, Union

from BCBio.GFF import parse
from Bio import SeqIO, SeqRecord, SeqFeature
from snakemake import main as smk


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


def prokka_gff_multithread(
    genomes: Iterable[Union[str, Path]],
    kingdom: prokka_kindom = "Bacteria",
    gff_out_dir: Union[str, Path] = "",
    threads: int = 8,
) -> Path:
    """WARNING: untested"""
    # if many genomes are provided, the file must exist
    _genome_files = [Path(file).expanduser().absolute() for file in genomes]
    for genome in _genome_files:
        if not str(genomes).endswith(".fa"):
            raise ValueError("genome file must endswith '.fa'")

    genome_files = []
    if gff_out_dir:
        gff_out_dir_ = Path(gff_out_dir).expanduser().absolute()
        gff_out_dir_.mkdir(parents=True, exist_ok=True)

        _genome_files_dict = {file.name: file for file in _genome_files}
        if len(_genome_files_dict) != len(_genome_files):
            raise ValueError("cannot collect genome_files with same name")

        for genome, genome_path in _genome_files_dict.items():
            new_file = gff_out_dir_/genome
            if new_file.exists():
                if new_file != genome_path:
                    raise FileExistsError(genome_path)
            else:
                shutil.copy(genome_path, new_file)
            genome_files.append(new_file)
    else:
        genome_files.extend(_genome_files)

    smk_workflow = Path(__file__).parent.parent / "workflow"
    target_smk_file = smk_workflow / "genome.smk"
    tpmf_outs = [Path(f"{genome.name[:-3]}-prokka.{kingdom}.gff") for genome in genomes]
    tpmf_outs_str = " ".join(tpmf_outs)
    smk_params = (
        f"-s {target_smk_file} "
        f"{tpmf_outs_str} "
        f"--use-conda "
        f"-c{threads} -rp"
    )
    try:
        print("params:", "snakemake", smk_params)
        smk(smk_params)
    except SystemExit as se:
        if se.code:
            print(se.code, se.with_traceback(None))
            raise RuntimeError("snakemake seems not run successfully.")

    # clean up temp added genomes
    if gff_out_dir:
        for new_file in genome_files:
            if new_file not in _genome_files:
                shutil.rmtree(new_file)

    return tpmf_outs


gff_out_format = Literal["faa", "fna"]

def prokka_gff_extract_protein_fa(gff_file, out_format:gff_out_format="faa"):
    with open(gff_file) as in_handle:
        rec: SeqRecord.SeqRecord = None
        for rec in parse(in_handle):
            fet: SeqFeature.SeqFeature = None
            for fet in rec.features:
                if fet.type != "CDS":
                    continue
                seq: SeqRecord.SeqRecord = fet.extract(rec)
                if out_format == "faa":
                    seq = seq.translate()
                seq.id = rec.id + "_" + str(int(fet.id.rsplit('_', 1)[1]))
                seq.description = " # ".join((str(i) for i in (
                    "", fet.location.start, fet.location.end, fet.location.strand,
                    ";".join(f"{k}={','.join(v)}" for k, v in fet.qualifiers.items())
                ))).strip()
                yield seq
