# -*- coding: utf-8 -*-
"""
* @Date: 2022-10-11 13:49:35
* @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
* @LastEditTime: 2025-07-09 15:31:47
* @FilePath: /genome/genome/prokka.py
* @Description:
"""

import os
import shutil
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Iterable, Literal, Union

from Bio import SeqIO, SeqRecord

from .pyrule import smk, rules_dir, smk_conda_env


PathLike = Union[str, Path]
prokka_kindom = Literal["Archaea", "Bacteria", "Mitochondria", "Viruses"]


def prokka_gff_onethread(
    genome: Union[PathLike, Iterable[SeqRecord.SeqRecord]],
    kingdom: prokka_kindom = "Bacteria",
    gff_out: PathLike = "",
) -> Path:
    # infer gff_out automatically if not given in some cases
    """
    Runs Prokka genome annotation on a single genome input using a single thread and outputs a GFF file.
    
    Parameters:
        genome: Either a file path to a FASTA file ending with `.fa` or an iterable of `SeqRecord` objects representing the genome sequence.
        kingdom: The kingdom classification for Prokka annotation (default is "Bacteria").
        gff_out: Optional output path for the resulting GFF file. If not provided, the output path is inferred from the genome filename.
    
    Returns:
        Path to the generated GFF file containing the Prokka annotation results.
    
    Raises:
        ValueError: If the genome input is not a valid file path or iterable of `SeqRecord`, or if the file does not end with `.fa`.
        FileNotFoundError: If the specified genome file does not exist.
        RuntimeError: If the Snakemake workflow fails to complete successfully.
        NotImplementedError: If the workflow does not complete as expected.
    """
    if not gff_out:
        if not isinstance(genome, str) and not isinstance(genome, Path):
            raise ValueError("initial filename must be provided")
        if not str(genome).endswith(".fa"):
            raise ValueError("genome file must endswith '.fa'")
        gff_out_ = Path(str(genome)[:-3] + f"-prokka_{kingdom}.gff")
    else:
        gff_out_ = Path(gff_out)

    with NamedTemporaryFile("w", suffix=".fa", delete=True) as tmpf:
        tpmf_out = Path(f"{tmpf.name[:-3]}-prokka_{kingdom}.gff")

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
        smk_params = f"-s {target_smk_file} " f"{tpmf_out} " f"--use-conda " f"-c1 -p"

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
    genomes: Iterable[PathLike],
    kingdom: prokka_kindom = "Bacteria",
    gff_out_dir: PathLike = "",
    threads: int = 8,
) -> Iterable[Path]:
    """
    Runs Prokka genome annotation on multiple FASTA files in parallel using Snakemake, generating GFF output files.
    
    If an output directory is specified, genome files are copied there before processing, and outputs are written to that directory. Duplicate genome filenames are not allowed. Cleans up temporary genome copies and log files after completion.
    
    Parameters:
        genomes: Iterable of input genome FASTA file paths (must end with `.fa`).
        kingdom: Kingdom classification for Prokka annotation (default is "Bacteria").
        gff_out_dir: Optional directory to store genome copies and GFF outputs. If not provided, outputs are written alongside input files.
        threads: Number of threads to use for parallel Snakemake execution (default is 8).
    
    Returns:
        Iterable of Paths to the generated GFF files.
    """
    # if many genomes are provided, the file must exist
    _genome_files = [Path(file).expanduser().absolute() for file in genomes]
    for genome in _genome_files:
        if not str(genome).endswith(".fa"):
            raise ValueError("genome file must endswith '.fa'")

    genome_files = []
    if gff_out_dir:
        gff_out_dir_ = Path(gff_out_dir).expanduser().absolute()
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
                os.system(f"touch -amcr {genome_path} {new_file}")
            genome_files.append(new_file)
    else:
        genome_files.extend(_genome_files)

    target_smk_file = rules_dir / "genome.smk"
    tpmf_outs = [f"{str(genome)[:-3]}-prokka_{kingdom}.gff" for genome in genome_files]
    tpmf_outs_str = " ".join(tpmf_outs)
    smk_params = (
        f"-s {target_smk_file} "
        f"{tpmf_outs_str} "
        f"--use-conda "
        f"--conda-prefix {smk_conda_env} "
        f"-c{threads} -p "
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
                Path(new_file).unlink()
    for tpmf_out in tpmf_outs:
        Path(f"{tpmf_out[:-3]}log").unlink()

    return [Path(tpmf_out) for tpmf_out in tpmf_outs]
