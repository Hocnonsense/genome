# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-12 16:35:45
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2025-01-13 16:39:12
 * @FilePath: /genome/genome/prodigal.py
 * @Description:
"""


import shutil
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Final, Iterable, Literal

from Bio import SeqIO, SeqRecord
import pyrodigal

from .pyrule import smk, rules_dir, smk_conda_env


PathLike = str | Path
prodigal_mode: Final = ["single", "meta", "gvmeta"]


def check_genome_length_prodigal(genome: PathLike | Iterable[SeqRecord.SeqRecord]):
    """
    Determine if a genome sequence meets the minimum length requirement for Prodigal single mode.
    
    Parameters:
    	genome: A file path to a FASTA file or an iterable of SeqRecord objects representing genome sequences.
    
    Returns:
    	bool: True if the total sequence length is at least 20,000 base pairs, otherwise False.
    """
    if isinstance(genome, str) or isinstance(genome, Path):
        genome_iter = SeqIO.parse(genome, "fasta")
    else:
        genome_iter = genome
    return sum(len(i) for i in genome_iter) >= 20000


def prodigal_gff_onethread(
    genome: PathLike | Iterable[SeqRecord.SeqRecord],
    mode: Literal["single", "meta", "gvmeta"] = "single",
    gff_out: PathLike = "",
    trans_table=11,
) -> Path:
    # infer gff_out automatically if not given in some cases
    """
    Run Prodigal gene prediction on a single genome and output results in GFF format.
    
    Parameters:
        genome: Path to a FASTA file or an iterable of SeqRecord objects representing the genome sequences.
        mode: Prodigal mode to use; one of "single", "meta", or "gvmeta". Determines the gene finder configuration.
        gff_out: Optional output path for the GFF file. If not provided, the output filename is inferred from the genome input.
        trans_table: Translation table number to use for gene prediction (default is 11).
    
    Returns:
        Path to the generated GFF output file containing gene predictions and input sequences.
    """
    if not gff_out:
        if not isinstance(genome, str) and not isinstance(genome, Path):
            raise ValueError("without gff output, initial filename must be provided")
        if not str(genome).endswith(".fa"):
            raise ValueError("without gff output, genome file must endswith '.fa'")
        gff_out_ = Path(str(genome)[:-3] + f"-prodigal_{mode}.gff")
    else:
        gff_out_ = Path(gff_out)

    if not isinstance(genome, str) and not isinstance(genome, Path):
        seqs: Iterable[SeqRecord.SeqRecord] = genome  # type: ignore [union-attr]
    else:
        seqs = SeqIO.parse(genome, "fasta")

    if mode == "gvmeta":
        import pyrodigal_gv

        gf: pyrodigal.GeneFinder = pyrodigal_gv.ViralGeneFinder(meta=True, mask=True)
    else:
        if mode == "meta":
            gf = pyrodigal.GeneFinder(meta=True, mask=True)
        elif mode == "single":
            gf = pyrodigal.GeneFinder(meta=False, mask=True, closed=True)
            seqs = list(seqs)
            gf.train(*(bytes(i.seq) for i in seqs), translation_table=trans_table)

    with NamedTemporaryFile("w+", suffix=".fa", delete=True) as tmpf:
        tpmf_out = Path(f"{tmpf.name[:-3]}-prodigal_{mode}.gff")
        with open(tpmf_out, "w") as fo:
            for i in seqs:
                gf.find_genes(bytes(i.seq)).write_gff(
                    fo, str(i.id), include_translation_table=True
                )
                SeqIO.write(i, tmpf, format="fasta-2line")
            print("##FASTA", file=fo)
            tmpf.flush()
            tmpf.seek(0)
            while True:
                # read 16 Kib words one time
                block = tmpf.read(65536)
                if not block:
                    break
                fo.write(block)

        shutil.move(tpmf_out, gff_out_)
        return gff_out_


def prodigal_multithread(
    genomes: Iterable[PathLike],
    mode: Literal["single", "meta", "gvmeta"] = "single",
    out_dir: PathLike = "",
    suffix="gff",
    threads: int = 8,
) -> Iterable[Path]:
    """
    Run Prodigal gene prediction on multiple genome files in parallel using Snakemake.
    
    Filters input genome files to those at least 20,000 bp in length (required for "single" mode), validates file extensions, and optionally copies genomes to an output directory. Constructs output filenames based on mode and suffix, then executes a Snakemake workflow for parallel gene prediction. Cleans up temporary genome copies if created.
    
    Parameters:
        genomes: Iterable of genome file paths to process.
        mode: Prodigal mode to use; one of "single", "meta", or "gvmeta".
        out_dir: Optional directory to copy genome files and store results.
        suffix: Output file suffix; supports "gff", "-ge33.faa", or "-ge33.fna".
        threads: Number of parallel threads to use.
    
    Returns:
        List of paths to the generated output files.
    """
    # if many genomes are provided, the file must exist
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

        for genome_path in _genome_files:
            new_file = gff_out_dir_ / genome_path.name
            if new_file.exists():
                if new_file != genome_path:
                    raise FileExistsError(new_file)
            else:
                shutil.copy(genome_path, new_file)
            genome_files.append(new_file)
    else:
        genome_files.extend(_genome_files)

    target_smk_file = rules_dir / "genome.smk"

    # region quick fix suffix
    if suffix in ["faa", "fna"]:
        suffix = "-ge33." + suffix
    if "." not in suffix:
        suffix = "." + suffix
    # endregion quick fix suffix

    tpmf_outs = [
        f"{str(genome)[:-3]}-prodigal_{mode}{suffix}" for genome in genome_files
    ]
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
    s_genome_files = [str(i) for i in _genome_files]
    if out_dir:
        for new_file in genome_files:
            if str(new_file) not in s_genome_files:
                new_file.unlink()

    return [Path(tpmf_out) for tpmf_out in tpmf_outs]
