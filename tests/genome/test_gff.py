# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-12 19:53:55
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-01-15 22:50:13
 * @FilePath: /genome/tests/genome/test_gff.py
 * @Description:
__file__ = "tests/genome/test_prokka.py"
"""

import os
from pathlib import Path
from genome.gff import Parse
from Bio import SeqIO

try:
    from _decorator import temp_output, test_temp, test_files
except (ModuleNotFoundError, ImportError):
    from tests.genome._decorator import temp_output, test_temp, test_files


@temp_output
def test_gff_extract_protein_fa(test_temp: Path):
    gff = test_files / "binny_contigs_4bins-top10-prodigal.gvmeta.gff"
    expect = test_files / "binny_contigs_4bins-top10-prodigal.gvmeta-ge33.faa"

    test_out = test_temp / "binny_contigs_4bins-top10-prodigal.gvmeta-ge33.faa"
    SeqIO.write(
        sorted(Parse(gff).extract(min_aa_length=33), key=lambda x: x.id),
        test_out,
        "fasta-2line",
    )

    os.system(f"md5sum {expect}")
    os.system(f"md5sum {test_out}")


@temp_output
def test_gff_reset_reference(test_temp: Path):
    gff = test_files / "binny_contigs_4bins-top10-prodigal.gvmeta.gff"
    genome = test_files / "binny_contigs_4bins.fa"
    expect = test_files / "binny_contigs_4bins-top10-prodigal.gvmeta-ge33.faa"

    test_out = test_temp / "binny_contigs_4bins-top10-prodigal.gvmeta-ge33.faa"
    SeqIO.write(
        sorted(
            Parse(gff, gff).extract(min_aa_length=33),
            key=lambda x: x.id,
        ),
        test_out,
        "fasta-2line",
    )

    os.system(f"md5sum {expect}")
    os.system(f"md5sum {test_out}")
