# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-12 19:53:55
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-06-18 11:14:22
 * @FilePath: /genome/tests/genome/test_gff.py
 * @Description:
__file__ = "tests/genome/test_prokka.py"
"""

from genome.gff import Parse
from Bio import SeqIO

from tests import Path, temp_output, test_files, test_temp


@temp_output
def test_gff_extract_protein_fa(test_temp: Path):
    gff = test_files / "binny_contigs_4bins-top10-prodigal.gvmeta.gff"

    test_out = test_temp / "binny_contigs_4bins-top10-prodigal.gvmeta-ge33.faa"
    SeqIO.write(
        sorted(Parse(gff).extract(min_aa_length=33), key=lambda x: str(x.id)),
        test_out,
        "fasta-2line",
    )


@temp_output
def test_gff_reset_reference(test_temp: Path):
    gff = test_files / "binny_contigs_4bins-top10-prodigal.gvmeta.gff"
    genome = test_files / "binny_contigs_4bins.fa"

    test_out = test_temp / "binny_contigs_4bins-top10-prodigal.gvmeta-ge33.faa"
    SeqIO.write(
        sorted(
            Parse(genome, gff).extract(min_aa_length=33),
            key=lambda x: str(x.id),
        ),
        test_out,
        "fasta-2line",
    )
