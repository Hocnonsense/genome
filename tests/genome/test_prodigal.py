# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-11 20:11:35
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-04-27 18:40:43
 * @FilePath: /genome/tests/genome/test_prodigal.py
 * @Description:
__file__ = "test/genome/test_prodigal.py"
"""

from pathlib import Path
from Bio import SeqIO
from genome.bin_statistic_ext import format_bin_input
from genome.prodigal import prodigal_gff_onethread, prodigal_multithread

try:
    from _decorator import temp_output, test_temp, test_files
except (ModuleNotFoundError, ImportError):
    from tests.genome._decorator import temp_output, test_temp, test_files


@temp_output
def test_prodigal_gff_onethread(test_temp: Path):
    genome = test_files / "binny_contigs_4bins.fa"

    SeqIO.write(
        (i for i, _ in zip(SeqIO.parse(genome, "fasta"), range(10))),
        test_out_fa := test_temp / "binny_contigs_4bins-top10.fa",
        "fasta",
    )
    expect = test_files / "binny_contigs_4bins-top10-prodigal_gvmeta.gff"

    gff = prodigal_gff_onethread(test_out_fa, "gvmeta")
    assert gff == (test_temp / expect.name).expanduser().absolute()

    # assert gff.read_text() == expect.read_text()


@temp_output
def test_prodigal_multithread(test_temp: Path):
    ctg2mag = test_files / "binny_unitem_unanimous.tsv"
    support = test_files / "binny_contigs_4bins.fa"
    bin_input_dir, binids, suffix = format_bin_input(
        bin_output=test_temp / "prodigal",
        bin_input=ctg2mag,
        support=support,
        keep_if_avail=False,
    )
    for bin_faa in prodigal_multithread(
        (bin_input_dir / (binid + suffix) for binid in binids),
        mode="single",
        out_dir=bin_input_dir,
        suffix="-ge33.faa",
        threads=2,
    ):
        bin_faa.rename(str(bin_faa)[:-25] + ".faa")
