# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-11 20:11:35
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-12-21 20:09:59
 * @FilePath: /genome/tests/genome/test_prodigal.py
 * @Description:
__file__ = "test/genome/test_prodigal.py"
"""

from pathlib import Path
import shutil
from genome.prodigal import prodigal_gff_onethread

try:
    from _decorator import temp_output, test_temp, test_files
except (ModuleNotFoundError, ImportError):
    from tests.genome._decorator import temp_output, test_temp, test_files


@temp_output
def test_prodigal_gff_onethread(test_temp: Path):
    genome = test_files / "binny_contigs_4bins.fa"
    expect = test_files / "binny_contigs_4bins-prodigal.single.gff"

    test_out = test_temp / "binny_contigs_4bins-prodigal.single.gff"
    gff = prodigal_gff_onethread(shutil.copy(genome, test_temp), "single")
    assert gff == test_out.expanduser().absolute()

    assert gff.read_text() == expect.read_text()


@temp_output
def test_pyrodigal(test_temp: Path):
    genome = test_files / "binny_contigs_4bins.fa"
    from Bio import SeqIO
    import pyrodigal
    import pyrodigal_gv

    seqs = SeqIO.to_dict(SeqIO.parse(genome, "fasta"))

    gf = pyrodigal.GeneFinder(meta=True)

    gfv = pyrodigal_gv.ViralGeneFinder(meta=True)

    with open("a", "w") as fa, open("b", "w") as fb:
        for i in seqs.values():
            gf.find_genes(bytes(i.seq)).write_gff(
                fa, i.id, include_translation_table=True
            )
            gfv.find_genes(bytes(i.seq)).write_gff(
                fb, i.id, include_translation_table=True
            )
    # gene1 = genes[0]
    # gene1.begin

    gf = pyrodigal.GeneFinder(meta=False)
    gf.train(*(bytes(i.seq) for i in seqs.values()))
