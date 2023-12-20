# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-15 22:24:10
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-12-20 19:56:44
 * @FilePath: /genome/test/genome/test_gene_clust.py
 * @Description:
__file__ = "test/genome/test_prokka.py"
"""

from pathlib import Path

from genome.gene_clust import mmseq_clust, MmseqOut
from genome.gff import Parse


try:
    from ._decorator import temp_output, test_files
except (ModuleNotFoundError, ImportError):
    from _decorator import temp_output, test_files


@temp_output
def test_mmseq_clust(test_temp: Path):
    expect_prefix = test_temp / "some_gene"
    gffs = [
        test_files / "WIB11_CH1_genome-prodigal.single.gff",
    ]
    faa_iter = (faa for gff in gffs for faa in Parse(gff).extract())
    mo = mmseq_clust(
        [test_files / "some_gene-clu_rep.faa"], faa_iter, expect_prefix, threads=10
    )
