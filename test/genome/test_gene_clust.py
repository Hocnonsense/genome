# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-15 22:24:10
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-10-15 22:37:09
 * @FilePath: /genome/test/genome/test_gene_clust.py
 * @Description:
__file__ = "test/genome/test_prokka.py"
"""

from pathlib import Path

from genome.gene_clust import mmseq_clust, MmseqOut
from genome.gff import Parse

test_temp = Path(__file__).parent.parent / "temp"
test_files = Path(__file__).parent.parent / "file"


def test_mmseq_clust():
    expect_prefix = test_files / "some_gene"
    gffs = [
        test_files / "metadecoder.2-prokka.Bacteria.gff",
        test_files / "metadecoder.3-prokka.Bacteria.gff",
        test_files / "metadecoder.4-prodigal.simple.gff",
    ]
    faa_iter = (faa for gff in gffs for faa in Parse(gff).extract())
    mo = mmseq_clust(
        [test_files / "metadecoder.1.prokka.faa"], faa_iter, expect_prefix, threads=10
    )
