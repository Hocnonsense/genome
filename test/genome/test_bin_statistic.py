# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-12 19:53:55
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-10-15 19:08:03
 * @FilePath: /genome/test/genome/test_bin_statistic.py
 * @Description:
__file__ = "/home/hwrn/software/genome/test/genome/test_prokka.py"
"""

from pathlib import Path

import pandas as pd

from genome.bin_statistic import BinStatisticContainer

test_temp = Path(__file__).parent.parent / "temp"
test_files = Path(__file__).parent.parent / "file"


def test_prokka_gff_bin_statistic():
    gff = test_files / "metadecoder.1-prokka.Bacteria.gff"
    bsc = BinStatisticContainer.read_gff(gff)
    bs = bsc.statistic()
    assert int(bs.gc * 100) == 56
    assert int(bs.gc_std * 10000) == 206
    assert bs.bp_size == 1040626
    assert bs.max_contig_len == 9951
    assert bs.contigs_num == 293
    assert bs.contig_n50 == 3411
    assert bs.ambiguous_bases_num == 0
    assert int(bs.coding_density * 100) == 78
    assert bs.genes_num == 1130


def test_genome_bin_statistic():
    genome = test_files / "metadecoder.4.fa"
    bsc = BinStatisticContainer.read_contig(genome, "fasta")
    test_contig_cutoffs = (
        0,
        500,
        1000,
        1500,
        2000,
        2500,
        3000,
        4000,
        5000,
        10000,
    )
    df = pd.DataFrame(
        {
            min_contig_len: bsc.statistic(min_contig_len)._asdict()
            for min_contig_len in test_contig_cutoffs
        }
    ).T
    assert df.shape == (10, 10)
