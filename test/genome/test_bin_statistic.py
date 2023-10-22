# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-12 19:53:55
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-09-17 18:16:16
 * @FilePath: /genome/test/genome/test_bin_statistic.py
 * @Description:
__file__ = "test/genome/test_bin_statistic.py"
"""

from pathlib import Path
from timeit import timeit

import pandas as pd

from genome.bin_statistic import BinStatisticContainer, contig2bin

test_temp = Path(__file__).parent.parent / "temp"
test_files = Path(__file__).parent.parent / "file"


def test_prokka_gff_bin_statistic():
    gff = test_files / "metadecoder.1-prokka.Bacteria.gff"
    bsc = BinStatisticContainer.read_gff(gff)
    bs = bsc.statistic()
    bsc.calculate_gc_std()
    assert int(bs.gc * 100) == 56
    assert round(bs.gc_std * 10000) == 206
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


def test_genome_stat_speed():
    genome = test_files / "metadecoder.4.fa"
    bsc_seq_quick = BinStatisticContainer.quick_read_contig(
        genome, "fasta"
    ).calculate_seq_stats()
    bsc_seq = BinStatisticContainer.read_contig(genome, "fasta").calculate_seq_stats()
    assert bsc_seq_quick == bsc_seq
    quick_time = timeit(
        lambda: BinStatisticContainer.quick_read_contig(
            genome, "fasta"
        ).calculate_seq_stats(),
        number=100,
    )
    normal_time = timeit(
        lambda: BinStatisticContainer.read_contig(
            genome, "fasta"
        ).calculate_seq_stats(),
        number=100,
    )
    print(f"if only calculate seq length, will spend {quick_time:.4f} seconds")
    print(f"if calculate more featuers, will spend {normal_time:.4f} seconds")


def test_contig2bin():
    test_fa = test_files / "02_assem..TY.041_cut.fa"
    test_c2b = test_files / "binsingle" / "dastool.tsv"
    temp_bin = test_temp / "dastool"

    contig2bin(temp_bin, test_c2b, test_fa)
