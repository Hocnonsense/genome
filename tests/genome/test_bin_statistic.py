# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-12 19:53:55
* @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
* @LastEditTime: 2025-07-09 17:12:16
* @FilePath: /genome/tests/genome/test_bin_statistic.py
 * @Description:
__file__ = "test/genome/test_bin_statistic.py"
"""

from timeit import timeit

from genome.bin_statistic import BinStatisticContainer, Contig2Bin
from tests import Path, temp_output, test_files, test_temp


@temp_output
def test_contig2bin(test_temp: Path):
    """
    Tests the Contig2Bin class by generating bin files from a contig-to-bin mapping and a FASTA file.
    
    This function verifies that Contig2Bin correctly processes a TSV mapping file and a corresponding FASTA file, producing bin files in a temporary output directory.
    """
    test_fa = test_files / "binny_contigs_4bins.fa"
    test_c2b = test_files / "binny_unitem_unanimous.tsv"
    temp_bin = test_temp / "binny_unitem_unanimous"

    Contig2Bin(test_c2b, test_fa)(temp_bin)


def test_prokka_gff_bin_statistic():
    """
    Test that bin statistics computed from a Prokka GFF annotation file match expected values.
    
    Validates GC content, GC standard deviation, base pair size, maximum contig length, number of contigs, N50 contig length, ambiguous base count, coding density, and gene count for a known test GFF file.
    """
    gff = test_files / "binny_contigs_4bins-top10-prodigal.gvmeta.gff"
    bsc = BinStatisticContainer.read_gff(gff, min_aa_len=0)
    bs = bsc.statistic()
    bsc.calculate_gc_std()
    assert int(bs.gc * 100) == 58
    assert round(bs.gc_std * 10000) == 235
    assert bs.bp_size == 101874
    assert bs.max_contig_len == 31552
    assert bs.contigs_num == 10
    assert bs.contig_n50 == 21642
    assert bs.ambiguous_bases_num == 0
    assert int(bs.coding_density * 100) == 86
    assert bs.genes_num == 105


def test_genome_bin_statistic():
    """
    Test computation of genome bin statistics from FASTA and GFF files.
    
    This test verifies that `BinStatisticContainer` correctly computes various statistics for genome bins, including GC content, GC standard deviation, total base pairs, maximum contig length, number of contigs, N50, ambiguous base count, coding density, and gene count. It checks results for multiple minimum contig length cutoffs and validates output formats for both DataFrame and dictionary representations.
    """
    genome = test_files / "binny_contigs_4bins.fa"
    bsc = BinStatisticContainer.read_contig(genome, "fasta")
    test_contig_cutoffs = (0, 500, 1000, 1500, 2000, 2500, 3000, 4000, 5000, 10000)
    df = BinStatisticContainer.to_data_frame(
        {
            min_contig_len: bsc.statistic(min_contig_len)
            for min_contig_len in test_contig_cutoffs
        }
    )
    assert df.to_csv() == (
        ",gc,gc_std,bp_size,max_contig_len,contigs_num,contig_n50,ambiguous_bases_num,contig_cutoff,coding_density,genes_num\n"
        "0,0.5729274134675137,0.04675647608437309,9281832.0,2862522.0,412.0,50324.0,0.0,0.0,0.0,0.0\n"
        "500,0.5729274134675137,0.04675647608437309,9281832.0,2862522.0,412.0,50324.0,0.0,500.0,0.0,0.0\n"
        "1000,0.5729580565857126,0.04677822669500688,9272588.0,2862522.0,400.0,50324.0,0.0,1000.0,0.0,0.0\n"
        "1500,0.5731026378071076,0.046286125146007766,9240289.0,2862522.0,374.0,50888.0,0.0,1500.0,0.0,0.0\n"
        "2000,0.5732389810676235,0.04619816129069593,9210307.0,2862522.0,357.0,50888.0,0.0,2000.0,0.0,0.0\n"
        "2500,0.5733409290071567,0.046714789166143564,9173643.0,2862522.0,341.0,50888.0,0.0,2500.0,0.0,0.0\n"
        "3000,0.5734988441163675,0.04661651188392546,9143654.0,2862522.0,330.0,51525.0,0.0,3000.0,0.0,0.0\n"
        "4000,0.5738066237265611,0.04671084550171524,9081021.0,2862522.0,312.0,51525.0,0.0,4000.0,0.0,0.0\n"
        "5000,0.5744143866969582,0.04682288888739633,8957831.0,2862522.0,285.0,52361.0,0.0,5000.0,0.0,0.0\n"
        "10000,0.5777346851059988,0.04906833937103942,8229293.0,2862522.0,187.0,54972.0,0.0,10000.0,0.0,0.0\n"
    )
    gff = test_files / "binny_contigs_4bins-top10-prodigal.gvmeta.gff"
    bsc = BinStatisticContainer.read_gff(gff)
    assert bsc.statistic()._asdict() == {
        "gc": 0.58939474252508,
        "gc_std": 0.023501846757011037,
        "bp_size": 101874,
        "max_contig_len": 31552,
        "contigs_num": 10,
        "contig_n50": 21642,
        "ambiguous_bases_num": 0,
        "contig_cutoff": 0,
        "coding_density": 0.8597581325951665,
        "genes_num": 104,
    }


def test_genome_stat_speed(report=False):
    """
    Test that quick and standard contig reading methods produce identical sequence statistics.
    
    Optionally reports the execution time for both methods over multiple iterations.
    """
    genome = test_files / "binny_contigs_4bins.fa"
    bsc_seq_quick = BinStatisticContainer.quick_read_contig(
        genome, "fasta"
    ).calculate_seq_stats()
    bsc_seq = BinStatisticContainer.read_contig(genome, "fasta").calculate_seq_stats()
    assert bsc_seq_quick == bsc_seq
    if report:
        quick_time = timeit(
            lambda: BinStatisticContainer.quick_read_contig(
                genome, "fasta"
            ).calculate_seq_stats(),
            number=20,
        )
        normal_time = timeit(
            lambda: BinStatisticContainer.read_contig(
                genome, "fasta"
            ).calculate_seq_stats(),
            number=20,
        )
        print(f"if only calculate seq length, will spend {quick_time:.4f} seconds")
        print(f"if calculate more features, will spend {normal_time:.4f} seconds")
