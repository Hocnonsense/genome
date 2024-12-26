# -*- coding: utf-8 -*-
"""
 * @Date: 2024-12-26 10:26:38
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-12-26 21:30:48
 * @FilePath: /genome/tests/genome/test_gene_statistic.py
 * @Description:
"""
# """


from genome.gene_statistic import aa_mw, ARSC, CodonTable, GeneStatisticContainer

from Bio.Seq import Seq

from tests import Path, temp_output, test_files, test_temp


def test_aa_mw():
    aa_seq = Seq("MATRKGFEPSTSGVTGRRSNQLNYLAEFMVGTTGLEPVTLCL*")
    assert round(aa_mw(aa_seq), 4) == 107.9561


def test_arsc():
    aa_seq = Seq("MATRKGFEPSTSGVTGRRSNQLNYLAEFMVGTTGLEPVTLCL*")
    arsc = ARSC.parse(aa_seq)
    assert arsc == (112, 13, 3, arsc.mw, 42)
    assert round(arsc.mw, 4) == 4534.1558
    assert arsc.scale()[:3] == tuple(i / 42 for i in (112, 13, 3))


def test_scu():
    na1 = Seq(
        "ATGGCGACCAGGAAGGGGTTCGAACCCTCGACCTCCGGCGTGACAGGCCGGCGTTCTAACCAGCTGAACTACCTGGCCGAATTTATGGTGGGAACAACAGGGCTCGAACCTGTGACCCTCTGCTTGTAA"
    )
    ctp = CodonTable.get().parse(na1)
    assert [round(ctp.wighted_F(i), 4) for i in range(1, 5)] == [1, 1.7053, 3, 3.2391]
    assert round(ctp.SCU, 4) == 51.3764
    assert round(ctp.gc_variability, 4) == 0.2313


def test_gene_stat():
    gff = test_files / "binny_contigs_4bins-top10-prodigal.gvmeta.gff"
    gsc = GeneStatisticContainer.read_gff(gff)
    gs = gsc.statistic()
    assert round(gs.scu, 4) == 41.1858
    assert round(gs.gc_variability, 4) == 0.2791
    assert round(gs.C_ARSC, 4) == 2.8870
    assert round(gs.N_ARSC, 4) == 0.3130
    assert round(gs.S_ARSC, 4) == 0.0496
    assert round(gs.avg_protein_mw, 4) == 30774.3951
    assert round(gs.avg_protien_len, 4) == 280.1923
    assert round(gs.genes_num, 4) == 104
