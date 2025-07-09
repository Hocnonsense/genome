# -*- coding: utf-8 -*-
"""
* @Date: 2024-12-26 10:26:38
* @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
* @LastEditTime: 2025-07-09 19:56:00
* @FilePath: /genome/tests/genome/test_gene_statistic.py
* @Description:
"""
# """


from genome.gene_statistic import aa_mw, ARSC, CodonTable, GeneStatisticContainer
from genome.gff import parse

from Bio.Seq import Seq

from tests import Path, temp_output, test_files, test_temp


def test_aa_mw():
    aa_seq = Seq("MATRKGFEPSTSGVTGRRSNQLNYLAEFMVGTTGLEPVTLCL*")
    assert round(aa_mw(aa_seq), 4) == 107.9561
    arsc = ARSC.parse(aa_seq)
    assert aa_mw(aa_seq) == arsc.mw / arsc.len


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
    """
    in GCA_019978365.1.gff, there is a transl_except of Sec. use this case to test gene_statistic
    """
    gff = parse(test_files / "GCA_019978365.1.gff", test_files / "GCA_019978365.1.fa")
    # (recds,) = [i for i in gff.extract(translate=False) if i.id == "BCX53216.1"]
    # (recaa,) = [i for i in gff.extract() if i.id == "BCX53216.1"]
    gsc = GeneStatisticContainer.read_gff_parser(gff)
    gs = gsc.statistic()
    assert gs.table == "11"
    print(gs)
    assert round(gs.scu, 4) == 41.7683
    assert round(gs.gc_variability, 4) == 0.2562
    assert round(gs.C_ARSC, 4) == 2.8389
    assert round(gs.N_ARSC, 4) == 0.3635
    assert round(gs.S_ARSC, 4) == 0.0369
    assert round(gs.avg_protein_mw, 4) == 34988.4176
    assert round(gs.avg_protein_len, 4) == 321.4920
    assert round(gs.genes_num, 4) == 4996
