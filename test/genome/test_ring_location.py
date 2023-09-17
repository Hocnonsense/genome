# -*- coding: utf-8 -*-
"""
 * @Date: 2023-09-15 18:39:44
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-09-17 17:49:37
 * @FilePath: /genome/test/genome/test_ring_location.py
 * @Description:
"""

from Bio import Seq
from genome.ring_location import RingSequenceData


def test_rsd_getitem():
    s1 = Seq.Seq(RingSequenceData("AAAAAAATCG"))
    s2 = Seq.Seq("AAAAAAATCG")

    assert s1[::-1] == s2[::-1]
    assert s1[1:10] == s2[1:10]
    assert s1[2:-1] == s2[2:-1]
    assert not s1[2 : -1 - len(s1)]
    assert s1[6:3] == s2[6:] + s2[:3]
    assert s1[-1 : -1 + len(s1)] == s2[-1:] + s2[: -1 + len(s2)]
    assert s1[2 : -1 + len(s1) * 0] == s1[2 : -1 + len(s1) * 1]
    assert s1[2 : -1 + len(s1) * 2] == s2[2:-1] + s1[-1 : -1 + len(s1)] * 1
    assert s1[2 : -1 + len(s1) * 3] == s2[2:-1] + s1[-1 : -1 + len(s1)] * 2
    assert s1[2 - len(s1) : 9 - len(s1)] == s2[2:9]
    assert s1[6:9] == s2[6:9]
    assert s1[6 : 9 + len(s1) * 0] == s2[6:9]
    assert s1[6 : 9 + len(s1) * 1] == s2[6:9] + s1[9 : 9 + len(s1)] * 1
    assert s1[6 : 9 + len(s1) * 2] == s2[6:9] + s1[9 : 9 + len(s1)] * 2
    assert s1[6 : 9 + len(s1) * 3] == s2[6:9] + s1[9 : 9 + len(s1)] * 3
    assert s1[6 - len(s1) : 9 - len(s1)] == s2[6:9]


def test_rsd_contains():
    s1 = Seq.Seq(RingSequenceData("AAAAAAATCG"))

    assert s1[2 : 2 + len(s1)] in s1
