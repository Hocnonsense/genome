# -*- coding: utf-8 -*-
"""
* @Date: 2025-01-25 14:23:55
* @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
* @LastEditTime: 2025-07-09 19:49:30
* @FilePath: /genome/tests/genome/test_clust.py
* @Description:
"""
# """

from tests import Path, temp_output, test_files, test_temp
from genome.clust import UnionFind, read_fastani, resolve_ungroup, run_pairwise_ani


def test_resolve_ungroup():
    """
    Test the resolve_ungroup function for correct handling of group label conflicts.
    
    Verifies that when group labels are consistent, the original groups remain unchanged, and when labels differ, new group identifiers are assigned to resolve conflicts.
    """
    assert list(
        resolve_ungroup(
            [1, 2, 3, 4, 4, 4],
            [1, 2, "b", "a", "a", "a"],
        )
    ) == [1, 2, 3, 4, 4, 4]
    assert list(
        resolve_ungroup(
            [1, 2, 3, 4, 4, 5],
            [2, 1, "b", "a", "a", "a"],
        )
    ) == [1, 2, 3, "4_a", "4_a", "5_a"]


def test_igroup_includes_isolated():
    """
    Test that UnionFind.igroup correctly includes isolated elements as singleton clusters.
    
    Verifies that elements without any links are returned as individual clusters alongside connected groups.
    """
    elements = ["a", "b", "c"]
    links = [("a", "b")]  # "c" has no connections
    clusters = {frozenset(group) for group in UnionFind.igroup(links, elements)}
    assert clusters == {frozenset({"a", "b"}), frozenset({"c"})}


def test_unionfind():
    """
    Test that UnionFind.igroup correctly clusters connected elements into groups.
    
    Verifies that elements linked directly or indirectly are grouped together, while separate clusters remain distinct.
    """
    assert set(
        frozenset(i) for i in UnionFind.igroup((("a", "b"), ("b", "c"), ("d", "e")))
    ) == {frozenset({"a", "b", "c"}), frozenset({"d", "e"})}


@temp_output
def test_run_pairwise_ani(test_temp: Path):
    """
    Test the pairwise ANI clustering workflow using a temporary FASTANI output file.
    
    Creates a temporary ANI data file, reads it, performs clustering with specified parameters, and asserts that the resulting cluster assignments and output columns match expected values.
    """
    with open(test_temp / "ani.txt", "w") as f:
        ani_txt = (
            "E_faecalis_TX0104.fa\tE_faecalis_YI6-1.fna\t98.4653\t889\t1028\n"
            "E_faecalis_TX0104.fa\tE_faecalis_TX0104.fa\t99.9985\t1019\t1028\n"
            "E_faecalis_TX0104.fa\tE_faecalis_T2.fna\t98.3323\t909\t1028\n"
            "E_faecalis_YI6-1.fna\tE_faecalis_YI6-1.fna\t100\t1004\t1007\n"
            "E_faecalis_YI6-1.fna\tE_faecalis_T2.fna\t99.9041\t981\t1007\n"
            "E_faecalis_YI6-1.fna\tE_faecalis_TX0104.fa\t98.3811\t885\t1007\n"
            "E_faecalis_T2.fna\tE_faecalis_T2.fna\t99.9998\t1073\t1080\n"
            "E_faecalis_T2.fna\tE_faecalis_YI6-1.fna\t99.9152\t979\t1080\n"
            "E_faecalis_T2.fna\tE_faecalis_TX0104.fa\t98.2885\t920\t1080\n"
        )
        f.write(ani_txt)
    ndb = read_fastani(f.name)
    assert list(ndb.columns) == ["reference", "query", "ani", "alignment_fraction"]
    cdb = run_pairwise_ani(ndb, "average", cov=0.5, ani=0.99)[0]

    s = cdb.to_dict()
    assert s == {
        "genome": {
            0: "E_faecalis_T2.fna",
            1: "E_faecalis_TX0104.fa",
            2: "E_faecalis_YI6-1.fna",
        },
        "cluster": {0: 1, 1: 2, 2: 1},
        "cluster_single": {0: 1, 1: 2, 2: 1},
    }
