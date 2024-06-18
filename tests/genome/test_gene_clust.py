# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-15 22:24:10
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-06-18 11:13:39
 * @FilePath: /genome/tests/genome/test_gene_clust.py
 * @Description:
"""

import pytest

from genome.gene_clust import mmseq_clust, MmseqOut, UniRefClu, extract
from genome.gff import Parse


from tests import Path, temp_output, test_files, test_temp
from tests.genome._decorator import pytest_mark_resource


def test_mo_from_in_faa():
    faa = "some_gene.faa"
    mo = MmseqOut.from_in_faa(faa)
    assert mo.all_100 == Path("some_gene-clu_100.tsv")
    assert mo.all_clu == Path("some_gene-clu.tsv")
    assert mo.all_clu_faa == Path("some_gene-clu_rep.faa")


def test_mo_from_in_faa_fail():
    faa = "some_gene.fa"
    with pytest.raises(AssertionError):
        mo = MmseqOut.from_in_faa(faa)


def test_mo_from_prefix():
    prefix = "some_gene"
    mo = MmseqOut.from_prefix(prefix)
    assert mo.all_100 == Path("some_gene-clu_100.tsv")
    assert mo.all_clu == Path("some_gene-clu.tsv")
    assert mo.all_clu_faa == Path("some_gene-clu_rep.faa")


def test_mo_from_aout():
    mo = MmseqOut.from_prefix("some_gene")
    for i in mo:
        assert MmseqOut.from_aout(i) == mo


def test_mo_from_aout_fail():
    with pytest.raises(KeyError):
        mo = MmseqOut.from_aout("some_gene-clu.faa")


def test_mo_from_prefix_modify():
    prefix = "some_gene"
    mo = MmseqOut.from_prefix_modify(prefix, all_clu=Path("some_gene_1/clu.tsv"))
    assert mo.all_100 == Path("some_gene-clu_100.tsv")
    assert mo.all_clu == Path("some_gene_1/clu.tsv")
    assert mo.all_clu_faa == Path("some_gene-clu_rep.faa")


@temp_output
def test_mo_load_rep2all(test_temp: Path):
    mo = MmseqOut.from_prefix(test_temp / "some_gene")
    with open(mo.all_100, "w") as f100:
        f100.write("A1\tA1\n" "A1\tA2\n" "A1\tA3\n" "B1\tB1\n" "C1\tC1\n")
    with open(mo.all_clu, "w") as f95:
        f95.write("A1\tA1\n" "C1\tB1\n" "A1\tC1\n")
    rep2all = mo.load_rep2all()
    assert rep2all.to_csv() == ",Rep,All\n0,A1,A1\n1,A1,A2\n2,A1,A3\n3,C1,B1\n4,A1,C1\n"
    assert (
        mo.load_rep2all(keep=True).to_csv()
        == ",Rep,Rep100,All\n0,A1,A1,A1\n1,A1,A1,A2\n2,A1,A1,A3\n3,C1,B1,B1\n4,A1,C1,C1\n"
    )
    assert (
        mo.load_rep2all(keep=["Rep100", "Rep"]).to_csv()
        == ",Rep100,Rep\n0,A1,A1\n1,A1,A1\n2,A1,A1\n3,B1,C1\n4,C1,A1\n"
    )


@pytest_mark_resource
@temp_output
def test_mmseq_clust(test_temp: Path):
    expect_prefix = test_temp / "some_gene"
    gffs = [
        test_files / "binny_contigs_4bins-top10-prodigal.gvmeta.gff",
    ]
    faa_iter = (faa for gff in gffs for faa in Parse(gff).extract())
    mo = mmseq_clust([], faa_iter, expect_prefix, threads=10)
    assert mo.all_100 == Path("some_gene-clu_100.tsv")
    assert mo.all_clu == Path("some_gene-clu.tsv")
    assert mo.all_clu_faa == Path("some_gene-clu_rep.faa")


@pytest_mark_resource
@temp_output
def test_uniref_clust(test_temp: Path):
    gffs = [
        test_files / "binny_contigs_4bins-top10-prodigal.gvmeta.gff",
    ]
    faa_iter = (faa for gff in gffs for faa in Parse(gff).extract())
    urc = UniRefClu.exec_rep2all([], faa_iter, threads=10)
    faa_iter = (faa for gff in gffs for faa in Parse(gff).extract())
    assert list(extract(urc["U50"].values, [], faa_iter))
