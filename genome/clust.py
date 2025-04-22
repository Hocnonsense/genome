# -*- coding: utf-8 -*-
"""
* @Date: 2025-01-25 14:18:14
* @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
* @LastEditTime: 2025-04-22 11:10:03
* @FilePath: /genome/genome/clust.py
* @Description:
"""
# """

import numpy as np
import pandas as pd
import scipy.cluster
from scipy.spatial import distance as ssd

from typing import Iterable, TypeVar

import numpy as np
import pandas as pd


class UnionFind:
    def __init__(self, *elements: str, dynamic=False):
        self.parent = {e: e for e in elements}
        self.dynamic = dynamic

    def __contains__(self, x):
        return x in self.parent

    def find(self, x):
        if self.dynamic and x not in self:
            return self.parent.setdefault(x, x)
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, x, y):
        if not self.dynamic:
            if not (x in self and y in self):
                return
        self.parent[self.find(y)] = self.find(x)

    def groups(self):
        groups: dict[str, set[str]] = {}
        for e in self.parent:
            groups.setdefault(self.find(e), set()).add(e)
        return groups.values()

    @classmethod
    def igroup(
        cls,
        links: Iterable[tuple[str, str]],
        elements: Iterable[str] = tuple(),
    ):
        uf = cls(*elements)
        uf.dynamic = len(uf.parent) == 0
        for x, y in links:
            uf.union(x, y)
        return uf.groups()


def cluster_hierarchical(
    db: pd.DataFrame, linkage_method="single", linkage_cutoff=0.10
):
    """
    Perform hierarchical clustering on a symmetrical distiance matrix

    Args:
        db: result of db.pivot usually
        linkage_method: passed to scipy.cluster.hierarchy.fcluster
        linkage_cutoff: distance to draw the clustering line (default = .1)

    Returns:
        list: [Cdb, linkage]
    """
    # Generate linkage dataframe
    arr = ssd.squareform(np.asarray(db))
    linkage = scipy.cluster.hierarchy.linkage(arr, method=linkage_method)
    # Form clusters
    fclust = scipy.cluster.hierarchy.fcluster(
        linkage, linkage_cutoff, criterion="distance"
    )
    # Make Cdb
    cdb = pd.DataFrame({"genome": pd.Series(db.columns)})
    cdb["cluster"] = pd.Series(fclust)
    return cdb, linkage


T = TypeVar("T")
K = TypeVar("K")


def resolve_ungroup(
    clu_main: Iterable[T],
    clu_compare: Iterable[K],
    diff_clu_namer=lambda c1, c2: f"{c1}_{c2}",
):
    check = pd.DataFrame({"c1": clu_main, "c2": clu_compare})
    uniq = check.drop_duplicates()
    warn1 = uniq["c1"].value_counts().pipe(lambda s: s[s > 1]).index
    warn2 = uniq["c2"].value_counts().pipe(lambda s: s[s > 1]).index
    warn_index = check["c1"].isin(warn1) | check["c2"].isin(warn2)
    clu0 = check.apply(lambda s: diff_clu_namer(s["c1"], s["c2"]), axis=1)
    clu0[~warn_index] = check["c1"]
    # print(warn1, warn2, clu0, flush=True)
    return clu0


def read_fastani(out_base):
    odb = pd.read_csv(
        out_base, names=["reference", "query", "ani", "j1", "j2"], sep="\t"
    ).assign(
        alignment_fraction=lambda df: [(j1 / j2) for j1, j2 in zip(df["j1"], df["j2"])],
        ani=lambda df: df["ani"].apply(lambda x: x / 100),
    )
    fdb = (
        odb.pivot(index="reference", columns="query", values="ani")
        .reset_index(level=0)
        .fillna(0)
        .melt(id_vars=["reference"])
        .rename(columns={"value": "ani"})
    )
    # Add back alignment coverage
    fodb = fdb.merge(
        odb[["reference", "query", "alignment_fraction"]],
        on=["reference", "query"],
        how="outer",
    ).assign(
        alignment_fraction=lambda df: df["alignment_fraction"].fillna(0),
    )
    return fodb


def run_pairwise_ani(
    ndb: pd.DataFrame, s_lmethod="single", link_single=True, cov=0.5, ani=0.99
):
    s_l_cutoff = 1 - ani
    _min_cov = float(cov)
    Ldb = ndb.assign(
        ani=lambda df: df.apply(
            lambda s: (_min_cov < s["alignment_fraction"]) * s["ani"],
            axis=1,
        ),
        av_ani=lambda df: df.set_index(["query", "reference"])["ani"].pipe(
            lambda s: [r == q or np.mean([s[q, r], s[r, q]]) for q, r in s.index]
        ),
        dist=lambda df: 1 - df["av_ani"].astype(float),
    )
    Gdb, linkage = cluster_hierarchical(
        Ldb.pivot(index="reference", columns="query", values="dist"),
        linkage_method=s_lmethod,
        linkage_cutoff=s_l_cutoff,
    )
    if not link_single:
        return Gdb, linkage
    m2_groups = Ldb.pipe(lambda df: df[df["dist"] <= s_l_cutoff])[
        ["query", "reference"]
    ].pipe(
        lambda df: UnionFind.igroup(df.itertuples(index=False), Ldb["query"].unique())
    )
    mdb = (
        pd.DataFrame({"genome": m2_groups})
        .assign(_cluster_single=lambda df: range(1, len(df["genome"]) + 1))
        .explode("genome")
    )
    cdb2 = (
        Gdb[["genome", "cluster"]]
        .merge(mdb)
        .assign(
            cluster_single=lambda df: resolve_ungroup(
                df["cluster"], df["_cluster_single"]
            )
        )
    )[["genome", "cluster", "cluster_single"]]
    return cdb2, linkage
