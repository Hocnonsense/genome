# -*- coding: utf-8 -*-
"""
* @Date: 2025-01-25 14:18:14
* @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
* @LastEditTime: 2025-07-09 19:26:17
* @FilePath: /genome/genome/clust.py
* @Description:
"""
# """

from typing import Iterable, TypeVar

import numpy as np
import pandas as pd
import scipy.cluster
from scipy.spatial import distance as ssd


class UnionFind:
    def __init__(self, *elements: str, dynamic=False):
        """
        Initialize the union-find data structure with the given elements.
        
        Parameters:
            elements (str): Elements to include in the initial disjoint sets.
            dynamic (bool): If True, allows dynamic addition of new elements during find/union operations.
        """
        self.parent = {e: e for e in elements}
        self.dynamic = dynamic

    def __contains__(self, x):
        """
        Check if the specified element is present in the union-find structure.
        
        Returns:
            bool: True if the element exists, False otherwise.
        """
        return x in self.parent

    def find(self, x):
        """
        Finds and returns the representative element (root) of the set containing `x`, applying path compression for efficiency.
        
        If `dynamic` is enabled and `x` is not present, adds `x` as its own parent.
        """
        if self.dynamic and x not in self:
            return self.parent.setdefault(x, x)
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, x, y):
        """
        Unites the sets containing elements x and y.
        
        If the union-find is not dynamic and either element is missing, the operation is skipped.
        """
        if not self.dynamic:
            if not (x in self and y in self):
                return
        self.parent[self.find(y)] = self.find(x)

    def groups(self):
        """
        Return all groups of connected elements as sets.
        
        Returns:
        	A collection of sets, each containing elements that are connected in the union-find structure.
        """
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
        """
        Create groups of connected elements from a list of pairwise links using the union-find algorithm.
        
        Parameters:
            links: Iterable of pairs representing connections between elements.
            elements: Optional iterable of elements to initialize the union-find structure.
        
        Returns:
            A collection of sets, each set containing elements that are connected.
        """
        uf = cls(*elements)
        uf.dynamic = len(uf.parent) == 0
        for x, y in links:
            uf.union(x, y)
        return uf.groups()


def cluster_hierarchical(
    db: pd.DataFrame, linkage_method="single", linkage_cutoff=0.10
):
    """
    Performs hierarchical clustering on a symmetric distance matrix and assigns cluster labels.
    
    Parameters:
        db (pd.DataFrame): Symmetric distance matrix with genomes as columns and rows.
        linkage_method (str): Linkage method for hierarchical clustering (default is "single").
        linkage_cutoff (float): Distance threshold for forming flat clusters (default is 0.10).
    
    Returns:
        tuple: A tuple containing:
            - pd.DataFrame: DataFrame mapping each genome to its assigned cluster.
            - np.ndarray: Linkage matrix produced by hierarchical clustering.
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
    """
    Resolve cluster label conflicts between two cluster assignments by combining labels for ambiguous cases.
    
    For each pair of cluster labels from `clu_main` and `clu_compare`, if either label is associated with multiple unique pairs, a combined label is generated using `diff_clu_namer`. Otherwise, the primary cluster label from `clu_main` is retained.
    
    Parameters:
        clu_main: Primary cluster labels.
        clu_compare: Secondary cluster labels.
        diff_clu_namer: Function to generate a combined label from two cluster labels (default: concatenation with underscore).
    
    Returns:
        pd.Series: Series of resolved cluster labels, with combined labels for conflicts.
    """
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
    """
    Read and process a FastANI output file, returning ANI and alignment fraction data in long format.
    
    Parameters:
        out_base (str): Path to the FastANI output file.
    
    Returns:
        pd.DataFrame: DataFrame with columns `reference`, `query`, `ani` (as a fraction), and `alignment_fraction` (fraction of aligned sequence).
    """
    odb = pd.read_csv(
        out_base, names=["reference", "query", "ani", "j1", "j2"], sep="\t"
    ).assign(
        alignment_fraction=lambda df: np.where(
            df["j2"] != 0, df["j1"] / df["j2"], np.nan
        ),
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
    """
    Performs pairwise ANI clustering on genome similarity data using hierarchical and optional single-linkage methods.
    
    Filters genome pairs by alignment coverage and computes a distance metric as 1 minus the average ANI. Performs hierarchical clustering on the resulting distance matrix. If `link_single` is True, additionally identifies connected components using single-linkage grouping and resolves cluster assignments between the two methods.
    
    Parameters:
        ndb (pd.DataFrame): DataFrame containing columns 'reference', 'query', 'ani', and 'alignment_fraction'.
        s_lmethod (str, optional): Linkage method for hierarchical clustering (default is "single").
        link_single (bool, optional): Whether to perform additional single-linkage grouping (default is True).
        cov (float, optional): Minimum alignment coverage threshold (default is 0.5).
        ani (float, optional): ANI threshold for clustering (default is 0.99).
    
    Returns:
        tuple: A tuple containing a DataFrame with columns 'genome', 'cluster', and 'cluster_single' (if `link_single` is True), or just 'genome' and 'cluster' (if False), and the hierarchical clustering linkage matrix.
    """
    s_l_cutoff = 1 - ani
    _min_cov = float(cov)
    Ldb = ndb.assign(
        ani=lambda df: df.apply(
            lambda s: (_min_cov < s["alignment_fraction"]) * s["ani"],
            axis=1,
        ),
        av_ani=lambda df: df.merge(
            df[["query", "reference", "ani"]].rename(
                columns={"query": "reference", "reference": "query", "ani": "ani_r"}
            ),
            how="left",
        )[["ani", "ani_r"]].mean(axis=1),
        dist=lambda df: np.where(
            df["query"] == df["reference"], 0, 1 - df["av_ani"].astype(float)
        ),
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
        .assign(_cluster_single=lambda df: df.reset_index().index + 1)
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
