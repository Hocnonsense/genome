# -*- coding: utf-8 -*-
"""
* @Date: 2022-04-15 13:56:44
* @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
* @LastEditTime: 2025-05-13 22:29:55
* @FilePath: /genome/genome/gene_annot.py
* @Description:
"""

import os
import re
from pathlib import Path
from typing import Generator, Iterable, NamedTuple, Optional, Union

import pandas as pd

try:
    from PyLib.PyLibTool.file_info import verbose_import

    logger = verbose_import(__name__, __doc__)
except ImportError:
    import logging

    logger = logging.getLogger(__name__)


def read_table(text, sep="\t", annot="#", title: Optional[list] = None, openit=False):
    if openit:
        text = open(text)  # type: ignore
    for line in text:
        if line.startswith(annot):
            if title is not None:
                title.clear()
                title.extend(line[len(annot) :].rstrip().split(sep))
            continue
        values = line.strip().split(sep)
        if values:
            yield values
    if openit:
        text.close()  # type: ignore


PathLike = Union[str, Path]


class Gene2KOOut(NamedTuple):
    faa: Path = Path("gene-clu_rep.faa")
    ghost: Path = Path("gene-clu_rep-ghost.faa")
    kofam: Path = Path("gene-clu_rep-kofam.faa")
    eggnog: Path = Path("gene-clu_rep-eggnog.faa")
    mantis: Path = Path("gene-clu_rep-mantis.faa")

    @classmethod
    def from_prefix(cls, prefix: PathLike):
        return cls(
            faa=Path(f"{prefix}.faa"),
            ghost=Path(f"{prefix}-ghost.faa"),
            kofam=Path(f"{prefix}-kofam.faa"),
            eggnog=Path(f"{prefix}-eggnog.faa"),
            mantis=Path(f"{prefix}-mantis.faa"),
        )

    @classmethod
    def from_faa(cls, faa: PathLike):
        assert str(faa).endswith(".faa")
        prefix = str(faa)[:-4]
        return cls(
            faa=Path(f"{prefix}.faa"),
            ghost=Path(f"{prefix}-ghost.faa"),
            kofam=Path(f"{prefix}-kofam.faa"),
            eggnog=Path(f"{prefix}-eggnog.faa"),
            mantis=Path(f"{prefix}-mantis.faa"),
        )


class Gene2KO:
    class Annoters(type):
        """"""

        iters: list[type["Gene2KO.GeneKoIter"]] = None  # type: ignore

        def __new__(cls, name, bases, attrs):
            new_cls = type.__new__(cls, name, bases, attrs)
            if cls.iters is None:
                cls.iters = []
            else:
                cls.iters.append(new_cls)  # type: ignore
            return new_cls

    class GeneKoIter(metaclass=Annoters):
        def __init__(self, filename: Path):
            self.filename = filename

        def __call__(self) -> Generator[tuple[str, str], None, None]:
            raise NotImplementedError

    class ghost(GeneKoIter):
        def __call__(self) -> Generator[tuple[str, str], None, None]:
            with open(self.filename) as text:
                for values in read_table(text):
                    if len(values) > 1:
                        gene, ko = values[0:2]
                        yield gene, ko

    class kofam(GeneKoIter):
        def __call__(self) -> Generator[tuple[str, str], None, None]:
            with open(self.filename) as text:
                for values in read_table(text):
                    if len(values) > 1:
                        gene, ko = values[1:3]
                        yield gene, ko

    class eggnog(GeneKoIter):
        def __call__(self) -> Generator[tuple[str, str], None, None]:
            i_KEGG_ko = 11
            with open(self.filename) as text:
                for values in read_table(text):
                    kos = values[i_KEGG_ko]
                    if len(kos) > 1:
                        gene = values[0]
                        for ko in kos.split(","):
                            yield gene, ko[3:]

    class mantis(GeneKoIter):
        # KO_PATTERN = re.compile("\\b(K\\d{5})(?=(;|$))")
        KO_PATTERN = re.compile("\\b(K\\d{5})\\b")

        def __call__(self) -> Generator[tuple[str, str], None, None]:
            with open(self.filename) as text:
                for values in read_table(text):
                    Ref_Hits = values[2]
                    kos: list[str] = re.findall(self.KO_PATTERN, Ref_Hits)
                    if len(kos) > 0:
                        gene = values[0]
                        for ko in kos:
                            yield gene, ko

    ## collect gene KO
    def get_gene_KOs(self) -> dict[str, str]:
        """Only keep the first match:
        >>> gene_KOs.setdefault(gene, ko)"""
        gene_KOs: dict[str, str] = {}
        for annoter, file in zip(self.Annoters.iters, self.ann_files):
            if not os.path.isfile(file):
                continue

            gene_KOs_: dict[str, list[str]] = {}
            for gene, ko in annoter(file)():
                gene_KOs_.setdefault(gene, []).append(ko)
            for gene in gene_KOs_:
                if gene not in gene_KOs:
                    gene_KOs[gene] = ":".join(sorted(gene_KOs_[gene]))

        return gene_KOs

    def get_gene_annots(self):
        return pd.Series(self.get_gene_KOs(), name="KO")

    def __init__(self, pattern: Union[Path, list[Path]]):
        if not isinstance(pattern, list):
            pattern = [pattern]
        ann_files_ = self._infer_ann_files(pattern)
        ann_files = [
            ann_files_.get(i, Path()) for i, _ in enumerate(self.Annoters.iters)
        ]

        if not any(ann_files):
            raise FileNotFoundError(f"pattern(s) '{pattern}' do not match any file!")

        self.ann_files = ann_files

    def _infer_ann_files(self, patterns: list[Path]):
        ann_files = {i: Path() for i, _ in enumerate(self.Annoters.iters)}

        for pattern in reversed(patterns):
            pattern_re = re.compile(pattern.name)
            for file in pattern.parent.iterdir():
                if pattern_re.search(file.name):
                    for i, source in enumerate(self.Annoters.iters):
                        if source.__name__ in file.name.lower():
                            logger.warning(
                                f"Detect {source.__name__} annotation: {file}"
                            )
                            ann_files[i] = file
        return ann_files


def get_all_gene_annots(gene_annots: pd.Series, rep2all: Optional[pd.Series] = None):
    ko_exploded = pd.DataFrame(
        gene_annots.apply(lambda x: x.split(":") if x.startswith("K") else [])
        .explode()
        .dropna()
        .rename("KO")
    )

    if rep2all is not None:
        all_gene_annots = ko_exploded.merge(
            rep2all, left_index=True, right_index=True
        ).set_index(rep2all.name)
    else:
        all_gene_annots = ko_exploded
    return all_gene_annots


def main(
    annot_prefix,
    all_100_path: Path,
    all_clu_path: Path,
    all_gene_annots_path: Optional[Path] = None,
):
    from .gene_clust import MmseqOut

    rep2all = MmseqOut(all_100_path, all_clu_path, Path()).load_rep2all()

    # gene_annots = pd.read_csv(gene_annots_path, index_col=0)["ko"]
    gene_annots = Gene2KO(annot_prefix).get_gene_annots()

    all_gene_annots = get_all_gene_annots(gene_annots, rep2all)
    if all_gene_annots_path:
        all_gene_annots.to_csv(all_gene_annots_path)

    return all_gene_annots


class MantisAnnot(NamedTuple):
    query: str
    ref_files: list[str]
    annots: dict[str, list[str]]

    @staticmethod
    def parse_links(links: Iterable[tuple[str, str]]):
        all_annots: dict[str, list[str]] = {}
        for k, v in links:
            all_annots.setdefault(k, []).append(v)
        return all_annots

    @classmethod
    def from_links(cls, line: str):
        (query, ref_files, ref_hits, consensus_hits, total_hits, _links) = (
            line.strip().split("\t", 5)
        )
        return cls(
            query,
            ref_files.split(";"),
            cls.parse_links(i.split(":", 1) for i in _links[2:].split("\t")),  # type: ignore[misc]
        )

    @classmethod
    def from_refs(cls, line: str):
        (query, ref_files, ref_hits, consensus_hits, total_hits, _links) = (
            line.strip().split("\t", 5)
        )
        return cls(
            query,
            ref_files.split(";"),
            cls.parse_links(zip(ref_files.split(";"), ref_hits.split(";"))),
        )

    @classmethod
    def get_link_annots(cls, mantis_file: Path):
        """
        >>> mantis_file = Path("results/genes-mantis.tsv")
        """
        total_annots: dict[str, dict[str, list[str]]] = {}
        with open(mantis_file) as fi:
            next(fi)
            for line in fi:
                if not line:
                    continue
                m_a = cls.from_links(line)
                query, annots = m_a.query, m_a.annots
                total_annots[query] = {
                    "KO": annots.pop("kegg_ko", []),
                    "COG": annots.pop("cog", []),
                    "arCOG": annots.pop("arcog", []),
                    # "description": annots.pop("description", []),
                    "PFAM": annots.pop("pfam", []),
                    "EC": annots.pop("enzyme_ec", []),
                    **annots,
                }
        return pd.DataFrame(total_annots).T

    @classmethod
    def get_ref_annots(cls, mantis_file: Path):
        """
        >>> mantis_file = Path("results/genes-mantis.tsv")
        """
        total_annots: dict[str, dict[str, list[str]]] = {}
        with open(mantis_file) as fi:
            next(fi)
            for line in fi:
                if not line:
                    continue
                m_a = cls.from_refs(line)
                query, annots = m_a.query, m_a.annots
                total_annots[query] = annots
        return pd.DataFrame(total_annots).sort_index().T
