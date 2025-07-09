# -*- coding: utf-8 -*-
"""
* @Date: 2022-10-15 21:29:41
* @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
* @LastEditTime: 2025-07-01 16:40:24
* @FilePath: /genome/genome/gene_clust.py
* @Description:
"""


import os
import shutil
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import TYPE_CHECKING, Collection, Iterable, Literal, NamedTuple

import pandas as pd
from Bio import SeqIO, SeqRecord

from .pyrule import rules_dir, smk

PathLike = str | Path

GENE_CLUST_SMK = rules_dir / "gene_clust.smk"


class _CluBase:
    @classmethod
    def in_faa(cls, prefix: PathLike = "{prefix}"):
        return Path(f"{prefix}.faa")

    @classmethod
    def from_in_faa(cls, faa: PathLike):
        "create UniRefClu from input faa file"
        assert f"{faa}".endswith(".faa")
        return cls.from_prefix(str(faa)[:-4])

    @classmethod
    def from_prefix(cls, prefix: PathLike):
        "create UniRefClu from prefix of output files"
        return cls(*(Path(f"{i}".format(prefix=prefix)) for i in cls()))

    @classmethod
    def from_aout(cls, aout: PathLike):
        "auto recognize prefix from output"
        aout_ = str(aout)
        for suffix in (i.name.format(prefix="") for i in cls()):
            if aout_.endswith(suffix):
                return cls.from_prefix(aout_[: -len(suffix)])
        raise KeyError("Cannot determine suffix, please check file name")

    def _modify(self, *modify: PathLike | None):
        return self.__class__(
            *(
                Path(modify) if modify else default
                for default, modify in zip(self, modify)
            )
        )

    TSV_COLS: tuple[str, ...] = ("All",)

    def _load_rep2all(self, *files):
        rep2all: pd.DataFrame
        for i, f in enumerate(files):
            names = self.TSV_COLS[i : i + 2][::-1]
            df = pd.read_csv(f, sep="\t", header=None, names=names)
            rep2all = df if i == 0 else rep2all.merge(df)
        return rep2all

    def load_rep2all(self, keep: Literal[True] = True) -> pd.DataFrame: ...

    @property
    def rep2all(self):
        return self.load_rep2all(keep=True)

    if TYPE_CHECKING:

        def __iter__(self):
            yield Path()

    @classmethod
    def exec_rep2all(
        cls,
        files: Iterable[PathLike] | None = None,
        faas: Iterable[SeqRecord.SeqRecord] | None = None,
        keep_prefix: PathLike | Literal[False] | None = None,
        threads=4,
        profile: PathLike | None = None,
    ):
        assert files or faas
        with NamedTemporaryFile("w", suffix=".faa", delete=True) as tmpf:
            with open(tmpf.name) as fi:
                if faas:
                    SeqIO.write(faas, tmpf.name, "fasta-2line")
                    keep_prefix = keep_prefix or False
                if files:
                    for i, file in enumerate(files):
                        with open(file) as fi:
                            while True:  # read 16 Kib words one time
                                if not (block := fi.read(65536)):
                                    break
                                tmpf.write(block)
                            tmpf.flush()
                    if keep_prefix is None and i == 0:
                        # infer filename automatically if not given in some cases
                        prefix = Path(file).with_suffix("")
                        if str(file) == cls.in_faa(prefix):
                            keep_prefix = False
                            tmpf.name = str(file)
            if keep_prefix:
                Path(keep_prefix).mkdir(parents=True, exist_ok=True)
                shutil.move(tmpf.name, cls.in_faa(keep_prefix))
                tpmf_out = cls.from_prefix(keep_prefix)
            else:
                tpmf_out = cls.from_prefix(tmpf.name[:-4])
            tpmf_out_str = " ".join([str(i) for i in tpmf_out])
            if profile:
                profile_str = f"--profile {profile}"
            smk_params = (
                f"-s {GENE_CLUST_SMK} "
                f"{tpmf_out_str} "
                f"--use-conda "
                f"{profile_str} "
                f"-c{threads} -p "
            )

            try:
                os.system(f"ls {tmpf.name}")
                print("params:", "snakemake", smk_params)
                smk(smk_params)
            except SystemExit as se:
                if se.code:
                    print(se.code, se.with_traceback(None))
                    raise RuntimeError("snakemake seems not run successfully.")
                else:
                    return tpmf_out.rep2all

        raise NotImplementedError("")


class _UniRefClu(NamedTuple):
    u100: Path = Path("{prefix}-uniref100.tsv")
    u90: Path = Path("{prefix}-uniref90.tsv")
    u50: Path = Path("{prefix}-uniref50.tsv")


class UniRefClu(_UniRefClu, _CluBase):
    @classmethod
    def from_prefix_modify(
        cls,
        prefix: PathLike,
        u100: PathLike | None = None,
        u90: PathLike | None = None,
        u50: PathLike | None = None,
    ):
        return cls.from_prefix(prefix)._modify(u100, u90, u50)

    TSV_COLS = "All", "U100", "U90", "U50"

    def load_rep2all(
        self,
        keep: Iterable[Literal["All", "U100", "U90", "U50"]] | Literal[True] = (
            "All",
            "U50",
        ),
    ):
        """
        keep:
            if True, keep all columns,
            else: keep the columns in the list
        """
        keep_ = list(self.TSV_COLS if keep is True else keep)
        return self._load_rep2all(self.u100, self.u90, self.u50)[keep_]


def extract(
    subset: Collection[str],
    files: Iterable[PathLike] | None = None,
    faas: Iterable[SeqRecord.SeqRecord] | None = None,
):
    assert files or faas
    if faas:
        for faa in faas:
            if faa.id in subset:
                yield faa
    if files:
        for file in files:
            for faa in SeqIO.parse(file, "fasta"):
                if faa.id in subset:
                    yield faa


class _MmseqOut(NamedTuple):
    all_100: Path = Path("{prefix}-clu_100.tsv")
    all_clu: Path = Path("{prefix}-clu.tsv")
    all_clu_faa: Path = Path("{prefix}-clu_rep.faa")


class MmseqOut(_MmseqOut, _CluBase):
    @classmethod
    def from_prefix_modify(
        cls,
        prefix: PathLike,
        all_100: PathLike | None = None,
        all_clu: PathLike | None = None,
        all_clu_faa: PathLike | None = None,
    ):
        return cls.from_prefix(prefix)._modify(all_100, all_clu, all_clu_faa)

    TSV_COLS = "All", "Rep100", "Rep"

    def load_rep2all(
        self,
        keep: Iterable[Literal["All", "Rep100", "Rep"]] | Literal[True] = (
            "All",
            "Rep",
        ),
    ):
        keep_ = list(self.TSV_COLS if keep is True else keep)
        return self._load_rep2all(self.all_100, self.all_clu)[keep_]


def mmseq_clust(
    files: Iterable[PathLike] | None = None,
    faas: Iterable[SeqRecord.SeqRecord] | None = None,
    out_prefix: PathLike | MmseqOut = "gene",
    threads=4,
) -> MmseqOut:
    # infer gff_out automatically if not given in some cases
    if not isinstance(out_prefix, tuple):
        _out_prefix = MmseqOut.from_prefix(out_prefix)
    else:
        _out_prefix = MmseqOut(*out_prefix)
    with NamedTemporaryFile("w", suffix=".faa", delete=True) as tmpf:
        MmseqOut.exec_rep2all(files, faas, tmpf.name[:-4], threads)
        tpmf_out = MmseqOut.from_prefix(tmpf.name[:-4])
        os.system(f"ls {tmpf.name}")
        for ffrom, fto in zip(tpmf_out, _out_prefix):
            Path(fto).parent.mkdir(parents=True, exist_ok=True)
            shutil.move(ffrom, fto)
        return _out_prefix


class _MmFamily(NamedTuple):
    mf100: Path = Path("{prefix}-mf100.tsv")
    mfamily: Path = Path("{prefix}-mfamily.tsv")


class MmFamily(_MmFamily, _CluBase):
    @classmethod
    def from_prefix_modify(
        cls,
        prefix: PathLike,
        mf100: PathLike | None = None,
        mfamily: PathLike | None = None,
    ):
        return cls.from_prefix(prefix)._modify(mf100, mfamily)

    TSV_COLS = "All", "F100", "Family"

    def load_rep2all(
        self,
        keep: Iterable[Literal["All", "F100", "Family"]] | Literal[True] = (
            "All",
            "Family",
        ),
    ):
        """
        keep:
            if True, keep all columns,
            else: keep the columns in the list
        """
        keep_ = list(self.TSV_COLS if keep is True else keep)
        return self._load_rep2all(self.mf100, self.mfamily)[keep_]


class _MmSpecies(NamedTuple):
    """https://www.nature.com/articles/s41559-024-02357-0

    All gene sequences were clustered using MMseqs2 (ref. 42) with minimum overlap of 50%,
        minimum identity threshold of 80% and clustering mode 0.
        The rest of the parameters were left as default.
    """

    mf100: Path = Path("{prefix}-mf100.tsv")
    mspecies: Path = Path("{prefix}-mspecies.tsv")


class MmSpecies(_MmSpecies, _CluBase):
    @classmethod
    def from_prefix_modify(
        cls,
        prefix: PathLike,
        mf100: PathLike | None = None,
        mspecies: PathLike | None = None,
    ):
        return cls.from_prefix(prefix)._modify(mf100, mspecies)

    TSV_COLS = "All", "F100", "Species"

    def load_rep2all(
        self,
        keep: Iterable[Literal["All", "F100", "Species"]] | Literal[True] = (
            "All",
            "Species",
        ),
    ):
        """
        keep:
            if True, keep all columns,
            else: keep the columns in the list
        """
        keep_ = list(self.TSV_COLS if keep is True else keep)
        return self._load_rep2all(self.mf100, self.mspecies)[keep_]
