# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-15 21:29:41
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-06-18 11:30:56
 * @FilePath: /genome/genome/gene_clust.py
 * @Description:
"""


import os
import shutil
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Collection, Iterable, Literal, NamedTuple

import pandas as pd
from Bio import SeqIO, SeqRecord

from .pyrule import smk, smk_conda_env, rules_dir

PathLike = str | Path


class UniRefClu(NamedTuple):
    u100: str = "{prefix}-uniref100.tsv"
    u90: str = "{prefix}-uniref90.tsv"
    u50: str = "{prefix}-uniref50.tsv"

    @classmethod
    def in_faa(cls, prefix: PathLike = "{prefix}"):
        return f"{prefix}.faa"

    @classmethod
    def from_in_faa(cls, faa: PathLike):
        "create UniRefClu from input faa file"
        assert str(faa).endswith(".faa")
        return cls.from_prefix(str(faa)[:-4])

    @classmethod
    def from_prefix(cls, prefix: PathLike):
        "create UniRefClu from prefix of output files"
        return cls(*(i.format(prefix=prefix) for i in cls()))

    @classmethod
    def from_aout(cls, aout: PathLike):
        "auto recognize prefix from output"
        aout_ = str(aout)
        for suffix in (i.format(prefix="") for i in cls()):
            if aout_.endswith(suffix):
                return cls.from_prefix(aout_[: -len(suffix)])
        raise KeyError("Cannot determine suffix, please check file name")

    @classmethod
    def from_prefix_modify(
        cls,
        prefix: PathLike,
        u100: PathLike | None = None,
        u90: PathLike | None = None,
        u50: PathLike | None = None,
    ):
        return cls(
            *(
                modify if modify else default
                for default, modify in zip(cls.from_prefix(prefix), (u100, u90, u50))
            )
        )

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
        if keep is True:
            keep_ = ["All", "U100", "U90", "U50"]
        else:
            keep_ = list(keep)
        u100 = pd.read_csv(self.u100, sep="\t", header=None, names=["U100", "All"])
        u90 = pd.read_csv(self.u90, sep="\t", header=None, names=["U90", "U100"])
        u50 = pd.read_csv(self.u50, sep="\t", header=None, names=["U50", "U90"])
        rep2all = u50.merge(u90).merge(u100)[keep_]
        return rep2all

    @classmethod
    def exec_rep2all(
        cls,
        files: Iterable[PathLike] | None = None,
        faas: Iterable[SeqRecord.SeqRecord] | None = None,
        threads=4,
    ):
        # infer gff_out automatically if not given in some cases
        assert files or faas
        with NamedTemporaryFile("w", suffix=".faa", delete=True) as tmpf:
            tpmf_out = cls.from_prefix(tmpf.name[:-4])
            tpmf_out_str = " ".join([str(i) for i in tpmf_out])

            with open(tmpf.name) as fi:
                if faas:
                    SeqIO.write(faas, tmpf.name, "fasta-2line")
                if files:
                    for file in files:
                        with open(file) as fi:
                            while True:
                                # read 16 Kib words one time
                                block = fi.read(65536)
                                if not block:
                                    break
                                tmpf.write(block)
                            tmpf.flush()

            target_smk_file = rules_dir / "gene_clust.smk"
            smk_params = (
                f"-s {target_smk_file} "
                f"{tpmf_out_str} "
                f"--use-conda "
                f"--conda-prefix {smk_conda_env} "
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
                    return tpmf_out.load_rep2all(keep=True)

        raise NotImplementedError("")


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


class MmseqOut(NamedTuple):
    all_100: Path = Path("gene-clu_100.tsv")
    all_clu: Path = Path("gene-clu.tsv")
    all_clu_faa: Path = Path("gene-clu_rep.faa")

    @classmethod
    def from_in_faa(cls, faa: PathLike):
        "create MmseqOut from input faa file"
        assert str(faa).endswith(".faa")
        return cls.from_prefix(str(faa)[:-4])

    @classmethod
    def from_prefix(cls, prefix: PathLike):
        "create MmseqOut from prefix of output files"
        return cls(
            Path(f"{prefix}-clu_100.tsv"),
            Path(f"{prefix}-clu.tsv"),
            Path(f"{prefix}-clu_rep.faa"),
        )

    @classmethod
    def from_aout(cls, aout: PathLike):
        "auto recognize prefix from output"
        aout_ = str(aout)
        for suffix in ("-clu_100.tsv", "-clu.tsv", "-clu_rep.faa"):
            if aout_.endswith(suffix):
                return cls.from_prefix(aout_[: -len(suffix)])
        raise KeyError("Cannot determine suffix, please check file name")

    @classmethod
    def from_prefix_modify(
        cls,
        prefix: PathLike,
        all_100: PathLike | None = None,
        all_clu: PathLike | None = None,
        all_clu_faa: PathLike | None = None,
    ):
        all_100_ = all_100 if all_100 else f"{prefix}-clu_100.tsv"
        all_clu_ = all_clu if all_clu else f"{prefix}-clu.tsv"
        all_clu_faa_ = all_clu_faa if all_clu_faa else f"{prefix}-clu_rep.faa"

        return cls(Path(all_100_), Path(all_clu_), Path(all_clu_faa_))

    def load_rep2all(
        self,
        keep: Iterable[Literal["Rep", "Rep100", "All"]] | Literal[True] = (
            "Rep",
            "All",
        ),
    ):
        """
        keep:
            if True, keep all columns,
            else: keep the columns in the list
        """
        if keep is True:
            keep_ = ["Rep", "Rep100", "All"]
        else:
            keep_ = list(keep)
        all_100 = pd.read_csv(
            self.all_100, sep="\t", header=None, names=["Rep100", "All"]
        )
        all_clu = pd.read_csv(
            self.all_clu, sep="\t", header=None, names=["Rep", "Rep100"]
        )
        rep2all = all_100.merge(all_clu)[keep_]
        return rep2all


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
    assert files or faas
    with NamedTemporaryFile("w", suffix=".faa", delete=True) as tmpf:
        tpmf_out = MmseqOut.from_prefix(tmpf.name[:-4])
        tpmf_out_str = " ".join([str(i) for i in tpmf_out])

        with open(tmpf.name) as fi:
            if faas:
                SeqIO.write(faas, tmpf.name, "fasta-2line")
            if files:
                for file in files:
                    with open(file) as fi:
                        while True:
                            # read 16 Kib words one time
                            block = fi.read(65536)
                            if not block:
                                break
                            tmpf.write(block)
                        tmpf.flush()

        target_smk_file = rules_dir / "gene_clust.smk"
        smk_params = (
            f"-s {target_smk_file} "
            f"{tpmf_out_str} "
            f"--use-conda "
            f"--conda-prefix {smk_conda_env} "
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
                for ffrom, fto in zip(tpmf_out, _out_prefix):
                    Path(fto).parent.mkdir(parents=True, exist_ok=True)
                    shutil.move(ffrom, fto)
                return _out_prefix

    raise NotImplementedError("")
