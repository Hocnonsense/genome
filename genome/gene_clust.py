# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-15 21:29:41
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-08-03 19:28:52
 * @FilePath: /genome/genome/gene_clust.py
 * @Description:
"""


import os
import shutil
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Iterable, Literal, Optional, Union, NamedTuple

from Bio import SeqIO, SeqRecord
from snakemake import main as smk
import pandas as pd


PathLike = Union[str, Path]
prodigal_mode = Literal["single", "meta"]


class MmseqOut(NamedTuple):
    all_100: Path = Path("gene-clu_100.tsv")
    all_clu: Path = Path("gene-clu.tsv")
    all_clu_faa: Path = Path("gene-clu_rep.faa")

    @classmethod
    def from_prefix(cls, prefix: PathLike):
        return cls(
            Path(f"{prefix}-clu_100.tsv"),
            Path(f"{prefix}-clu.tsv"),
            Path(f"{prefix}-clu_rep.faa"),
        )

    @classmethod
    def from_prefix_modify(
        cls,
        prefix: PathLike,
        all_100: Optional[PathLike] = None,
        all_clu: Optional[PathLike] = None,
        all_clu_faa: Optional[PathLike] = None,
    ):
        all_100_ = all_100 if all_100 else f"{prefix}-clu_100.tsv"
        all_clu_ = all_clu if all_clu else f"{prefix}-clu.tsv"
        all_clu_faa_ = all_clu_faa if all_clu_faa else f"{prefix}-clu_rep.faa"

        return cls(Path(all_100_), Path(all_clu_), Path(all_clu_faa_))

    def load_rep2all(self):
        all_100 = pd.read_csv(
            self.all_100, sep="\t", header=None, names=["Rep100", "All"]
        )
        all_clu = pd.read_csv(
            self.all_clu, sep="\t", header=None, names=["Rep", "Rep100"]
        )
        rep2all = all_100.merge(all_clu)[["Rep", "All"]]
        return rep2all


def mmseq_clust(
    files: Optional[Iterable[PathLike]] = None,
    faas: Optional[Iterable[SeqRecord.SeqRecord]] = None,
    out_prefix: Union[PathLike, MmseqOut] = "gene",
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

        smk_workflow = Path(__file__).parent.parent / "workflow"
        smk_conda_env = Path(__file__).parent.parent / ".snakemake" / "conda"
        target_smk_file = smk_workflow / "gene_clust.smk"
        smk_params = (
            f"-s {target_smk_file} "
            f"{tpmf_out_str} "
            f"--use-conda "
            f"--conda-prefix {smk_conda_env} "
            f"-c{threads} -rp "
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
