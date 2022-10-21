# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-15 21:29:41
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-10-21 10:06:07
 * @FilePath: /genome/genome/gene_clust.py
 * @Description:
"""


import os
import shutil
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Iterable, Literal, Union, NamedTuple

from Bio import SeqIO, SeqRecord
from snakemake import main as smk


PathLike = Union[str, Path]
prodigal_mode = Literal["single", "meta"]


class MmseqOut(NamedTuple):
    all_100: Path = Path("gene-clu_100.tsv")
    all_clu: Path = Path("gene-clu.tsv")
    all_clu_faa: Path = Path("gene-clu_rep.faa")

    @classmethod
    def from_prefix(cls, prefix: PathLike):
        return cls(
            Path(str(prefix) + "-clu_100.tsv"),
            Path(str(prefix) + "-clu.tsv"),
            Path(str(prefix) + "-clu_rep.faa"),
        )


def mmseq_clust(
    files: Iterable[PathLike] = None,
    faas: Iterable[SeqRecord.SeqRecord] = None,
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
        target_smk_file = smk_workflow / "gene_clust.smk"
        smk_params = (
            f"-s {target_smk_file} "
            f"{tpmf_out_str} "
            f"--use-conda "
            f"-c{threads} -rp"
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
                    fto.parent.mkdir(parents=True, exist_ok=True)
                    shutil.move(ffrom, fto)
                return _out_prefix

    raise NotImplementedError("")
