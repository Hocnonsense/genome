# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-23 17:07:05
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-01-14 17:17:26
 * @FilePath: /genome/genome/create_conda_env.py
 * @Description:
"""

from pathlib import Path
from tempfile import TemporaryDirectory

from .pyrule import smk, smk_workflow, smk_conda_env


def create_conda_env(*envs: str):
    with TemporaryDirectory() as tmpd:
        rule_files = [f"{tmpd}/create_conda_env-{env}-finished" for env in envs]
        target_smk_file = smk_workflow / "create_conda_env.smk"
        smk_params = (
            f"-s {target_smk_file} "
            f"{' '.join(rule_files)} "
            f"--use-conda "
            f"--conda-create-envs-only "
            f"--conda-prefix {smk_conda_env} "
            f"-c1 -p "
        )
        try:
            print("params:", "snakemake", smk_params)
            smk(smk_params)
        except SystemExit as se:
            if se.code:
                print(se.code, se.with_traceback(None))
                raise RuntimeError("snakemake seems not run successfully.")
            else:
                return True
    raise NotImplementedError("")


def list_envs():
    avail_envs_dir = Path(__file__).parent.parent / "envs"
    return [i.name[:-5] for i in avail_envs_dir.glob("*.yaml")]


def create_conda_env_gene_clust():
    return create_conda_env("gene_clust")


def create_conda_env_prokka():
    return create_conda_env("prokka")


def create_conda_env_binning():
    return create_conda_env("binning", "metadecoder", "vamb")


def create_conda_env_all():
    return create_conda_env(*list_envs())


if __name__ == "__main__":
    print("Creating all envs...")
    create_conda_env_all()
