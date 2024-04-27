# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-23 21:35:48
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-05-19 16:47:53
 * @FilePath: /genome/test/genome/test_create_conda_env.py
 * @Description:
__file__ = "test/genome/test_create_conda_env.py"
"""

from genome.create_conda_env import list_envs, create_conda_env_gene_clust

try:
    from _decorator import temp_output, test_temp, test_files, pytest_mark_resource
except (ModuleNotFoundError, ImportError):
    from tests.genome._decorator import (
        temp_output,
        test_temp,
        test_files,
        pytest_mark_resource,
    )

# def test_list_envs():
#    assert sorted(list_envs()) == [
#        "binning",
#        "binunion",
#        "concoct",
#        "gene_clust",
#        "genome",
#        "gunc",
#        "metadecoder",
#        "prokka",
#        "vamb",
#    ]


@pytest_mark_resource
def test_create_conda_env_gene_clust():
    create_conda_env_gene_clust()
