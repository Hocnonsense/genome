# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-23 21:35:48
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-06-18 11:14:14
 * @FilePath: /genome/tests/genome/test_create_conda_env.py
 * @Description:
__file__ = "test/genome/test_create_conda_env.py"
"""

from genome.create_conda_env import list_envs, create_conda_env_gene_clust

from tests.genome._decorator import pytest_mark_resource

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
