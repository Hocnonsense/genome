# -*- coding: utf-8 -*-
"""
* @Date: 2023-09-05 11:40:31
* @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
* @LastEditTime: 2025-05-21 22:27:36
* @FilePath: /genome/genome/pyrule/gene_clust.py
* @Description:
"""

from . import general_register, rules_dir

target_smk_file = rules_dir / "gene_clust.smk"

register = general_register(
    snakefile=target_smk_file, module_name=__name__.replace(".", "_DOT_")
)
