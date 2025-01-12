# -*- coding: utf-8 -*-
"""
 * @Date: 2023-09-05 11:40:31
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-01-14 17:44:49
 * @FilePath: /genome/genome/pyrule/mmseq_clust_95.py
 * @Description:
"""

from . import general_register, rules_dir


register = general_register(
    snakefile=rules_dir / "gene_clust.smk",
    module_name=__name__.replace(".", "_DOT_"),
    default_config=dict(),
    # "test include 'rules.smk'",
)
