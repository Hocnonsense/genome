# -*- coding: utf-8 -*-
"""
 * @Date: 2023-12-21 21:05:25
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-06-17 23:29:42
 * @FilePath: /genome/genome/pyrule/gene_predict.py
 * @Description:
"""

from . import general_register, rules_dir

target_smk_file = rules_dir / "genome.smk"

register = general_register(
    snakefile=target_smk_file,
    module_name=__name__.replace(".", "_DOT_"),
    default_config=dict(
        prokka_output_suffixes="gff",
    ),
)
