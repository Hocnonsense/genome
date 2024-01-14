# -*- coding: utf-8 -*-
"""
 * @Date: 2023-12-21 21:05:25
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-01-14 17:21:49
 * @FilePath: /genome/genome/pyrule/gene_predict.py
 * @Description:
"""


from . import general_register, smk_workflow


target_smk_file = smk_workflow / "genome.smk"

register = general_register(
    snakefile=target_smk_file,
    module_name=__name__.replace(".", "_DOT_"),
    default_config=dict(
        prokka_output_suffixes="gff",
    ),
)
