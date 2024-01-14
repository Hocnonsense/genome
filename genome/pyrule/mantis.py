# -*- coding: utf-8 -*-
"""
 * @Date: 2023-07-22 15:34:34
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-12-28 14:03:24
 * @FilePath: /genome/genome/pyrule/mantis.py
 * @Description:
"""
from . import general_register, smk_workflow

register = general_register(
    snakefile=smk_workflow / "mantis.smk",
    module_name=__name__.replace(".", "_DOT_"),
    default_config=dict(
        threads=64,
        # mantis_config=str(mantis_config),
    ),
    # "test include 'rules.smk'",
)
