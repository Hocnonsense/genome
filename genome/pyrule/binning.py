# -*- coding: utf-8 -*-
"""
 * @Date: 2023-12-20 20:58:03
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-01-14 17:20:55
 * @FilePath: /genome/genome/pyrule/binning.py
 * @Description:
"""
# """

from . import general_register, rules_dir
from .. import binning


target_smk_file = rules_dir / "binning" / "__init__.smk"

register = general_register(
    snakefile=target_smk_file,
    module_name=__name__.replace(".", "_DOT_"),
    default_config=dict(
        MIN_BIN_CONTIG_LEN=2500,
        bin_methods=binning.default_bin_methods,
    ),
    # "test include 'rules.smk'",
)
