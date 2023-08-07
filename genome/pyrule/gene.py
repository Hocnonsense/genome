# -*- coding: utf-8 -*-
"""
 * @Date: 2023-07-22 15:34:34
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-08-07 15:17:31
 * @FilePath: /genome/genome/pyrule/gene.py
 * @Description:
"""

from . import cache, mantis

mantis.register(cache["workflow"], cache["config"]["mantis_config"])
