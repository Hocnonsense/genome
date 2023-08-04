# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-13 09:49:04
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-08-04 11:07:42
 * @FilePath: /genome/setup.py
 * @Description:
"""

import os
from pathlib import Path

os.chdir(Path(__file__).parent)


from setuptools import setup, find_packages

setup(
    name="genome",
    version="0.1.1",
    author="hwrn.aou",
    author_email="hwrn.aou@sjtu.edu.cn",
    description="genome storage and analysis unit",
    # 你要安装的包，通过 setuptools.find_packages 找到当前目录下有哪些包
    packages=find_packages(),
)
