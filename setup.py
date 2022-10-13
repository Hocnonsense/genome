# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-13 09:49:04
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-10-13 10:17:07
 * @FilePath: /genome/setup.py
 * @Description:
"""

from setuptools import setup, find_packages

setup(
    name="genome",
    version="0.0.1",
    author="hwrn.aou",
    author_email="hwrn.aou@sjtu.edu.cn",
    description="genome storage and analysis unit",
    # 你要安装的包，通过 setuptools.find_packages 找到当前目录下有哪些包
    packages=find_packages(),
)
