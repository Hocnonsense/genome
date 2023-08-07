# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-13 09:49:04
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-08-07 15:25:45
 * @FilePath: /genome/setup.py
 * @Description:
"""

import os
from pathlib import Path

repo_path = Path(__file__).parent

os.chdir(repo_path)


def get_version(file: str):
    with open(file) as f:
        for line in f:
            if line.startswith("## changelog"):
                break
        line = next(f)
        version = line.strip().rsplit(maxsplit=1)[1].rstrip(":")
    return version


if __name__ == "__main__":
    from setuptools import setup, find_packages

    setup(
        name="genome",
        version=get_version(repo_path / "changelog.md"),
        author="hwrn.aou",
        author_email="hwrn.aou@sjtu.edu.cn",
        description="genome storage and analysis unit",
        # 你要安装的包，通过 setuptools.find_packages 找到当前目录下有哪些包
        packages=find_packages(),
    )
