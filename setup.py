# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-13 09:49:04
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-04-25 11:13:16
 * @FilePath: /genome/setup.py
 * @Description:
"""

import os
from pathlib import Path

repo_path = Path(__file__).parent

os.chdir(repo_path)


def get_version(file: Path):
    with open(file) as f:
        for line in f:
            if line.startswith("## changelog"):
                break
        for line in f:
            v = line.strip().rsplit(maxsplit=1)
            if len(v) > 1:
                break
        else:
            v = ["", "0+unknown"]
        version = v[1].rstrip(":")
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
