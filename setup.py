# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-13 09:49:04
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-04-22 22:54:53
 * @FilePath: /genome/setup.py
 * @Description:
"""

import os
from pathlib import Path

repo_path = Path(__file__).parent

os.chdir(repo_path)


def get_version(file: str | Path):
    "first try get version via versionneer, then try to get version from changelog.md"
    try:
        from tests.genome._version import get_versions

        version = get_versions()["version"]
        if version:
            return version
    except ImportError:
        pass
    with open(file) as f:
        for line in f:
            if line.startswith("## changelog"):
                break
        for line in f:
            if line.strip():
                break
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
        packages=find_packages(include=["genome*"]),
        include_package_data=True,
    )
