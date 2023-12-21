# -*- coding: utf-8 -*-
"""
 * @Date: 2023-10-31 11:34:55
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-10-31 11:40:48
 * @FilePath: /genome/test/genome/anchor_yaml.py
 * @Description:
"""
# """


from pathlib import Path

import yaml


PathLike = str | Path


class AnchorYaml:
    def __init__(self, yaml_file: PathLike) -> None:
        self._yaml_file = Path(yaml_file)
        self._work_dir = self._yaml_file.parent
        with open(self._yaml_file, "r") as yi:
            self.content: dict | list = yaml.safe_load(yi)

    def get(self, name):
        return self.content[name]
