# -*- coding: utf-8 -*-
"""
 * @Date: 2023-07-22 15:34:50
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-07-22 16:13:21
 * @FilePath: /genome/pyrule/__init__.py
 * @Description:
"""

from pathlib import Path


class __Workflow:
    def __getattribute__(self, __name: str):
        def decorate(*nargs, **kwargs):
            def decorate1(f):
                return f

            return decorate1

        return decorate


cache: dict = {"workflow": __Workflow()}
envs_dir = Path(__file__).parent.parent / "envs"


def register(**kwargs):
    cache.update(**kwargs)
