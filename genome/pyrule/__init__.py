# -*- coding: utf-8 -*-
"""
 * @Date: 2023-07-22 15:34:50
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-12-20 20:18:36
 * @FilePath: /genome/genome/pyrule/__init__.py
 * @Description:
"""

from pathlib import Path
from typing import Any, Callable

import snakemake.workflow as _wf


envs_dir = Path(__file__).parent.parent.parent / "envs"


def general_register(
    snakefile: str | Path,
    module_name: str,
    config: dict[str, Any] | None = None,
    ruleinfo: str | Callable[[], None] = lambda: None,
):
    def register(
        workflow: _wf.Workflow,
        name=None,
    ):
        name = name or module_name
        workflow.module(name, snakefile=snakefile, config=config)

        def userule(
            rules=("*",),
            exclude_rules=(),
            name_modifier: str | None = None,
        ):
            return workflow.userule(
                rules=rules,
                from_module=name,
                exclude_rules=exclude_rules,
                name_modifier=name_modifier or None,
            )(ruleinfo)

        return userule

    return register


class __Workflow:
    def __getattribute__(self, __name: str):
        def decorate(*nargs, **kwargs):
            def decorate1(f):
                return f

            return decorate1

        return decorate


cache: dict = {"workflow": __Workflow()}


def register(**kwargs):
    cache.update(**kwargs)
