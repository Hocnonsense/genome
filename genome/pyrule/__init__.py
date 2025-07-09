# -*- coding: utf-8 -*-
"""
* @Date: 2023-07-22 15:34:50
* @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
* @LastEditTime: 2025-05-21 22:25:34
* @FilePath: /genome/genome/pyrule/__init__.py
* @Description:
"""

from pathlib import Path
from typing import Any, Callable
import importlib_resources

try:
    from snakemake import main as smk
except ImportError:
    from snakemake.cli import main as smk

import snakemake.workflow as _wf

envs_dir = Path(importlib_resources.files("genome.pyrule")) / "envs"  # type: ignore[reportArgumentType]
smk_conda_env = Path(__file__).parent.parent.parent / ".snakemake" / "conda"
rules_dir = Path(importlib_resources.files("genome.pyrule")) / "workflow"  # type: ignore[reportArgumentType]


def general_register(
    snakefile: str | Path,
    module_name: str,
    default_config: dict[str, Any] | None = None,
):
    """
    Creates a registration function for loading and applying Snakemake modules with optional configuration.
    
    Parameters:
        snakefile (str or Path): Path to the Snakemake file to register.
        module_name (str): Name to assign to the registered module.
        default_config (dict, optional): Default configuration to merge with user-provided config.
    
    Returns:
        Callable: A function that registers the module with a workflow and returns a rule application helper.
    """
    def register(workflow: _wf.Workflow, name=None, config=None):
        name = name or module_name
        workflow.module(
            name,
            snakefile=snakefile,
            config=(default_config or {}) | (config or {}),
        )

        def userule(
            rules=("*",),
            exclude_rules=(),
            name_modifier: str | None = None,
            ruleinfo: str | Callable[[], None] = lambda: None,
        ):
            return workflow.userule(
                rules=rules,
                from_module=name,
                exclude_rules=exclude_rules,
                name_modifier=name_modifier or None,
            )(ruleinfo)

        return userule

    return register


cache: dict[str, _wf.Workflow] = {}


def register(**kwargs):
    cache.update(**kwargs)
