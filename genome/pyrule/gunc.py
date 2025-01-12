# -*- coding: utf-8 -*-
"""
 * @Date: 2023-08-06 18:29:50
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-01-14 17:27:47
 * @FilePath: /genome/genome/pyrule/gunc.py
 * @Description:
"""

from . import general_register, _wf, rules_dir


DEFAULT_GUNC_DB_PATH = "data/database/gunc_db"


target_smk_file = rules_dir / "gunc.smk"

register = general_register(
    snakefile=target_smk_file,
    module_name=__name__.replace(".", "_DOT_"),
    default_config=dict(
        gunc_db_path=DEFAULT_GUNC_DB_PATH,
    ),
)


def register_binning(
    workflow: _wf.Workflow, rules: _wf.Rules, GUNC_DB=DEFAULT_GUNC_DB_PATH
):
    _userule = register(
        workflow, name="gunc_workflow", config=dict(gunc_db_path=GUNC_DB)
    )
    _userule(rules=["gunc_download_db"])

    @(lambda ruleinfo: _userule(rules=["gunc_run_ctg2mag"], ruleinfo=ruleinfo))
    @workflow.input(
        bins_faa="{any}-bins/union/{union_method}{marker}-binsfaa.tsv",
        GUNC_DB=rules.gunc_download_db.output.GUNC_DB,
    )
    @workflow.output(mag2gunc="{any}-bins/filter/{union_method}{marker}-gunc.tsv")
    @workflow.params(bins_faa="{any}-bins/union/{union_method}{marker}-binsfaa")
    @workflow.run
    def _():
        pass
