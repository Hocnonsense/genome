# -*- coding: utf-8 -*-
"""
 * @Date: 2023-08-06 18:29:50
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-12-22 20:53:38
 * @FilePath: /genome/genome/pyrule/gunc.py
 * @Description:
"""

from . import general_register, _wf, Path


target_smk_file = Path(__file__).parent / "gunc.smk"

register = general_register(
    snakefile=target_smk_file,
    module_name=__name__.replace(".", "_DOT_"),
    default_config=dict(
        gunc_db_path="data/database/gunc_db",
    ),
)


def register_binning(workflow: _wf.Workflow, rules: _wf.Rules, GUNC_DB: str):
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
