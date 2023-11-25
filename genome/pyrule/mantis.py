# -*- coding: utf-8 -*-
"""
 * @Date: 2023-07-22 15:34:34
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-11-05 15:08:45
 * @FilePath: /meta-snakemake-minimal/src/libs/genome/genome/pyrule/mantis.py
 * @Description:
"""

import snakemake.workflow as _wf
from snakemake import shell
from snakemake.io import touch

from . import envs_dir

mantis_setup_shellcmd = """
mantis setup \
    --mantis_config {input.mantis_config} \
    -c {threads}
"""

mantis_run_shellcmd = """
rm -f smk-mantis

mantis \
    run \
    --mantis_config {params.mantis_config} \
    -c {threads} \
    -i {input.faa} \
    -o smk-mantis

mv smk-mantis/consensus_annotation.tsv \
    {output.annot}
mv smk-mantis \
    {input.faa}-mantis/
"""


def register(workflow: _wf.Workflow, mantis_config: str):
    mantis_config_check = str(mantis_config) + ".mantis_setup.done"

    @workflow.rule(name="annotate_gene_mantis_check")
    @workflow.input(mantis_config=mantis_config)
    @workflow.output(mantis_config_check=touch(mantis_config_check))
    @workflow.threads(64)
    @workflow.conda(envs_dir / "mantis.yaml")
    @workflow.shadow("shallow")
    @workflow.shellcmd(mantis_setup_shellcmd)
    @workflow.run
    def __rule_annotate_gene_mantis_check(
        input,
        output,
        params,
        wildcards,
        threads,
        resources,
        log,
        version,
        rule,
        conda_env,
        container_img,
        singularity_args,
        use_singularity,
        env_modules,
        bench_record,
        jobid,
        is_shell,
        bench_iteration,
        cleanup_scripts,
        shadow_dir,
        edit_notebook,
        conda_base_path,
        basedir,
        runtime_sourcecache_path,
        __is_snakemake_rule_func=True,
    ):
        shell(
            mantis_setup_shellcmd,
            bench_record=bench_record,
            bench_iteration=bench_iteration,
        )

    @workflow.rule(name="annotate_gene_mantis")
    @workflow.input(faa="{any}.faa", mantis_config_check=mantis_config_check)
    @workflow.output(annot="{any}-{method}.tsv")
    @workflow.params(mantis_config=mantis_config)
    @workflow.threads(64)
    @workflow.conda(envs_dir / "mantis.yaml")
    @workflow.register_wildcard_constraints(method="mantis")
    @workflow.shadow("shallow")
    @workflow.shellcmd(mantis_run_shellcmd)
    @workflow.run
    def __rule_annotate_gene_mantis(
        input,
        output,
        params,
        wildcards,
        threads,
        resources,
        log,
        version,
        rule,
        conda_env,
        container_img,
        singularity_args,
        use_singularity,
        env_modules,
        bench_record,
        jobid,
        is_shell,
        bench_iteration,
        cleanup_scripts,
        shadow_dir,
        edit_notebook,
        conda_base_path,
        basedir,
        runtime_sourcecache_path,
        __is_snakemake_rule_func=True,
    ):
        shell(
            mantis_run_shellcmd,
            bench_record=bench_record,
            bench_iteration=bench_iteration,
        )
