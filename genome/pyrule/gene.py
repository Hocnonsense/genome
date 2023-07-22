# -*- coding: utf-8 -*-
"""
 * @Date: 2023-07-22 15:34:34
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-07-22 20:29:47
 * @FilePath: /genome/genome/pyrule/gene.py
 * @Description:
"""

from snakemake import shell
from snakemake.workflow import Workflow
from snakemake.io import touch

from . import cache, envs_dir

workflow: Workflow = cache["workflow"]
config: dict = cache["config"]

mantis_config = cache["config"]["mantis_config"]
mantis_config_check = mantis_config + ".check"


@workflow.rule(name="annotate_gene_mantis_check")
@workflow.input(mantis_config=mantis_config)
@workflow.output(mantis_config_check=touch(mantis_config_check))
@workflow.threads(64)
@workflow.conda(envs_dir / "mantis.yaml")
@workflow.shadow("shallow")
@workflow.shellcmd(
    """
        mantis setup \
            --mantis_config {input.mantis_config} \
            -c {threads}
        """
)
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
        """
        mantis setup \
            --mantis_config {input.mantis_config} \
            -c {threads}
        """,
        bench_record=bench_record,
        bench_iteration=bench_iteration,
    )


@workflow.rule(name="annotate_gene_mantis")
@workflow.input(faa="{any}.faa", mantis_config_check=mantis_config_check)
@workflow.output(annot="{any}-{method}.tsv")
@workflow.params(mantis_config=mantis_config)
@workflow.threads(64)
@workflow.conda(envs_dir / "mantis.yaml")
@workflow.wildcard_constraints(method="mantis")
@workflow.shadow("shallow")
@workflow.shellcmd(
    """
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
)
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
        """
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
        """,
        bench_record=bench_record,
        bench_iteration=bench_iteration,
    )
