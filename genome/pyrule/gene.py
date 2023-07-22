# -*- coding: utf-8 -*-
"""
 * @Date: 2023-07-22 15:34:34
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-07-22 16:46:34
 * @FilePath: /genome/pyrule/gene.py
 * @Description:
"""

from snakemake import shell
from snakemake.workflow import Workflow

from . import cache, envs_dir

workflow: Workflow = cache["workflow"]
config: dict = cache["config"]


@workflow.rule(name="annotate_gene_mantis")
@workflow.input(faa="{any}.faa")
@workflow.output(annot="{any}-{method}.tsv")
@workflow.params(mantis_config=config.get("mantis_config", ""))
@workflow.threads(64)
@workflow.conda(envs_dir / "mantis.yaml")
@workflow.wildcard_constraints(method="mantis")
@workflow.shadow("shallow")
@workflow.message(
    """
        # if you didn't create the database, please run:
        mantis setup \
            --mantis_config {params.mantis_config} \
            -c {threads}
        """
)
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
