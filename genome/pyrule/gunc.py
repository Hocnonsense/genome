# -*- coding: utf-8 -*-
"""
 * @Date: 2023-08-06 18:29:50
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-08-07 13:41:21
 * @FilePath: /genome/genome/pyrule/gunc.py
 * @Description:
"""

import snakemake.workflow as _wf
from snakemake import shell
from snakemake.io import directory

from . import envs_dir

gunc_download_db = """
mdkir -p `dirname {output.GUNC_DB}`

gunc download_db {output.GUNC_DB}
"""

gunc_run_shellcmd = """
rm -f smk-gunc
mkdir smk-gunc

gunc run \
    --db_file {input.GUNC_DB} \
    --input_dir {input.bins_faa} \
    --file_suffix .faa \
    --gene_calls \
    --temp_dir smk-gunc \
    --out_dir smk-gunc \
    --threads {threads} \
    --detailed_output

cp `ls smk-gunc/GUNC.*maxCSS_level.tsv|head -n1` {output.gunc_out_tsv}
mv smk-gunc {output.gunc_out_dir}
"""


def register(workflow: _wf.Workflow, GUNC_DB: str):
    @workflow.rule(name="gunc_download_db")
    @workflow.output(GUNC_DB=GUNC_DB)
    @workflow.threads(64)
    @workflow.conda(envs_dir / "gunc.yaml")
    @workflow.shadow("shallow")
    @workflow.shellcmd(gunc_download_db)
    @workflow.run
    def __rule_gunc_download_db(
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
            gunc_download_db,
            bench_record=bench_record,
            bench_iteration=bench_iteration,
        )

    @workflow.rule(name="gunc_run")
    @workflow.input(bins_faa="{any}-bins_faa", GUNC_DB=GUNC_DB)
    @workflow.output(
        gunc_out_tsv="{any}-gunc.tsv", gunc_out_dir=directory("{any}-gunc-dir")
    )
    @workflow.threads(64)
    @workflow.conda(envs_dir / "gunc.yaml")
    @workflow.shadow("shallow")
    @workflow.shellcmd(gunc_run_shellcmd)
    @workflow.run
    def __rule_gunc_run(
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
            gunc_run_shellcmd,
            bench_record=bench_record,
            bench_iteration=bench_iteration,
        )
