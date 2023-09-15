# -*- coding: utf-8 -*-
"""
 * @Date: 2023-09-05 11:40:31
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-09-05 11:55:33
 * @FilePath: /genome/genome/pyrule/mmseq_clust_95.py
 * @Description:
"""

import snakemake.workflow as _wf
from snakemake import shell

from . import envs_dir
from ..gene_clust import MmseqOut

mmseq_clust_95_shellcmd = """
rm -f smk-mmseq smk-mmseq-2
mkdir smk-mmseq smk-mmseq-2
declare DB=smk-mmseq/gene
declare DB2=smk-mmseq-2/gene

mmseqs createdb {input.protein} ${{DB}}

mmseqs cluster ${{DB}} ${{DB}}_clu_100 ${{DB2}} --cov-mode 1 -c 1 --min-seq-id 1 `#-k 10` --threads {threads}
mmseqs createtsv ${{DB}} ${{DB}} ${{DB}}_clu_100 {output.all_100} --threads {threads}
mmseqs createsubdb ${{DB}}_clu_100 ${{DB}} ${{DB}}_100

mmseqs cluster ${{DB}}_100 ${{DB}}_clu ${{DB2}} `#--cov-mode 1` -c 0.9 --min-seq-id 0.95 `#-k 10` --threads {threads}
mmseqs createtsv ${{DB}}_100 ${{DB}}_100 ${{DB}}_clu {output.all_clu} --threads {threads}
mmseqs createsubdb ${{DB}}_clu ${{DB}}_100 ${{DB}}_clu_rep
mmseqs convert2fasta ${{DB}}_clu_rep {output.all_clu_faa}

#mmseqs createseqfiledb ${{DB}} ${{DB}}_clu ${{DB}}_clu_seq --threads {threads}
#mmseqs result2flat ${{DB}} ${{DB}} ${{DB}}_clu_seq ${{clust_out}}-clu_seq.faa.clu
"""


def register(workflow: _wf.Workflow):
    protein = "{any}.faa"
    mo = MmseqOut.from_in_faa(protein)

    @workflow.rule(name="annotate_gene_mantis")
    @workflow.input(protein=protein)
    @workflow.output(all_100=mo.all_100, all_clu=mo.all_clu, all_clu_faa=mo.all_clu_faa)
    @workflow.threads(64)
    @workflow.conda(envs_dir / "gene_clust.yaml")
    @workflow.shadow("shallow")
    @workflow.shellcmd(mmseq_clust_95_shellcmd)
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
            mmseq_clust_95_shellcmd,
            bench_record=bench_record,
            bench_iteration=bench_iteration,
        )
