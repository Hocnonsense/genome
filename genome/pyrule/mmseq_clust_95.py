# -*- coding: utf-8 -*-
"""
 * @Date: 2023-09-05 11:40:31
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-12-20 21:09:36
 * @FilePath: /genome/genome/pyrule/mmseq_clust_95.py
 * @Description:
"""

from . import envs_dir, general_register, _wf, Path
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


mmseq_clust_95_input_protein = "{any}.faa"
mo = MmseqOut.from_in_faa(mmseq_clust_95_input_protein)


register = general_register(
    snakefile=Path(__file__).parent / "shell_conda_rule.smk",
    module_name=__name__.replace(".", "_DOT_"),
    default_config=dict(
        name="annotate_gene_mantis",
        inputs=dict(protein=mmseq_clust_95_input_protein),
        outputs=mo._asdict(),
        params=dict(protein="{any}"),
        conda=envs_dir / "gene_clust.yaml",
        threads=64,
        shadow="shallow",
        shell=mmseq_clust_95_shellcmd,
    ),
    # "test include 'rules.smk'",
)
