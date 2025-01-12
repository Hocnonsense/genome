"""
 * @Date: 2022-10-10 15:30:31
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2025-01-12 14:47:05
 * @FilePath: /genome/genome/pyrule/workflow/gene_clust.smk
 * @Description:
    use mmseq to cluster genes

"""

# from genome.pyrule import mmseq_clust_95
#
# mmseq_clust_95.register(workflow, name="gene_clust_workflow")(
#    rules=["_py"], exclude_rules=[], name_modifier="mmseq_clust_95"
# )


from genome.gene_clust import MmseqOut, UniRefClu
from genome.gene_clust import extract as gene_clust_extract

mmseq_clust_95_input_prefix = "{any}"
mo = MmseqOut.from_prefix(mmseq_clust_95_input_prefix)


rule mmseq_clust_95:
    input:
        protein=mmseq_clust_95_input_prefix + ".faa",
    output:
        **mo._asdict(),
    params:
        protein=mmseq_clust_95_input_prefix,
    conda:
        "../envs/gene_clust.yaml"
    threads: 64
    shadow:
        "shallow"
    shell:
        """
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


urc = UniRefClu(mmseq_clust_95_input_prefix)


rule mmseq_uniref_cluster:
    input:
        protein=UniRefClu.in_faa(mmseq_clust_95_input_prefix),
    output:
        uniref100=urc.u100,
        uniref90=urc.u90,
        uniref50=urc.u50,
    conda:
        "../envs/gene_clust.yaml"
    threads: 64
    shadow:
        "shallow"
    shell:
        """
        # https://www.uniprot.org/help/uniref
        rm -f smk-mmseq smk-mmseq-2
        mkdir smk-mmseq smk-mmseq-2
        declare DB=smk-mmseq/gene
        declare DB2=smk-mmseq-2/gene

        mmseqs createdb {input.protein} ${{DB}}

        # region 100% nr
        # https://github.com/soedinglab/MMseqs2/issues/601
        mmseqs cluster ${{DB}} ${{DB}}_clu_100 ${{DB2}} \
            --threads {threads} \
            -c 1.0 --cov-mode 1 --min-seq-id 1.0 --exact-kmer-matching 1
        mmseqs createtsv ${{DB}} ${{DB}} ${{DB}}_clu_100 {output.uniref100} --threads {threads}
        mmseqs createsubdb ${{DB}}_clu_100 ${{DB}} ${{DB}}_100
        # endregion 100% nr

        # region 90% nr
        # https://zhuanlan.zhihu.com/p/643405258
        mmseqs cluster ${{DB}}_100 ${{DB}}_clu_90 ${{DB2}} \
            --threads {threads} \
            --cov-mode 0 -c 0.8 --min-seq-id 0.9
        mmseqs createtsv ${{DB}}_100 ${{DB}}_100 ${{DB}}_clu_90 {output.uniref90} --threads {threads}
        mmseqs createsubdb ${{DB}}_clu_90 ${{DB}}_100 ${{DB}}_90
        # endregion 90% nr

        # region 50% nr
        mmseqs cluster ${{DB}}_90 ${{DB}}_clu_50 ${{DB2}} \
            --threads {threads} \
            --cov-mode 0 -c 0.8 --min-seq-id 0.9
        mmseqs createtsv ${{DB}}_90 ${{DB}}_90 ${{DB}}_clu_50 {output.uniref50} --threads {threads}
        # endregion 50% nr
        #mmseqs createsubdb ${{DB}}_clu_50 ${{DB}}_90 ${{DB}}_50
        #mmseqs convert2fasta ${{DB}}_50 {output.uniref50}.faa
        """


rule mmseq_uniref_cluster_extract:
    input:
        protein="{any}.faa",
        uniref_clu="{any}-uniref{thr}.tsv",
    output:
        uniref_clu="{any}-uniref{thr}.faa",
    run:
        import pandas as pd
        from Bio import SeqIO

        df = pd.read_csv(input.uniref_clu, sep="\t", header=None)
        SeqIO.write(
            gene_clust_extract(df[1].values, [input.protein], []),
            output.uniref_clu,
            "fasta-2line",
        )
