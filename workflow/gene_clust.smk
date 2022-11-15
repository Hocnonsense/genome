"""
 * @Date: 2022-10-10 15:30:31
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-11-13 19:53:48
 * @FilePath: /genome/workflow/gene_clust.smk
 * @Description:
    use mmseq to cluster genes
"""


rule mmseq_clust_95:
    input:
        protein="{any}.faa",
    output:
        all_100="{any}-clu_100.tsv",
        all_clu="{any}-clu.tsv",
        all_clu_faa="{any}-clu_rep.faa",
    params:
        protein="{any}",
    threads: 40
    shadow:
        "shallow"
    conda:
        "../envs/gene_clust.yaml"
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
