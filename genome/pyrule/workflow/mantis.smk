"""
 * @Date: 2023-12-28 13:54:59
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-01-14 17:30:51
 * @FilePath: /genome/genome/pyrule/workflow/mantis.smk
 * @Description:
"""


rule annotate_gene_mantis_checkdb:
    input:
        mantis_config=config["mantis_config"],
    output:
        mantis_config_check=touch(str(config["mantis_config"]) + ".mantis_setup.done"),
    threads: 64
    shadow:
        "shallow"
    conda:
        "../envs/mantis.yaml"
    shell:
        """
        mantis setup \
            --mantis_config {input.mantis_config} \
            -c {threads}
        """


rule annotate_gene_mantis:
    input:
        faa="{any}.faa",
        mantis_config_check=str(config["mantis_config"]) + ".mantis_setup.done",
    output:
        annot="{any}-{method}.tsv",
    params:
        mantis_config=config["mantis_config"],
    wildcard_constraints:
        method="mantis",
    threads: 64
    shadow:
        "shallow"
    conda:
        "../envs/mantis.yaml"
    shell:
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