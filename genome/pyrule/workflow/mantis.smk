"""
 * @Date: 2023-12-28 13:54:59
* @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
* @LastEditTime: 2025-05-21 22:38:27
* @FilePath: /genome/genome/pyrule/workflow/mantis.smk
 * @Description:
"""


rule annotate_gene_mantis_check_db:
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
        mantis_config_check=rules.annotate_gene_mantis_check_db.output.mantis_config_check,
    output:
        annot="{any}-{method}.tsv",
        shadow_folder=directory("{any}-{method}"),
    params:
        mantis_config=config["mantis_config"],
    wildcard_constraints:
        method="mantis",
    threads: 64
    shadow:
        "minimal"
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
            smk-mantis.tsv
        (
            head -n 1 smk-mantis.tsv
            tail -n +2 smk-mantis.tsv | sort -k1,1 -V
        ) > smk-mantis.tsv.sorted
        mv smk-mantis.tsv.sorted \
            {output.annot}
        rm -rf {input.faa}-mantis/
        gzip -r smk-mantis/

        if [[ -e "{output.shadow_folder}" ]]; then
            rm -rf "{output.shadow_folder}"
        fi
        mv smk-mantis {output.shadow_folder}/
        """
