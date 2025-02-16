"""
 * @Date: 2024-01-10 16:05:33
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2025-02-16 11:29:24
 * @FilePath: /genome/genome/pyrule/workflow/checkm2.smk
 * @Description:
"""


rule checkm2_download_db:
    output:
        checkm2_db=directory(config["checkm2_db_path"]),
    conda:
        "../envs/checkm2.yaml"
    shadow:
        "shallow"
    shell:
        """
        checkm2 database --download --path {output.checkm2_db}
        """


rule ctg2faa_checkm2:
    input:
        checkm2_db=rules.checkm2_download_db.output.checkm2_db,
        ctg2faa="{any}-binsfaa.tsv",
    output:
        mag2checkm="{any}-checkm2.tsv",
    params:
        ctg2faa="{any}-binsfaa",
    conda:
        "../envs/checkm2.yaml"
    shadow:
        "shallow"
    threads: 64
    shell:
        """
        export CHECKM2DB=`ls {input.checkm2_db}/*/*.dmnd -t|head -n 1|xargs realpath`

        rm -f smk-checkm-input smk-checkm
        mkdir smk-checkm2-input smk-checkm
        cp `cat {input.ctg2faa}` smk-checkm2-input

        checkm2 predict \
            --threads {threads} \
            --input smk-checkm2-input \
            -x .faa \
            --genes \
            --output-directory smk-checkm2 \
            --force

        mv smk-checkm2/quality_report.tsv {output.mag2checkm}
        """
