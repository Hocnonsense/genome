"""
 * @Date: 2023-12-22 15:23:06
* @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
* @LastEditTime: 2025-07-09 17:31:31
* @FilePath: /genome/genome/pyrule/workflow/gunc.smk
 * @Description:
"""


rule gunc_download_db:
    output:
        GUNC_DB=Path(config["gunc_db_path"]) / "gunc_db_progenomes2.1.dmnd",
    params:
        GUNC_DB=config["gunc_db_path"],
    conda:
        "../envs/gunc.yaml"
    shadow:
        "shallow"
    shell:
        """
        mkdir -p {params.GUNC_DB}

        gunc download_db {params.GUNC_DB} \
        || (
            wget \
                https://swifter.embl.de/~fullam/gunc/gunc_db_progenomes2.1.dmnd.gz \
                -O {output.GUNC_DB}.gz;
            gunzip {output.GUNC_DB}.gz
        )
        """


rule gunc_run_ctg2mag:
    input:
        bins_faa="{any}",
        GUNC_DB=rules.gunc_download_db.output.GUNC_DB,
    output:
        mag2gunc="{any}-gunc.tsv",
    params:
        bins_faa="{any}",
    conda:
        "../envs/gunc.yaml"
    shadow:
        "shallow"
    threads: 64
    shell:
        """
        rm -f smk-gunc
        mkdir smk-gunc

        gunc run \
            --db_file {input.GUNC_DB} \
            --input_dir {params.bins_faa} \
            --file_suffix .faa \
            --gene_calls \
            --temp_dir smk-gunc \
            --out_dir smk-gunc \
            --threads {threads} \
            --detailed_output

        cp `ls smk-gunc/GUNC.*maxCSS_level.tsv|head -n1` {output.mag2gunc}
        """
