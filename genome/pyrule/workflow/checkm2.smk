"""
 * @Date: 2024-01-10 16:05:33
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-01-20 17:41:05
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
        ctg2faa="{any}-bins/union/{union_method}{marker}-binsfaa.tsv",
    output:
        mag2checkm="{any}-bins/filter/{union_method}{marker}-checkm2.tsv",
    params:
        ctg2faa="{any}-bins/union/{union_method}{marker}-binsfaa",
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


rule filter_fa_via_checkm2:
    input:
        contig="{any}-bins/input/" f"filter_lt.{MIN_BIN_CONTIG_LEN}.fa",
        ctg2mag="{any}-bins/union/{method}{marker}.tsv",
        mag2checkm="{any}-bins/filter/{method}{marker}-checkm2.tsv",
    output:
        lsmags="{any}-bins/filter/{method}{marker}-checkm2-bins.ls",
        mags_tsv="{any}-bins/filter/{method}{marker}-checkm2-bins.tsv",
    params:
        mags="{any}-bins/filter/{method}{marker}-checkm2-bins",
    shadow:
        "shallow"
    run:
        shell("/bin/rm -f smk-fliter smk-fliter.tsv")

        from genome.bin_statistic_ext import format_bin_input
        import pandas as pd

        bin_input_dir, binids, suffix = format_bin_input(
            bin_output=f"smk-fliter/discard",
            bin_input=input.ctg2mag,
            support=input.contig,
        )

        checkm2 = pd.read_csv(input.mag2checkm, sep="\t").rename(
            columns={"Name": "Bin Id"}
        )
        checkm2_filter = checkm2[
            (checkm2["Completeness"] >= 50) & (checkm2["Contamination"] <= 10)
        ]
        for bin_fa in checkm2_filter["Bin Id"]:
            shell(f"mv smk-fliter/discard/{bin_fa}.fa smk-fliter/")

        checkm2_filter.to_csv("smk-fliter.tsv", sep="\t", index=False)
        shell(
            """
            mv smk-fliter.tsv {output.mags_tsv}
            mv smk-fliter {params.mags}
            realpath {params.mags}/*.fa > {output.lsmags}
            """
        )
