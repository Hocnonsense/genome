"""
 * @Date: 2023-12-21 21:28:10
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-12-22 16:56:24
 * @FilePath: /genome/workflow/binning/filter.smk
 * @Description:
"""

from pathlib import Path


rule checkm_union_tsv:
    input:
        contig="{any}-bins/input/" f"filter_lt.{MIN_BIN_CONTIG_LEN}.fa",
        ctg2mag="{any}-bins/union/{method}-{marker}.tsv",
    output:
        mag2checkm="{any}-bins/filter/{method}-{marker}-checkm.tsv",
    threads: 64
    run:
        from genome.bin_statistic_ext import checkm

        checkm(
            bin_input=input.ctg2mag,
            support=input.contig,
            threads=threads,
        ).to_csv(output.mag2checkm, sep="\t", index=False)


# if config.get("GUNC_DB", ""):
rule gunc_download_db:
    output:
        GUNC_DB=Path(config["gunc_db_path"]) / "gunc_db_progenomes2.1.dmnd",
    params:
        GUNC_DB=config["gunc_db_path"],
    conda:
        "../../envs/gunc.yaml"
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
        contig="{any}-bins/input/" f"filter_lt.{MIN_BIN_CONTIG_LEN}.fa",
        ctg2mag="{any}-bins/union/{method}{marker}.tsv",
    output:
        gunc_out_tsv="{any}-bins/filter/{method}{marker}-gunc.tsv",


rule filter_union_to_fa:
    input:
        contig="{any}-bins/input/" f"filter_lt.{MIN_BIN_CONTIG_LEN}.fa",
        ctg2mag="{any}-bins/union/{method}{marker}.tsv",
        mag2checkm="{any}-bins/filter/{method}-{marker}-checkm.tsv",
        gunc_out_tsv="{any}-bins/filter/{method}{marker}-gunc.tsv",
    output:
        lsmags="{any}-bins/filter/{method}{marker}-bins.ls",
        mags_tsv="{any}-bins/filter/{method}{marker}-bins.tsv",
    params:
        mags="{any}-bins/filter/{method}{marker}-bins",
        GUNC_DB=config["gunc_db_path"],
    threads: 64
    wildcard_constraints:
        marker="-[^-]+|",
        method="dastool|unitem_greedy|unitem_consensus|unitem_unanimous",
    shadow:
        "shallow"
    run:
        shell("/bin/rm -f smk-fliter smk-fliter.tsv")

        mags_tsv = bin_filter(
            bin_out_dir="smk-fliter",
            bin_input=input.ctg2mag,
            support=input.contig,
            checkm_tsv_file=input.mag2checkm,
            gunc_tsv_file=input.gunc_out_tsv,
        )
        mags_tsv.to_csv("smk-fliter.tsv", sep="\t", index=False)

        shell(
            """
            mv smk-fliter.tsv {output.mags_tsv}
            mv smk-fliter {params.mags}
            realpath {params.mags}/*.fa > {output.lsmags}
            """
        )
