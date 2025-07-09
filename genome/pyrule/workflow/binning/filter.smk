"""
 * @Date: 2023-12-21 21:28:10
* @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
* @LastEditTime: 2025-07-09 17:32:09
* @FilePath: /genome/genome/pyrule/workflow/binning/filter.smk
 * @Description:
"""

try:
    MIN_BIN_CONTIG_LEN
except NameError:
    MIN_BIN_CONTIG_LEN = int(config.get("MIN_BIN_CONTIG_LEN", 1500))


rule ctg2faa:
    input:
        contig="{any}-bins/input/" f"filter_GE{MIN_BIN_CONTIG_LEN}.fa",
        ctg2mag="{any}-bins/union/{union_method}{marker}.tsv",
    output:
        ctg2faa="{any}-bins/union/{union_method}{marker}-binsfaa.tsv",
    params:
        ctg2faa="{any}-bins/union/{union_method}{marker}-binsfaa",
    wildcard_constraints:
        marker="-[^-]+|",
    shadow:
        "shallow"
    threads: 16
    run:
        shell("rm -rf smk-gene {params.ctg2faa} .snakemake")

        from genome.bin_statistic_ext import format_bin_input
        from genome.prodigal import prodigal_multithread

        bin_input_dir, binids, suffix = format_bin_input(
            bin_output=f"smk-gene/raw",
            bin_input=input.ctg2mag,
            support=input.contig,
        )
        for bin_faa in prodigal_multithread(
            (bin_input_dir / (binid + suffix) for binid in binids),
            mode="single",
            out_dir=bin_input_dir,
            suffix="-ge33.faa",
            threads=threads,
        ):
            bin_faa.rename(str(bin_faa)[:-25] + ".faa")
        shell(
            """
            mkdir smk-gene/genes
            mv smk-gene/raw/*.faa smk-gene/genes
            mv smk-gene/genes {params.ctg2faa}
            realpath {params.ctg2faa}/*.faa > {output.ctg2faa}
            """
        )


rule ctg2faa_checkm:
    input:
        ctg2faa="{any}-bins/union/{union_method}{marker}-binsfaa.tsv",
    output:
        mag2checkm="{any}-bins/filter/{union_method}{marker}-checkm.tsv",
    params:
        ctg2faa="{any}-bins/union/{union_method}{marker}-binsfaa",
    shadow:
        "shallow"
    threads: 64
    run:
        shell("rm -rf smk-checkm")
        from genome.bin_statistic_ext import CheckMFakeOptions

        CheckMFakeOptions(
            file=output.mag2checkm,
            bin_input=str(params.ctg2faa),
            output_dir="{any}-bins/filter/{union_method}{marker}-checkm",
            extension="faa",
            bCalledGenes=True,
            threads=threads,
        ).run()


if config.get("checkm2_db_path"):

    include: "../checkm2.smk"


if config.get("gunc_db_path"):
    from genome.pyrule import gunc

    gunc.register_binning(workflow, rules, config.get("gunc_db_path"))


rule filter_fa_via_bin_filter:
    input:
        contig="{any}-bins/input/" f"filter_GE{MIN_BIN_CONTIG_LEN}.fa",
        ctg2mag="{any}-bins/union/{method}{marker}.tsv",
        mag2checkm="{any}-bins/filter/{method}{marker}-checkm.tsv",
        gunc_out_tsv="{any}-bins/filter/{method}{marker}-gunc.tsv",
    output:
        lsmags="{any}-bins/filter/{method}{marker}-checkmgunc_bins.ls",
        mags_tsv="{any}-bins/filter/{method}{marker}-checkmgunc_bins.tsv",
    params:
        mags="{any}-bins/filter/{method}{marker}-bins",
    threads: 64
    wildcard_constraints:
        marker="-[^-]+|",
        method="dastool|unitem_greedy|unitem_consensus|unitem_unanimous",
    shadow:
        "shallow"
    run:
        shell("/bin/rm -f smk-fliter smk-fliter.tsv")

        from genome.bin_statistic_ext import bin_filter

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


rule filter_fa_via_checkm2:
    input:
        contig="{any}-bins/input/" f"filter_GE{MIN_BIN_CONTIG_LEN}.fa",
        ctg2mag="{any}-bins/union/{method}{marker}.tsv",
        mag2checkm="{any}-bins/union/{method}{marker}-checkm2.tsv",
    output:
        lsmags="{any}-bins/filter/{method}{marker}-checkm2_bins.ls",
        mags_tsv="{any}-bins/filter/{method}{marker}-checkm2_bins.tsv",
    params:
        mags="{any}-bins/filter/{method}{marker}-checkm2_bins",
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


rule rename_filtered_ls_tsv:
    input:
        lsmags="{any}-bins/filter/{method}{marker}-{check}_bins.ls",
        mags_tsv="{any}-bins/filter/{method}{marker}-{check}_bins.tsv",
    output:
        lsmags="{any}-bins/filter/{method}{marker}-{check}_bins-rename/{prefix}_bins.ls",
        mags_tsv="{any}-bins/filter/{method}{marker}-{check}_bins-rename/{prefix}_bins.tsv",
    params:
        mags="{any}-bins/filter/{method}{marker}-{check}_bins-rename/{prefix}_bins",
        new_bin_name=lambda _: (lambda binid: f"{_.prefix}-{binid:>04}",),
    run:
        import os
        import pandas as pd


        def readonly_and_ln_to(file_source: Path, file_target: Path):
            file_source.chmod(0o444)  # make the file readonly
            file_target.parent.mkdir(parents=True, exist_ok=True)
            if file_target.is_file() and file_source.lstat() == file_target.lstat():
                pass
            else:
                file_target.hardlink_to(file_source)
            return file_target


        bin_ls = pd.read_csv(input.lsmags, names=["Path"])
        bin_ls.index = bin_ls["Path"].apply(lambda s: Path(s).name.rsplit(".fa", 1)[0])
        checkm = pd.read_csv(input.mags_tsv, sep="\t")
        common_prefix_len = len(os.path.commonprefix(list(checkm["Bin Id"])))
        checkm_sort = checkm.assign(
            genome=lambda df: df["Bin Id"].apply(
                lambda x: params.new_bin_name[0](x[common_prefix_len:])
            )
        )
        bin_ls_new = bin_ls.merge(
            checkm_sort[["Bin Id", "genome"]],
            left_index=True,
            right_on="Bin Id",
        ).assign(
            New_path=lambda df: df["genome"].apply(
                lambda x: Path(params.mags) / x / (x + ".fa")
            )
        )
        bin_ls_new.apply(
            lambda s: readonly_and_ln_to(Path(s["Path"]), s["New_path"]), axis=1
        )
        bin_ls_new[["New_path"]].to_csv(output.lsmags, index=False, header=False)
        checkm_sort.set_index("genome").reset_index().drop(columns=["Bin Id"]).rename(
            columns={"genome": "Bin Id"}
        ).to_csv(output.mags_tsv, sep="\t", index=False)
