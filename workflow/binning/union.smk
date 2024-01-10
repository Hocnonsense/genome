
rule DASTool_create:
    input:
        contig="{any}-bins/input/" f"filter_lt.{MIN_BIN_CONTIG_LEN}.fa",
        union_methods=[
            ("{any}-bins/single/" f"{bin_method}.tsv") for bin_method in bin_methods
        ],
    output:
        ctg2mag="{any}-bins/union/{union_method}{marker}.tsv",
        out_dir=directory("{any}-bins/union/{union_method}{marker}-dir"),
    log:
        log="{any}-bins/union/{union_method}{marker}.log",
    wildcard_constraints:
        union_method="dastool",
        marker="-[^-]+|",
    threads: 64
    conda:
        "../../envs/binunion.yaml"
    shadow:
        "shallow"
    shell:
        """
        rm -f smk-dastool
        mkdir smk-dastool
        declare bins=`cat {input.union_methods_ls}|tr "$IFS" ','`
        declare labels=$(seq -s, 1 1 $(cat {input.union_methods_ls} |wc -l))

        DAS_Tool \
            -i $bins -l $labels \
            -c {input.contig} \
            -o smk-dastool/DASTool \
            --write_bin_evals \
            --search_engine diamond --score_threshold 0 \
            -t {threads} \
            --debug \
        |tee {log}

        cp smk-dastool/DASTool_DASTool_contig2bin.tsv \
            {output.ctg2mag}
        mv smk-dastool \
            {output.out_dir}
        """


rule UniteM_profile_bin_temp:
    input:
        contig="{any}-bins/input/" f"filter_lt.{MIN_BIN_CONTIG_LEN}.fa",
        ctg2bin="{any}-bins/single/{method}.tsv",
    output:
        profile=temp(
            directory(
                "{any}-bins/union/{union_method}{marker}.profile.bin_temp/{method}"
            )
        ),
    params:
        method="{method}",
    wildcard_constraints:
        marker="-[^-]+|",
    run:
        from genome.bin_statistic import contig2bin

        out_dir = Path(output.profile)
        contig2bin(
            outdir=out_dir,
            contig2bin_tsv=input.ctg2bin,
            contigs=input.contig,
        )


rule UniteM_profile:
    input:
        profile=[
            (
                "{any}-bins/union/{union_method}{marker}.profile.bin_temp/"
                f"{bin_method}"
            )
            for bin_method in bin_methods
        ],
    output:
        profile=directory("{any}-bins/union/{union_method}{marker}.profile"),
    params:
        profile="{any}-bins/union/{union_method}{marker}.profile.bin_temp",
    log:
        log="{any}-bins/union/{union_method}{marker}.profile.log",
    wildcard_constraints:
        union_method="unitem",
        marker="-[^-]+|",
    threads: 64
    conda:
        "../../envs/binunion.yaml"
    shadow:
        "shallow"
    shell:
        """
        rm -f smk-unitem
        mkdir -p smk-unitem

        mv {params.profile} smk-unitem/bins

        find smk-unitem/bins/ \
            -size -5k \
            -type f \
            -exec rm {{}} \\;

        unitem profile \
            -b smk-unitem/bins/* \
            -c {threads} \
            smk-unitem
        mv smk-unitem \
            {output.profile}
        """


rule UniteM_refine:
    input:
        profile="{any}-bins/union/{union_method}{marker}.profile",
    output:
        ctg2mag="{any}-bins/union/{union_method}_{selection}{marker}.tsv",
        out_dir=directory("{any}-bins/union/{union_method}_{selection}{marker}-dir"),
    params:
        selection="{selection}",
    log:
        log="{any}-bins/union/{union_method}_{selection}{marker}.log",
    wildcard_constraints:
        union_method="unitem",
        selection="greedy|consensus|unanimous",
        marker="-[^-]+|",
    threads: 1
    conda:
        "../../envs/binunion.yaml"
    shadow:
        "shallow"
    shell:
        """
        rm -f smk-unitem-out

        unitem {params.selection} \
            -b {input.profile}/bins/* \
            -p unitem_{params.selection}- \
            {input.profile} \
            smk-unitem-out \
        |tee {log}

        for i in smk-unitem-out/bins/*.fna.gz
        do
            binname=$(echo $(basename $i) | sed "s/\\\\.fna.gz//g")
            zcat $i \
            |grep ">" \
            |awk '{{print $1}}' \
            |perl -pe "s/\\n/\\t$binname\\n/g" \
            |perl -pe "s/>//g"
        done \
        > {output.ctg2mag}

        mv smk-unitem-out \
            {output.out_dir}
        """
