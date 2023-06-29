
rule DASTool_create:
    input:
        contig=contig,
        union_methods_ls="/".join([bin_union_dir, "bin_union{marker}-methods.csv"]),
    output:
        ctg2mag="/".join([bin_union_dir, "{union_method}{marker}.tsv"]),
        out_dir=directory("/".join([bin_union_dir, "{union_method}{marker}-dir"])),
    params:
        bin_dir="/".join([bin_union_dir, "{union_method}{marker}-dir"]),
    log:
        "/".join([bin_union_dir, "{union_method}{marker}-dir.log"]),
    wildcard_constraints:
        union_method="dastool",
        marker="-.+|",
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
        contig=contig,
        union_methods_ls="/".join([bin_union_dir, "bin_union{marker}-methods.csv"]),
    output:
        profile=temp(
            directory(
                "/".join([bin_union_dir, "{union_method}{marker}.profile.bin_temp"])
            )
        ),
    wildcard_constraints:
        marker="-.+|",
    run:
        from genome.bin_statistic_ext import contig2bin

        with open(input.union_methods_ls) as f1:
            for i, line in enumerate(f1):
                out_dir = Path(output.profile) / f"{i}"
                contig2bin(
                    outdir=out_dir, contig2bin_tsv=line.strip(), contigs=input.contig
                )


rule UniteM_profile:
    input:
        profile="/".join([bin_union_dir, "{union_method}{marker}.profile.bin_temp"]),
    output:
        profile=directory("/".join([bin_union_dir, "{union_method}{marker}.profile"])),
    params:
        bin_dir="/".join([bin_union_dir, "{union_method}{marker}", "bins"]),
    log:
        log="/".join([bin_union_dir, "{union_method}{marker}.profile.log"]),
    wildcard_constraints:
        union_method="unitem",
        selection="greedy|consensus|unanimous",
        marker="-.+|",
    threads: 64
    conda:
        "../../envs/binunion.yaml"
    shadow:
        "shallow"
    shell:
        """
        rm -f smk-unitem
        mkdir -p smk-unitem

        mv {input.profile} smk-unitem/bins

        find smk-unitem/bins/ \
            -size -5k \
            -type f \
            -exec rm {{}} \;

        unitem profile \
            -b smk-unitem/bins/* \
            -c {threads} \
            smk-unitem
        mv smk-unitem \
            {output.profile}
        """


rule UniteM_refine:
    input:
        profile="/".join([bin_union_dir, "{union_method}{marker}.profile"]),
    output:
        ctg2mag="/".join([bin_union_dir, "{union_method}_{selection}{marker}.tsv"]),
        out_dir=directory(
            "/".join([bin_union_dir, "{union_method}_{selection}{marker}-dir"])
        ),
    params:
        bin_dir="/".join([bin_union_dir, "{union_method}_{selection}{marker}", "bins"]),
        selection="{selection}",
    log:
        log="/".join([bin_union_dir, "{union_method}_{selection}{marker}-dir.log"]),
    wildcard_constraints:
        union_method="unitem",
        selection="greedy|consensus|unanimous",
        marker="-.+|",
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
