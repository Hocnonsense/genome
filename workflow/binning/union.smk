
rule DASTool_create:
    input:
        contig=ancient(contig),
        union_methods_ls="/".join(
            [bin_union_dir, "{union_method}{marker}", "methods.csv"]
        ),
    output:
        ctg2mag="/".join([bin_union_dir, "{union_method}{marker}.tsv"]),
        out_dir=directory("/".join([bin_union_dir, "{union_method}{marker}-dir"])),
    params:
        bin_dir="/".join([bin_union_dir, "{union_method}{marker}", "bins"]),
    log:
        "/".join([bin_union_dir, "{union_method}{marker}", "log"]),
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


rule UniteM_profile:
    input:
        contig=ancient(contig),
        union_methods_ls="/".join([bin_union_dir, "{union_method}{marker}-methods.csv"]),
    output:
        profile=directory("/".join([bin_union_dir, "{union_method}{marker}.profile"])),
    params:
        bin_dir="/".join([bin_union_dir, "{union_method}{marker}", "bins"]),
    log:
        log="/".join([bin_union_dir, "{union_method}{marker}", "log"]),
    wildcard_constraints:
        union_method="unitem",
        selection="greedy|consensus|unanimous",
        marker="-.+|",
    threads: 64
    conda:
        "../../envs/binunion.yaml"
    # shadow: "shallow"
    shell:
        """
        rm -f smk-unitem
        mkdir smk-unitem

        declare labels_i=1
        for i in `cat {input.union_methods_ls}`
        do
            Contigs2Bin_to_Fasta.sh \
                -i $i \
                -a {input.contig} \
                -o smk-unitem/bins/$labels_i

            find smk-unitem/bins/$labels_i \
                -size -5k \
                -type f \
                -exec rm {{}} \\;
            let labels_i+=1
        done

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
        log="/".join([bin_union_dir, "{union_method}_{selection}{marker}", "log"]),
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
            smk-unitem-out

        cp smk-unitem-out/bins/*.gz \

        for i in smk-unitem-out/bins/*.fna.gz
        do
            binname=$(echo $(basename $i) | sed "s/\\\\.fna.gz//g")
            zcat $i |grep ">" | perl -pe "s/\\n/\\t$binname\\n/g" | perl -pe "s/>//g"
        done \
        > {output.ctg2mag}

        mv smk-unitem-out \
            {output.out_dir}
        """
