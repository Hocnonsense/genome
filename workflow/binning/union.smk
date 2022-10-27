
rule DASTool_create:
    input:
        contig  = ancient(contig),
        union_methods_ls = "/".join([bin_union_dir, "{union_method}{marker}", "methods.csv"]),
    output:
        scf2bin = "/".join([bin_union_dir, "{union_method}{marker}.tsv"]),
        bineval = "/".join([bin_union_dir, "{union_method}{marker}", "allBins.eval"]),
        summary = "/".join([bin_union_dir, "{union_method}{marker}", "summary.tsv"]),
    params:
        bin_dir = "/".join([bin_union_dir, "{union_method}{marker}", "bins"]),
    log:
        log     = "/".join([bin_union_dir, "{union_method}{marker}", "log"]),
    wildcard_constraints:
        union_method = "dastool",
        marker = "-.+|",
    threads: 40
    conda: "../../envs/binunion.yaml"
    #shadow: "shallow"
    shell:
        """
        declare bins=`cat {input.union_methods_ls}|tr "$IFS" ','`
        declare labels=$(seq -s, 1 1 $(cat MANIFEST.in |wc -l))

        DAS_Tool \
            -i $bins -l $labels \
            -c {input.contig} \
            -o DASTool \
            --write_bin_evals \
            --search_engine diamond --score_threshold 0 \
            -t {threads} \
            --debug

        cp DASTool_DASTool.log {log.log}
        cp DASTool_DASTool_contig2bin.tsv {output.scf2bin}
        cp DASTool_allBins.eval {output.bineval}
        cp DASTool_DASTool_summary.tsv {output.summary}
        """
