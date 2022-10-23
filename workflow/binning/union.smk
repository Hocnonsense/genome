
rule DASTool_create:
    input:
        contig  = contig,
        ctg2mag = ["/".join([bin_single, f"{method}.tsv"]) for method in bin_methods],
    output:
        scf2bin = "/".join([bin_single, f"dastool.tsv"]),
        summary = "/".join([bin_single, f"dastool/summary.tsv"]),
        bin_dir = directory("/".join([bin_single, f"dastool/bins"])),
    params:
        ctg2mag = ",".join([str("/".join([bin_single, f"{method}.tsv"])) for method in bin_methods]),
        methods = ",".join(bin_methods),
    log:
        log     = "/".join([bin_single, f"dastool/log"]),
    threads: 8
    conda: "../../envs/binning.yaml"
    #shadow: "shallow"
    shell:
        """
        DAS_Tool \
            -i {params.ctg2mag} -l {params.methods} \
            -c {input.contig} \
            -o DASTool \
            --write_bins 1 --search_engine diamond --score_threshold 0 \
            -t {threads} \
            --debug

        cp DASTool_DASTool.log {log.log}
        cp DASTool_DASTool_scaffolds2bin.txt {output.scf2bin}
        cp DASTool_DASTool_summary.txt {output.summary}
        mv DASTool_DASTool_bins {output.bin_dir}
        """
