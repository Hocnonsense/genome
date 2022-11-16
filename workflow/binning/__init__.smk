MIN_BIN_CONTIG_LEN = config.get("MIN_BIN_CONTIG_LEN", "1500")
contig_raw = config.get("contig", "{any}.fa")
bams = config.get("bams", [])
lsbams = config.get("lsbams", "{any}.bams.ls")
jgi = config.get("jgi", "{any}-jgi.tsv")
bin_single = config.get("bin_single", "{any}-bins/single")
bin_methods = config.get("bin_methods", [])
bin_union_dir = config.get("bin_union_dir", "{any}-bins/union")

contig = "/".join([bin_single, f"contig.{MIN_BIN_CONTIG_LEN}.fa"])


include: "./single.smk"
include: "./union.smk"


rule bams_to_ls:
    input:
        bams=bams,
    output:
        lsbams=lsbams,
    shell:
        """
        ls {input.bams} > {output.lsbams}
        """


rule generate_union_methods_ls:
    input:
        ctg2mags=["/".join([bin_single, f"{method}.tsv"]) for method in bin_methods],
    output:
        union_methods_ls="/".join([bin_union_dir, "bin_union{marker}-methods.csv"]),
    shell:
        """
        ls {input.ctg2mags} > {output.union_methods_ls}
        """
