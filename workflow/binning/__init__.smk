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


from pathlib import Path

if not Path(lsbams).is_file():

    rule bams_to_ls:
        input:
            bams=bams,
            bais=[f"{bam}.bai" for bam in bams],
        output:
            lsbams=lsbams,
        shell:
            """
            ls {input.bams} > {output.lsbams}
            """

    rule bam_bai:
        input:
            bam="{any}.bam",
        output:
            bai="{any}.bam.bai",
        conda:
            "../../envs/metadecoder.yaml"
        threads: 8
        shell:
            """
            samtools index {input.bam} -@ {threads}
            """


rule generate_union_methods_ls:
    input:
        ctg2mags=["/".join([bin_single, f"{method}.tsv"]) for method in bin_methods],
    output:
        union_methods_ls="/".join([bin_union_dir, "bin_union{marker}-methods.csv"]),
    wildcard_constraints:
        marker="-.+|",
    shell:
        """
        ls {input.ctg2mags} > {output.union_methods_ls}
        """
