MIN_BIN_CONTIG_LEN = config.get("MIN_BIN_CONTIG_LEN", "1500")


module binning_workflow:
    snakefile:
        "./__init__.smk"
    config:
        {
            "MIN_BIN_CONTIG_LEN": MIN_BIN_CONTIG_LEN,
            "bin_methods": config.get("bin_methods", []),
        }


use rule * from binning_workflow as binning_*


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


use rule binning_clean_input_references as clean_input_references_python with:
    output:
        contig="{any}-bins/input/" f"filter_lt.{MIN_BIN_CONTIG_LEN}.fa",
        lsbams="{any}-bins/input/bams_uncheck.ls",
        jgi="{any}-bins/input/jgi.tsv",


checkpoint check_bam_bai:
    input:
        lsbams="{any}-bins/input/bams_uncheck.ls",
    output:
        lsbams="{any}-bins/input/bams_tocheck.ls",
    shell:
        "cp {input} {output}"


def get_bam_to_check_bai(wildcards):
    with open(checkpoints.check_bam_bai.get(**wildcards).output.lsbams) as fi:
        bami = [i.strip() + ".bai" for i in fi.read().split() if i.strip()]
    return bami


rule checked_bam_bai:
    input:
        lsbams="{any}-bins/input/bams_tocheck.ls",
        lsbais=get_bam_to_check_bai,
    output:
        lsbams="{any}-bins/input/bams.ls",
    shell:
        "cp {input.lsbams} {output.lsbams}"


ruleorder: checked_bam_bai > clean_input_references_python > binning_clean_input_references
