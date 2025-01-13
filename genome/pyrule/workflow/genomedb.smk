"""
 * @Date: 2025-01-13 17:27:32
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2025-01-13 20:49:07
 * @FilePath: /genome/genome/pyrule/workflow/genomedb.smk
 * @Description:
"""


rule gff_2_fa_label:
    input:
        gff="{any}/{genome}/{genome}{predictor}.gff",
    output:
        faa="{any}/{genome}/{genome}{predictor}-chr{marker}ge{min_aa_len}.{suffix}",
    params:
        min_aa_len="{min_aa_len}",
        chr="{genome}{marker}",
        suffix="{suffix}",
    wildcard_constraints:
        genome=r"[^/]+",
        predictor=r"[^/]+",
        suffix="faa|fna",
        #marker=r"(+|-|_|@|~|^|:|=|\||?)+",
        min_aa_len=r"\d+",
    threads: 1
    run:
        from Bio import SeqIO
        from genome import gff

        with open(output.faa, "w") as f:
            for faa in gff.Parse(input.gff).extract(
                translate=params.suffix == "faa",
                min_aa_length=int(params.min_aa_len),
            ):
                if faa.seq.starswith("*"):
                    continue
                    # faa.id = f"{faa.id}_partial"
                faa.id = f"{params.chr}{faa.id}"
                SeqIO.write(faa, f, "fasta-2line")
