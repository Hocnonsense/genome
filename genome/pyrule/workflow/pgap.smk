"""
* @Date: 2025-06-27 21:32:58
* @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
* @LastEditTime: 2025-07-09 16:34:44
* @FilePath: /genome/genome/pyrule/workflow/pgap.smk
* @Description: https://github.com/ncbi/pgap
"""

import re


def encode_id(s: str):
    return "".join(c if re.match(r"[A-Za-z0-9]", c) else f"_{ord(c):02X}" for c in s)


rule pgap_run_auto:
    input:
        fna="{any}.fa",
        PGAP_INPUT_DIR=Path(config.get("PGAP_DB", "data/database/PGAP"))
        / "PGAP_INPUT_DIR",
    output:
        gbk="{any}-pgap.gbk",
        gff="{any}-pgap.gff",
        tax="{any}-pgap.tax.txt",
    params:
        sif=Path(config.get("PGAP_DB", "data/database/PGAP")) / "pgap.sif",
        script_py=Path(config.get("PGAP_DB", "data/database/PGAP")) / "pgap.py",
    shadow:
        "minimal"
    threads: 16
    shell:
        """
        fna=`realpath {input.fna}` && rm smk-pgap -f

        PGAP_INPUT_DIR={input.PGAP_INPUT_DIR} \
        python {params.script_py} \
            --docker apptainer --container-path {params.sif} \
            --report-usage-false --no-internet --no-self-update \
            --cpus {threads} \
            --output smk-pgap \
            --ignore-all-errors \
            --genome {input.fna} \
            --taxcheck --auto-correct-tax --organism Thermoplasma \
            --verbose

        mv smk-pgap/annot.gbk {output.gbk}
        mv smk-pgap/annot_with_genomic_fasta.gff {output.gff}
        mv smk-pgap/ani-tax-report.txt {output.tax}
        """


rule pgap_run_organism:
    input:
        fna="{any}.fa",
        PGAP_INPUT_DIR=Path(config.get("PGAP_DB", "data/database/PGAP"))
        / "PGAP_INPUT_DIR",
    output:
        gbk="{any}-pgap_{organism}.gbk",
        gff="{any}-pgap_{organism}.gff",
    params:
        sif=Path(config.get("PGAP_DB", "data/database/PGAP")) / "pgap.sif",
        script_py=Path(config.get("PGAP_DB", "data/database/PGAP")) / "pgap.py",
        organism="{organism}",
    shadow:
        "minimal"
    threads: 16
    shell:
        """
        fna=`realpath {input.fna}` && rm smk-pgap -f

        PGAP_INPUT_DIR={input.PGAP_INPUT_DIR} \
        python {params.script_py} \
            --docker apptainer --container-path {params.sif} \
            --report-usage-false --no-internet --no-self-update \
            --cpus {threads} \
            --output smk-pgap \
            --ignore-all-errors \
            --genome {input.fna} \
            --organism {params.organism} \
            --verbose

        mv smk-pgap/annot.gbk {output.gbk}
        sed 's/## FASTA/##FASTA/g' smk-pgap/annot_with_genomic_fasta.gff > {output.gff}
        """


rule pgap_download_db:
    output:
        PGAP_INPUT_DIR=Path(config.get("PGAP_DB", "data/database/PGAP"))
        / "PGAP_INPUT_DIR",
    params:
        sif=Path(config.get("PGAP_DB", "data/database/PGAP")) / "pgap.sif",
        script_py=Path(config.get("PGAP_DB", "data/database/PGAP")) / "pgap.py",
    shell:
        """
        PGAP_INPUT_DIR={output.PGAP_INPUT_DIR} \
        python {params.script_py} \
            --docker apptainer --container-path {params.sif}
        """


rule pgap_clean_fa_input:
    input:
        fa="{any}.fa",
    output:
        fa="{any}.ncbi_valid.fa",
        tsv="{any}.ncbi_valid.tsv.gz",
    run:
        from Bio import SeqIO
        import pandas as pd

        name_mappng = {}
        for record in SeqIO.parse(input.fa, "fasta"):
            old_id = record.id
            record.id = encode_id(old_id)
            record.description = ""
            name_mappng[old_id] = record
        SeqIO.write(name_mappng.values(), output.fa, "fasta-2line")
        pd.DataFrame(
            {k: v.id for k, v in zip(("old_id", "new_id"), zip(*name_mappng.items()))}
        ).to_csv(output.tsv, sep="\t", index=False, header=False)
