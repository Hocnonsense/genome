"""
 * @Date: 2022-10-10 16:48:56
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-03-29 15:25:40
 * @FilePath: /genome/genome/pyrule/workflow/genome.smk
 * @Description:
"""

import os


prokka_output_suffixes = config.get(
    "prokka_output_suffixes", os.environ.get("prokka_output_suffixes", "gff").split()
)


rule gff_2_fa:
    input:
        gff="{any}.gff",
    output:
        faa="{any}-ge{min_aa_len}.{suffix}",
    params:
        min_aa_len="{min_aa_len}",
        suffix="{suffix}",
    wildcard_constraints:
        suffix="faa|fna",
        min_aa_len="\\d+",
    threads: 1
    run:
        from Bio import SeqIO
        from genome import gff

        SeqIO.write(
            sorted(
                gff.Parse(input.gff).extract(
                    translate=params.suffix == "faa",
                    min_aa_length=int(params.min_aa_len),
                ),
                key=lambda x: x.id,
            ),
            output.faa,
            "fasta-2line",
        )


rule prokka_raw:
    input:
        genome="{any}.fa",
    output:
        filetypes=[
            "{any}-prokka_{kingdom}_raw." f"{suffix}"
            for suffix in prokka_output_suffixes
        ],
        log="{any}-prokka_{kingdom}.log",
    conda:
        "../envs/prokka.yaml"
    params:
        filetypes=prokka_output_suffixes,
        output_prefix="{any}-prokka_{kingdom}_raw.",
        kingdom="{kingdom}",
    wildcard_constraints:
        kingdom="Archaea|Bacteria|Mitochondria|Viruses",
    threads: 8
    shadow:
        "shallow"
    shell:
        """
        rm -f smk-prokka

        prokka \
            --outdir smk-prokka \
            --kingdom {params.kingdom} \
            --cpus {threads} \
            --prefix genome \
            {input.genome} \
        2>&1 \
        |tee {output.log}

        for i in {params.filetypes}
        do
            mv smk-prokka/genome.$i \
               {params.output_prefix}$i
        done
        """


rule prokka_fix:
    input:
        filetype="{any}-prokka_{kingdom}_raw.{suffix}",
        log="{any}-prokka_{kingdom}.log",
    output:
        filetype="{any}-prokka_{kingdom}.{suffix}",
    run:
        from pathlib import Path

        rename_dict: dict[str, str] = {}
        i = 0
        with open(input.log) as log_in:
            for line in log_in:
                if "Changing illegal '|' to '_' in sequence name:" in line:
                    real_name = line.strip().rsplit(" ", 1)[1]
                    rename_dict[real_name.replace("|", "_")] = real_name
                    i += 1
        if len(rename_dict) != i:
            exit(1, "prokka cannot used to annot this genome.")
        with open(input.filetype) as fi, open(output.filetype, "w") as fo:
            for line in fi:
                for fake_name in rename_dict:
                    if fake_name in line:
                        line = line.replace(fake_name, rename_dict[fake_name])
                fo.write(line)

        Path(input.filetype).unlink()


rule prodigal_raw:
    input:
        genome="{any}.fa",
    output:
        gff="{any}-prodigal_{mode}.gff",
    params:
        mode="{mode}",
    wildcard_constraints:
        mode="single|meta|gvmeta",
    threads: 1
    run:
        from genome.prodigal import prodigal_gff_onethread

        prodigal_gff_onethread(input.genome, params.mode, output.gff)
