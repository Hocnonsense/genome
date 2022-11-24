"""
 * @Date: 2022-10-10 16:48:56
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-11-24 23:33:39
 * @FilePath: /genome/workflow/genome.smk
 * @Description:
"""
import os


def get_env_var(env_var: str, default: str) -> str:
    return os.environ.get(env_var, default)


prokka_output_suffixes: list[str] = get_env_var("prokka_output_suffixes", "gff").split()


rule prokka_raw:
    input:
        genome="{any}.fa",
    output:
        filetypes=[
            "{any}-prokka_raw.{kingdom}." f"{suffix}"
            for suffix in prokka_output_suffixes
        ],
        log="{any}-prokka.{kingdom}.log",
    conda:
        "../envs/prokka.yaml"
    params:
        filetypes=prokka_output_suffixes,
        output_prefix="{any}-prokka_raw.{kingdom}.",
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
        filetype="{any}-prokka_raw.{kingdom}.{suffix}",
        log="{any}-prokka.{kingdom}.log",
    output:
        filetype="{any}-prokka.{kingdom}.{suffix}",
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
        faa="{any}-prodigal.{mode}.faa",
        fna="{any}-prodigal.{mode}.fna",
        gff="{any}-prodigal.{mode}.gff",
    conda:
        "../envs/prokka.yaml"
    params:
        mode="{mode}",
    wildcard_constraints:
        mode="single|meta",
    threads: 1
    shell:
        """
        prodigal \
            -i {input.genome} \
            -p {params.mode} \
            -d {output.fna} \
            -a {output.faa} \
            -f gff \
            -o {output.gff} \
            -q

        echo '##FASTA' >> {output.gff}
        cat {input.genome} >> {output.gff}
        """


rule gunc_bins_faa:
    input:
        bins_faa="{any}-bins_faa",
    output:
        gunc_out_tsv="{any}-gunc.tsv",
        gunc_out_dir="{any}-gunc-dir",
    params:
        GUNC_DB=config.get("GUNC_DB", ""),
    conda:
        "../envs/gunc.yaml"
    threads: 32
    shadow:
        "shallow"
    shell:
        """
        rm -f smk-gunc
        mkdir smk-gunc

        gunc run \
            --db_file $GUNC_DB \
            --input_dir {input.bins_faa} \
            --file_suffix .faa \
            --gene_calls \
            --temp_dir smk-gunc \
            --out_dir smk-gunc \
            --threads {threads} \
            --detailed_output

        cp smk-gunc/GUNC.maxCSS_level.tsv {output.gunc_out_tsv}
        mv smk-gunc {output.gunc_out_dir}
        """
