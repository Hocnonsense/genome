"""
 * @Date: 2022-10-10 16:48:56
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-10-11 15:36:24
 * @FilePath: /genome/workflow/genome.smk
 * @Description:
"""
import os


def get_env_var(env_var: str, default: str) -> str:
    return os.environ.get(env_var, default)


prokka_output_suffixes = get_env_var("prokka_output_suffixes", "gff").split()


rule prokka:
    input:
        genome = "{any}.fa",
    output:
        filetypes = ["{any}-prokka.{kingdom}."f"{suffix}" for suffix in prokka_output_suffixes],
    conda: "../envs/prokka.yaml"
    params:
        filetypes = prokka_output_suffixes,
        output_prefix = "{any}-prokka.",
        kingdom = "{kingdom}",
    wildcard_constraints:
        kingdom = "Archaea|Bacteria|Mitochondria|Viruses",
    shadow: "shallow"
    shell:
        """
        rm -f smk-prokka

        prokka \
            --outdir smk-prokka \
            --kingdom {params.kingdom}
            --prefix genome \
            {input.genome}

        for i in {params.filetypes}
        do
            mv smk-prokka/prokka.$i \
               {params.output_prefix}$i
        done
        """
