"""
 * @Date: 2022-10-10 16:48:56
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-10-10 16:54:54
 * @FilePath: /genome/workflow/genome.smk
 * @Description:
"""

rule prokka:
    input:
        genome = "{any}.fa"
    output:
        gff = "{any}-prokka.gff"
    conda: "../envs/genome.yaml"
    shadow: "shallow"
    shell:
        """
        rm -f smk-prokka

        prokka \
            --outdir smk-prokka \
            --prefix genome \
            {input.genome}

        mv smk-prokka/prokka.gff \
           {output.gff}
        """
