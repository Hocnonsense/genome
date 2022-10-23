

rule create_conda_env_gene_clust:
    output: "{any}reate_conda_env-{gene_clust}-finished"
    wildcard_constraints:
        any = ".*/c|c",
    conda: "../envs/{gene_clust}.yaml"
    shell: "touch {output}"
