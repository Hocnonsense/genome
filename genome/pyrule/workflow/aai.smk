"""
 * @Date: 2025-01-08 17:39:43
 * @Authors: Antonio Camargo antoniop.camargo@gmail.com
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2025-01-17 23:51:23
 * @FilePath: /genome/genome/pyrule/workflow/aai.smk
 * @Description:

 retrieved from:
    https://github.com/apcamargo/bioinformatics-snakemake-pipelines/tree/main/genome-aai-pipeline/genome-aai-pipeline.smk
"""


rule default_aai_ln:
    input:
        faa="{ID}.faa",
        p2g="{ID}-prot2genome.tsv",
        aai="{ID}-diamond_e3i0c0_rbh-ani_s0f30.tsv",
    output:
        aai="{ID}-aai_default.tsv",
    shell:
        "ln {input.aai} {output.aai}"


rule diamond:
    input:
        faa="{ID}.faa",
    output:
        dmd="{ID}-diamond_e{max_eval}i{min_seq_id}.tsv",
    params:
        max_eval="1e-{max_eval}",
        min_seq_id=lambda _: 100 * float(f"0.{_.min_seq_id}"),
    wildcard_constraints:
        max_eval=r"\d+",
        min_seq_id=r"\d+",
    threads: config.get("diamond_threads", 8)
    conda:
        "../envs/tree.yaml"
    shell:
        """
        diamond blastp --threads {threads} \
            -e {params.max_eval} --sensitive \
            --id {params.min_seq_id} \
            -q {input} -d {input} -o {output} \
            --max-target-seqs 25000 \
            --outfmt 6 qseqid sseqid pident qcovhsp scovhsp length
        """


rule filter_diamond:
    input:
        dmd="{ID}-diamond_e{max_eval}i{min_seq_id}.tsv",
        tsv="{ID}-prot2genome.tsv",
    output:
        dmd="{ID}-diamond_e{max_eval}i{min_seq_id}c{min_cov}.tsv",
        rbh="{ID}-diamond_e{max_eval}i{min_seq_id}c{min_cov}_rbh.tsv",
    params:
        min_cov=lambda _: 100 * float(f"0.{_.min_cov}"),
    wildcard_constraints:
        max_eval=r"\d+",
        min_seq_id=r"\d+",
        min_cov=r"\d+",
    run:
        current_qseqid = None
        prot2genome = {}
        with open(input.tsv) as fin:
            for line in fin:
                protein, genome = line.strip().split()
                prot2genome[protein] = genome
        best_hits = {}
        with open(input.dmd) as fin, open(output.dmd, "w") as fout:
            for line in fin:
                qseqid, sseqid, _, qcovhsp, scovhsp, _ = line.strip().split()[:6]
                qseqgen, sseqgen = prot2genome[qseqid], prot2genome[sseqid]
                if qseqgen == sseqgen:
                    continue
                if (params["min_cov"] <= float(qcovhsp)) and (
                    params["min_cov"] <= float(scovhsp)
                ):
                    if qseqid != current_qseqid:
                        current_qseqid, current_qseqid_matches = qseqid, set()
                    if sseqgen not in current_qseqid_matches:
                        current_qseqid_matches.add(sseqgen)
                        fout.write(line)
                        best_hits.setdefault(qseqid, {})[sseqgen] = sseqid
        with open(output.dmd) as fin, open(output.rbh, "w") as fout:
            for line in fin:
                qseqid, sseqid, _, qcovhsp, scovhsp, _ = line.strip().split()[:6]
                qseqgen, sseqgen = prot2genome[qseqid], prot2genome[sseqid]
                if best_hits.get(sseqid, {}).get(qseqgen) == qseqid:
                    fout.write(line)


rule aaicalc:
    input:
        p2g="{ID}-prot2genome.tsv",
        rbh="{ID}-diamond_{label}_rbh.tsv",
    output:
        aai="{ID}-diamond_{label}_rbh-ani_s{min_n_shared_genes}f{min_frac_shared_genes}.tsv",
    params:
        min_cov=lambda _: float(f"0.{_.min_frac_shared_genes}"),
        min_n_shared_genes="{min_n_shared_genes}",
    run:
        gene_count: dict[str, int] = {}
        prot2genome: dict[str, str] = {}
        with open(input.p2g) as fin:
            for line in fin:
                gene, genome_id = line.strip().split()[:2]
                prot2genome[gene] = genome_id
                if genome_id not in gene_count:
                    gene_count[genome_id] = 0
                gene_count[genome_id] += 1
        genome_pairs_genes: dict[frozenset[str], set[frozenset[str]]] = dict()
        gene_pair_info: dict[frozenset[str], tuple[float, int, int]] = dict()
        with open(input.rbh) as fin:
            for line in fin:
                qseqid, sseqid, pident_, _, _, length_ = line.split()
                gene_pair = frozenset({qseqid, sseqid})
                genome_pairs_genes.setdefault(
                    frozenset((prot2genome[i] for i in gene_pair)), set()
                ).add(gene_pair)
                pident, length = float(pident_), int(length_)
                gene_pair_ = gene_pair_info.get(gene_pair, (0, 0, 0))
                gene_pair_info[gene_pair] = (
                    gene_pair_[0] + pident,
                    gene_pair_[1] + length,
                    gene_pair_[2] + 1,
                )
        with open(output.aai, "w") as fout:
            print(
                *("genome_1", "genome_2", "n_genes_genome_1", "n_genes_genome_2"),
                *("n_shared_genes", "aai"),
                sep="\t",
                file=fout,
            )
            for genome_pair, gene_pairs in genome_pairs_genes.items():
                genome_1, genome_2 = sorted(genome_pair)
                min_n_genes = min(gene_count[genome_1], gene_count[genome_2])
                n_shared_genes = len(gene_pairs)
                if (params["min_n_shared_genes"] <= n_shared_genes) and (
                    params["min_cov"] <= n_shared_genes / min_n_genes
                ):
                    pair_aai = 0.0
                    pair_total_length = 0.0
                    for gene_pair in gene_pairs:
                        pidents, lengths, counts = gene_pair_info[gene_pair]
                        pair_aai += (pidents / counts) * (lengths / counts)
                        pair_total_length += lengths / counts
                    print(
                        *(genome_1, genome_2),
                        *(gene_count[genome_1], gene_count[genome_2]),
                        *(n_shared_genes, pair_aai / pair_total_length / 100),
                        sep="\t",
                        file=fout,
                    )
