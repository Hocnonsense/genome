"""
 * @Date: 2025-01-13 17:27:32
* @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
* @LastEditTime: 2025-07-09 17:31:28
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
            seqs = gff.to_dict(
                gff.Parse(input.gff).extract(
                    translate="faa" == params.suffix,
                    min_aa_length=int(params.min_aa_len),
                )
            )
            for faa in sorted(seqs.values(), key=lambda x: x.id):
                if faa.seq.startswith("*"):
                    continue
                    # faa.id = f"{faa.id}_partial"
                faa.id = f"{params.chr}{faa.id}"
                SeqIO.write(faa, f, "fasta-2line")


rule fa_label2genome:
    input:
        faa="{any}{predictor}-chr{marker}ge{min_aa_len}.faa",
    output:
        gff="{any}{predictor}-chr{marker}ge{min_aa_len}-prot2genome.tsv",
    params:
        marker="{marker}",
    run:
        with open(output.gff, "w") as f:
            for line in open(input.faa):
                if line.startswith(">"):
                    gene = line[1:].split()[0]
                    genome = gene.split(params.marker)[0]
                    print(gene, genome, sep="\t", file=f)


rule extract_fna_ko:
    input:
        db_gffs=rules.genome_pan_transporter.params.expand_genomes(
            "{genome_prefix}{predicter}.gff"
        ),
        db_mant_s=rules.genome_pan_transporter.params.expand_genomes(
            "{genome_prefix}{predicter}-chr{marker}ge{min_aa_len}-mantis.tsv"
        ),
    output:
        select_f_a="{any}{bins_seperator}bins/union/pan{predicter}-chr{marker}ge{min_aa_len}-{db}={identifier}.{f_a}",
    params:
        genomes=rules.genome_pan_transporter.params.expand_genomes("{genome}"),
        suffix="{f_a}",
        marker="{marker}",
        identifier="{identifier}",
        min_aa_len="{min_aa_len}",
    wildcard_constraints:
        f_a="fna|faa",
        db="ko",
    localrule: True
    run:
        import pandas as pd
        from genome.gene_annot import Gene2KO
        from genome import gff
        import tqdm

        with open(output.select_f_a, "w") as fi:
            for g, gf, mts in zip(params.genomes, input.db_gffs, input.db_mant_s):
                select_gene = frozenset(
                    {
                        i.split(params.marker, 1)[1]
                        for i in Gene2KO(Path(mts))
                        .get_gene_annots()
                        .pipe(lambda s: s[s == params.identifier])
                        .index
                    }
                )
                if len(select_gene) == 0:
                    continue
                seqs = gff.to_dict(
                    gff.Parse(gf).extract(
                        translate="faa" == params.suffix,
                        min_aa_length=int(params.min_aa_len),
                    )
                )
                for seq in seqs.values():
                    if seq.id in select_gene:
                        seq.id = f"{g}{params.marker}{seq.id}"
                        print(seq.format("fasta-2line"), file=fi)
