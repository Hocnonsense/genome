"""
 * @Date: 2024-01-11 20:38:07
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-12-20 14:08:21
 * @FilePath: /meer-omics/meer_omics/rules/genome_pan_concat.smk
 * @Description:
"""


def expand_genomes(f_string="{genome}"):
    def _f(wildcards):
        id2mag_prefix_file = checkpoints.ls_genome_files.get(**wildcards).output.ls
        with open(id2mag_prefix_file) as fi:
            genomes = {
                k: v for k, v in (i.strip().split("\t") for i in fi if i.strip())
            }
        return [
            f_string.format(genome_prefix=genome_prefix, genome=genome, **wildcards)
            for genome, genome_prefix in genomes.items()
        ]

    return _f


def extract_genome_prefix(f_basename="{genome}", suffix=".fa"):
    def _f(wildcards):
        id2mag_prefix_file = checkpoints.ls_genome_files.get(**wildcards).output.ls
        with open(id2mag_prefix_file) as fi:
            for line in fi:
                if not line.strip():
                    continue
                genome_name, genome_prefix = line.strip().split()
                if genome_name == f_basename.format(**wildcards):
                    return genome_prefix + suffix

    return _f


rule genome_pan_transporter:
    params:
        expand_genomes=expand_genomes,
        extract_genome_prefix=extract_genome_prefix,


checkpoint ls_genome_files:
    input:
        ls="{any}{bins_seperator}bins.ls",
    output:
        ls="{any}{bins_seperator}bins/pan.id2prefix.tsv",
    wildcard_constraints:
        bins_seperator=r"\.|-|/|_",
    params:
        suffixes=(".fa", ".fna", ".fasta"),
        # expand_genomes=lambda _: expand_genomes,
        # extract_genome_prefix=lambda _: extract_genome_prefix,
    localrule: True
    run:
        import os

        with open(input.ls) as fi:
            genomes = {
                os.path.splitext(os.path.basename(v))[0]: (
                    os.path.splitext(v.split()[0])[0]
                )
                for v in (i.strip() for i in fi)
                if os.path.splitext(v)[1] in params.suffixes
            }
        with open(output.ls, "w") as fo:
            for g, p in genomes.items():
                fo.write(f"{g}\t{p}\n")


rule concat_bin_fa:
    input:
        ls="{any}{bins_seperator}bins.ls",
    output:
        all_fa="{any}{bins_seperator}bins/pan.concat.fa",
    shell:
        """
        for i in `cat {input.ls}`; do cat $i; done \
        > {output.all_fa}
        """


rule collect_all_f_a:
    input:
        db_f_as=expand_genomes("{genome_prefix}{predicter}.{f_a}"),
    output:
        all_f_a="{any}{bins_seperator}bins/pan{predicter}.concat.{f_a}",
    wildcard_constraints:
        f_a="faa|fna",
    threads: 1
    run:
        with open(output.all_f_a, "w") as fo:
            for genome_ls in input.db_f_as:
                with open(genome_ls) as fi:
                    for line in fi:
                        fo.write(line)


rule collect_all_faa_2binsfaa:
    input:
        db_f_as=expand_genomes("{genome_prefix}{predicter}.faa"),
    output:
        all_faa_tsv="{any}{bins_seperator}bins/union/pan{predicter}-binsfaa.tsv",
        all_faa=directory("{any}{bins_seperator}bins/union/pan{predicter}-binsfaa"),
    params:
        predicter="{predicter}",
    threads: 1
    run:
        shell("mkdir {output.all_faa} -p")
        for i in input:
            j = "".join(Path(i).name.rsplit(params.predicter, 1))
            shell("ln -sr {i} {output.all_faa}/{j}")
        shell("realpath -s {output.all_faa}/*.faa > {output.all_faa_tsv}")


rule collect_all_gff:
    input:
        db_gffs=expand_genomes("{genome_prefix}{predicter}.gff"),
    output:
        all_gff="{any}{bins_seperator}bins/pan{predicter}.concat.gff",
    run:
        with open(output.all_gff, "w") as go:
            # Warning: this is used to avoid duplicated name
            i, last_id = 0, ""
            for genome_ls in input.db_gffs:
                with open(genome_ls) as fi:
                    for line in fi:
                        if line.startswith("##FASTA"):
                            break
                        if not line.startswith("#"):
                            line_values = line.split("\t")
                            line_values_8 = line_values[8].split(";", 1)
                            id_values = line_values_8[0][3:].rsplit("_", 1)
                            if id_values[0] != last_id:
                                i += 1
                                last_id = id_values[0]
                            line = "\t".join(
                                [
                                    *line_values[:8],
                                    f"ID={i}_{id_values[1]};{line_values_8[1]}",
                                ]
                            )
                        go.write(line)
            for genome_ls in input.db_gffs:
                with open(genome_ls) as fi:
                    for line in fi:
                        if line.startswith("##FASTA"):
                            break
                    for line in fi:
                        go.write(line)


rule collect_all_bin_statistic:
    input:
        db_gffs=expand_genomes("{genome_prefix}{predicter}.gff"),
    output:
        bin_statistic="{any}{bins_seperator}bins/pan{predicter}.concat.statistic.tsv",
    params:
        ls="{any}{bins_seperator}bins/pan.id2prefix.tsv",
        f_string=lambda _: "{genome_prefix}{predicter}.gff",
    run:
        import pandas as pd
        from genome.bin_statistic import BinStatisticContainer

        with open(params.ls) as fi:
            genomes = {k: v for k, v in (i.strip().split() for i in fi)}

        bs = pd.DataFrame(
            {
                i: BinStatisticContainer.read_gff(
                    params.f_string.format(genome_prefix=genome_prefix, **wildcards)
                )
                .statistic()
                ._asdict()
                for i, genome_prefix in genomes.items()
            },
        ).T

        bs[
            [
                *("gc", "gc_std", "bp_size"),
                *("max_contig_len", "contigs_num", "contig_n50"),
                *("coding_density", "genes_num"),
            ]
        ].to_csv(output.bin_statistic, sep="\t")


rule collect_all_genome2gene:
    input:
        db_gffs=expand_genomes("{genome_prefix}{predicter}.gff"),
    output:
        genome2gene="{any}{bins_seperator}bins/pan{predicter}-ge33.concat.genome2gene.tsv",
    params:
        ls="{any}{bins_seperator}bins/pan.id2prefix.tsv",
        f_string=lambda _: "{genome_prefix}{predicter}.gff",
    run:
        from genome import gff

        with open(params.ls) as fi:
            genomes = {k: v for k, v in (i.strip().split() for i in fi)}


        with open(output.genome2gene, "w") as fo:
            fo.write("Genome\tGene\n")
            for i, genome_prefix in genomes.items():
                gff_file = params.f_string.format(
                    genome_prefix=genome_prefix, **wildcards
                )
                aas = sorted(
                    gff.Parse(gff_file).extract(),
                    key=lambda x: x.id,
                )
                for aa in aas:
                    fo.write(f"{i}\t{aa.id}\n")


# region fetchMG
rule collect_genome_fetchMG_align:
    input:
        markers=expand_genomes("{genome_prefix}{predicter}-fetchMGs/{marker}.faa"),
    output:
        marker="{any}{bins_seperator}bins/pan{predicter}/fetchMGs_{marker}.collect.faa",
    run:
        from Bio import SeqIO
        from pathlib import Path


        def iter_genome_fetchMG(markers):
            for marker in markers:
                try:
                    record = SeqIO.read(marker, "fasta")
                except ValueError:
                    continue
                record.id = Path(marker).parent.parent.name
                yield record


        SeqIO.write(iter_genome_fetchMG(input.markers), output.marker, "fasta")


# endregion fetchMG


# region genome gene annotation
rule collect_mantis_ko:
    input:
        db_mant_s=expand_genomes("{genome_prefix}{predicter}-mantis.tsv"),
    output:
        gene_kos="{any}{bins_seperator}bins/pan{predicter}.concat.all_KO.tsv",
    params:
        genomes=expand_genomes("{genome}"),
    run:
        import pandas as pd
        from genome.gene_clust import MmseqOut
        from genome.gene_annot import Gene2KO

        g2m = {}

        for genome, mantis in zip(params.genomes, input.db_mant_s):
            gene_annots = pd.concat([Gene2KO(Path(mantis)).get_gene_annots()])
            ko_exploded = (
                gene_annots.apply(lambda x: x.split(":") if x.startswith("K") else [])
                .explode()
                .dropna()
                .reset_index()
                .rename({"index": "All"}, axis=1)
            )
            g2m[genome] = ko_exploded

        pd.concat(g2m, names=["Genome"]).reset_index(["Genome"])[
            ["All", "KO", "Genome"]
        ].to_csv(output.gene_kos, sep="\t", index=False)


rule get_drep_gmodule:
    input:
        gene_kos="{any}.all_KO.tsv",
    output:
        genomeko="{any}.genomeko.csv",
        gmodule="{any}.gmodule.csv",
    localrule: True
    run:
        import pandas as pd
        from kegg_manual import load
        from kegg_manual.data import cache

        genomeko = (
            pd.read_csv(input.gene_kos, sep="\t", index_col=0)
            .dropna()
            .pivot_table(
                values="KO", index="KO", columns="Genome", aggfunc=len, fill_value=0
            )
        )
        genomeko.to_csv(output.genomeko)
        gmodule = load.brite_ko00002_gmodule(genomeko, cache.db_kegg_manual_data)
        gmodule.to_csv(output.gmodule)


# endregion genome gene annotation
