"""
 * @Date: 2022-10-08 11:54:54
* @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
* @LastEditTime: 2025-07-09 16:08:58
* @FilePath: /genome/genome/pyrule/workflow/tree.smk
 * @Description:
    draw tree of mags
"""

from pathlib import Path


FETCHMG_40_AA = [
    *("COG0012", "COG0016", "COG0018", "COG0048", "COG0049", "COG0052"),
    *("COG0080", "COG0085", "COG0087", "COG0088", "COG0090", "COG0091"),
    *("COG0092", "COG0093", "COG0094", "COG0096", "COG0097", "COG0098"),
    *("COG0099", "COG0100", "COG0102", "COG0103", "COG0124", "COG0172"),
    *("COG0184", "COG0185", "COG0186", "COG0197", "COG0200", "COG0201"),
    *("COG0202", "COG0215", "COG0256", "COG0495", "COG0522", "COG0525"),
    *("COG0533", "COG0541", "COG0552"),
]


rule genome_fetchMG:
    input:
        faa="{any}.faa",
        fna="{any}.fna",
    output:
        markers=["{any}-fetchMGs/" f"{marker}" ".faa" for marker in FETCHMG_40_AA],
        scores="{any}-fetchMGs/marker_genes_scores.tsv",
    params:
        out_dir="{any}-fetchMGs",
    shadow:
        "shallow"
    conda:
        "../envs/tree.yaml"
    shell:
        """
        rm -f smk-fetchMGs

        fetchMGs \
            -m extraction \
            {input.faa} \
            -o smk-fetchMGs \
            -d {input.fna} \
            -v -i \
            -t 1
        cp smk-fetchMGs/*.faa {params.out_dir}
        cp smk-fetchMGs/*.marker_genes_scores.table {output.scores}
        """


rule filter_genomes_by_n_markers:
    input:
        markers=[
            "{any}/pan{predicter}/fetchMGs_" f"{marker}" ".collect.faa"
            for marker in FETCHMG_40_AA
        ],
    output:
        genomes="{any}/pan{predicter}.fetchMGs_filter{min_markers}.tsv",
    params:
        min_markers="{min_markers}",
    wildcard_constraints:
        min_markers="\\d+",
    run:
        import pandas as pd
        from typing import Iterable
        from Bio import SeqIO


        def fetchMG_2_alignment_dir(alignment_paths: Iterable[str]):
            markers_: dict[int, dict[str, int]] = {}
            for i, file in enumerate(alignment_paths):
                markers_[i] = {}
                for record in SeqIO.parse(file, "fasta"):
                    markers_[i][record.id] = 1
            markers = pd.DataFrame(markers_).fillna(0)
            return markers


        markers = fetchMG_2_alignment_dir(input.markers)


        def passed_genomes(min_markers):
            return set(markers.index[markers.notna().sum(axis=1) >= min_markers])


        genomes = passed_genomes(int(params.min_markers))
        with open(output.genomes, "w") as go:
            go.write("\n".join(genomes))
            go.flush()


rule filter_markers_by_genome:
    input:
        marker="{any}/pan{predicter}/fetchMGs_{marker}.collect.faa",
        genomes="{any}/pan{predicter}.fetchMGs_filter{min_markers}.tsv",
    output:
        marker="{any}/pan{predicter}/fetchMGs_{marker}.filter{min_markers}.faa",
    run:
        from Bio import SeqIO

        with open(input.genomes) as gi:
            genomes = gi.read().strip().split()

        SeqIO.write(
            (i for i in SeqIO.parse(input.marker, "fasta") if i.id in genomes),
            output.marker,
            "fasta",
        )


rule markers_align_trim_collect:
    input:
        markers=[
            "{any}/pan{predicter}/fetchMGs_" f"{marker}" ".filter{min_markers}.afa"
            for marker in FETCHMG_40_AA
        ],
        genomes="{any}/pan{predicter}.fetchMGs_filter{min_markers}.tsv",
    output:
        marker="{any}/pan{predicter}.fetchMGs_filter{min_markers}.afa",
    run:
        from Bio import AlignIO, SeqIO

        with open(input.genomes) as gi:
            genomes = gi.read().strip().split()
        cat_seqs = {genome: "" for genome in genomes}

        for marker in input.markers:
            for i in SeqIO.parse(marker, "fasta"):
                break
            else:
                continue

            tmp = {k: "" for k in cat_seqs.keys()}
            for record in AlignIO.read(marker, "fasta"):
                cat_seqs[record.id] += record.seq
                tmp.pop(record.id)
            for record_id in tmp:
                cat_seqs[record_id] += "-" * len(record.seq)

        with open(output.marker, "w") as fout:
            for genome, seq in cat_seqs.items():
                print(f">{genome}\n{seq}", file=fout)


ruleorder: faa_mafft > markers_align_trim_collect


# region generic tree
rule faa_mafft:
    input:
        marker="{any}.faa",
    output:
        marker="{any}.afa",
    conda:
        "../envs/tree.yaml"
    threads: 8
    shell:
        """
        mafft --thread {threads} \
            --maxiterate 1000 --localpair \
            {input.marker} \
        >   {output.marker}
        """


rule fna_muscle:
    input:
        marker="{any}.fna",
    output:
        marker="{any}-muscle.afa",
    conda:
        "../envs/tree.yaml"
    threads: 8
    shadow:
        "minimal"
    shell:
        """
        mkdir smk-muscle
        cp {input.marker} smk-muscle/input.fna
        muscle \
            -threads {threads} \
            -align smk-muscle/input.fna \
            -output smk-muscle/output.afa
        mv smk-muscle/output.afa {output.marker}
        """


rule afa_trimal:
    input:
        afa="{any}.afa",
    output:
        trimal="{any}.afa.trimal",
    conda:
        "../envs/tree.yaml"
    shell:
        """
        trimal \
            -automated1 \
            -in {input.afa} \
            -out {output.trimal}
        """


rule trimal_fasttree_wag:
    input:
        trimal="{any}.{aln}",
    output:
        fasttree="{any}.{aln}.fasttree",
    wildcard_constraints:
        aln="trimal|aln",
    conda:
        "../envs/tree.yaml"
    shell:
        """
        FastTree -wag {input.trimal} > {output.fasttree}
        """


rule trimal_iqtree_mfp1kbnni:
    input:
        trimal="{any}.{aln}",
    output:
        iqtree="{any}.{aln}.iqtree",
        contree="{any}.{aln}.contree",
    wildcard_constraints:
        aln="trimal|aln",
    threads: 16
    conda:
        "../envs/tree.yaml"
    shell:
        """
        iqtree \
            -s {input.trimal} \
            -m MFP -B 1000 --bnni \
            -redo \
            -T {threads}
        """


# endregion generic tree


# region phylophlan
rule phylophlan_write_config_file:
    input:
        binsfaa="{any}-binsfaa",
    output:
        cfg="{any}-phylophlan_div.cfg",
    conda:
        "../envs/tree.yaml"
    localrule: True
    shell:
        """
        phylophlan_write_config_file \
            -o {output.cfg} \
            \
            -d a \
            --db_aa diamond \
            --map_dna diamond \
            --map_aa diamond \
            --msa mafft \
            --trim trimal \
            --tree1 iqtree
        """


phylophlan_db_folder = Path(
    config.get("phylophlan_databases", "data/database/phylophlan_databases")
)


rule phylophlan_download_db:
    output:
        db=directory(str(phylophlan_db_folder / "phylophlan")),
    params:
        url_tar="http://cmprod1.cibio.unitn.it/databases/PhyloPhlAn/phylophlan.tar",
        url_md5="http://cmprod1.cibio.unitn.it/databases/PhyloPhlAn/phylophlan.md5",
        db_folder=str(phylophlan_db_folder),
    shell:
        """
        mkdir {db_folder.db}
        cd    {db_folder.db}
            wget {params.url_tar}
            wget {params.url_md5}
            tar -xvf `basename {params.url_tar}`
            bunzip2 -k phylophlan/phylophlan.bz2
        """


rule phylophlan_clade_database:
    output:
        db=str(phylophlan_db_folder / "{clade}"),
    params:
        clade="{clade}",
    conda:
        "../envs/tree.yaml"
    shell:
        """
        phylophlan_setup_database \
            -g {params.clade} \
            -o {output.db} \
            --verbose
        """


rule phylophlan:
    input:
        binsfaa="{any}-binsfaa",
        cfg="{any}-phylophlan_{marker}.cfg",
        db=str(phylophlan_db_folder / "phylophlan"),
    output:
        aln="{any}-phylophlan_{marker}_{speed}.aln",
        tre="{any}-phylophlan_{marker}_{speed}.tre",
    params:
        basename=lambda _: Path(_.any).name,
        speed="{speed}",
        outdir="{any}-phylophlan_{marker}_{speed}",
    wildcard_constraints:
        speed="fast|accurate",
    conda:
        "../envs/tree.yaml"
    threads: 64
    shadow:
        "shallow"
    shell:
        """
        rm -f smk-phylophlan smk-phylophlan-input
        mkdir -p smk-phylophlan-input/{params.basename} smk-phylophlan
        cp {input.binsfaa}/* smk-phylophlan-input/{params.basename}

        phylophlan \
            -d `basename {input.db}` \
            --databases_folder `dirname {input.db}` \
            \
            --diversity high --{params.speed} \
            -f {input.cfg} \
            --min_num_markers 1 \
            \
            -i smk-phylophlan-input/{params.basename} \
            -o smk-phylophlan/phylophlan \
            --nproc {threads}

        ln smk-phylophlan/phylophlan/{params.basename}_concatenated.aln \
            {output.aln}
        ln smk-phylophlan/phylophlan{params.basename}.tre.treefile \
            {output.tre}
        mv smk-phylophlan/phylophlan {params.outdir}
        """


# endregion phylophlan
