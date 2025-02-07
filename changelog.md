<!--
 * @Date: 2023-08-07 15:18:41
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2025-02-07 11:46:16
 * @FilePath: /genome/changelog.md
 * @Description:
-->

changelog for genome
====================

---

## changelog

- 0.2.4+:
  - feat:
    - `rename_filtered_ls_tsv` to rename bin name after filtering
    - `gff.parse` may make things clearer
      - separate `extract`
      - separate `translate` (as well as alias `_translate`)
    - `gene_statistic` that estimate feature of gene,
      - i.e.:
        - [scu](https://doi.org/10.1093/molbev/mss201)
        - [gc_variability](https://www.nature.com/articles/s41564-017-0008-3)
        - [N-ARSC](https://www.nature.com/articles/s41564-017-0008-3)
        - [C-ARSC](https://www.nature.com/articles/s41564-017-0008-3)
      - the values can be calculated for entire genome using `GeneStatisticContainer.statistic`
    - update `bin_statistic.contig2bin` to `Contig2Bin` and `Binput`
      - `bin_statistic.contig2bin` alias to `bin_statistic.Contig2Bin(contig2bin_tsv, contigs)(outdir)`
      - `bin_statistic_ext.format_bin_input` alias to `bin_statistic.Binput.parse`
      - old APIs will be removed in next few versions
    - update `call_gene_id` usage in `gff.Parse.extract`
      - if `call_gene_id` is a function, the second parameter should be `SeqFeature`
    - rename `smk_workflow` to `rules_dir`,
    - add `rules_dir / "pan_concat.smk"` from meer_omics
    - add `rules_dir / "include_genomedb.smk"` to handle all methods other than binning methods for quick usage
    - add `mmseq_family` to cluster genes
    - `gff_2_fa_label` to extract gene with genome label from gff
    - add `InferGeneId` in gff to better handle gene id inference
    - add `mmseq_species` to cluster genes and related functions in `gene_clust`
    - `BinStatisticContainer` encoding filter using `min_aa_len=33`
    - `aai.smk` to calculate AAI
    - muscle in `tree.smk`
    - extract_fna_ko in `genomedb.smk`
    - drep methods
    - update `gff` and `gene_statistic`  to handle transl_execpt case such as `Sec`
- 0.2.4:
  - feat:
    - `UniRefClu` method to cluster genes
      - `gene_clust.UniRefClu` to clust gene in UniRef standard
      - `mmseq_uniref_cluster` and `mmseq_uniref_cluster_extract` in `smk_workflow / "gene_clust.smk"`
      - `mmseq_clust_95` will be removed in next few versions
  - fix:
    - update `bin_statistic.contig2bin`,
      previously it will open a lot of file and may panic if there are more than 1024 bins.
    - update `pyrule/workflow/binning/single.smk`,
      now it will just touch output and output.fail if nothing binned
    - update `prodigal.prodigal_gff_onethread`, now mask by default
    - update `gff.Parse`
      - a clearer `__init__` function, and silence at that time
      - can extract feature across end of (a circular) genome
  - chore:
    - update `metadecoder` from `1.0.17` to `1.0.19`
    - update tests imports
- 0.2.3:
  - feat!:
    - rename API for snakemake output:
      - `-`: indicate the software
      - `_`: indicate the param select in given software
      - `.`: indicate a format, and is preferred than `-`
    - `_version` for report version when installing this repository
    - MmseqOut
      - allow `MmseqOut.from_aout` to auto recognize prefix from output
      - allow `MmseqOut.load_rep2all` to output "Rep100" if needed
    - tests and docs
- 0.2.2:
  - fix:
    - update gff.Parse
    - bug in binning/filter
- 0.2.1:
  - fix:
    - make mantis more robust
    - use gff.Parse with a gff and another genome easier.
- 0.2.0:
  - feat: change dir of workflow rules and envs and make `genome` capable to be install via git
- 0.1.6:
  - feat:
    - checkm2 that can be enabled by `config["checkm2_db_path"]`
    - support snakemake>=8
- 0.1.5:
  - fix:
    - fix metabat2 min contig length
    - fix gff to faa/fa start
    - fix mantis marker
  - feat:
    - add anchor_yaml to handle snakemake files
    - `pyrule.mmseq_clust_95.register`
    - `pyrule.binning.register`
    - change from prodigal to pyrodigal
      - mode: `meta`, `single`, `gvmeta`
      - change `suffix` allowed values in `prodigal_multithread`
    - change bin_filter
      - the monkey fixes will be removed in next few versions (tests required)
    - change contig2bin and format_bin_input
      - contig2bin return `out_dir (str), bin_id (list[str], without suffix)`
      - format_bin_input return
          `out_dir (str), bin_id (list[str], without suffix), suffix`
- 0.1.4:
  - remove `genome.pyrule.gene`.
  - feat
    - update `genome.binning`. Now it no longer accept old api, and related workflow can be referred directly in snakemake module style!
- 0.1.3:
  - feat
    - add `infer_prodigal_gene_id` and `infer_refseq_gene_id` for gff parser
    - new `BinStatisticContainer.read_gff_parser` for generate `BinStatistic`
- 0.1.2:
  - feat
    - change mantis to `genome.pyrule.mantis`, new include method are welcome
    - update gunc in `genome.pyrule`
  - fix
    - fix bug when running checkm and gunc
- 0.1.1:
  - add mantis rules for snakemake in `genome.pyrule`, only for snakemake>=7.31
  - remove `checkm` and `contig2bin` from `genome.binning`
- 0.1.0:
  - [Version match with meta-snakemake-minimal v0.1.0](http://202.120.45.162:12080/Metabolic_Modeling/genome/releases/tag/version-0.1.0)
- 0.0.3:
  - remove dependence of bcbio-gff (no other change)
- 0.0.2:
  - move functions to new packages
    - changing:
      - `checkm` from `genome.binning` to `genome.bin_statistic_ext`
      - `contig2bin` from `genome.binning` to `genome.bin_statistic`
- 0.0.1: init

# [***$\not$`<!-- @Hwrn -->`*~~`\`~~**](README.md)
