<!--
 * @Date: 2023-08-07 15:18:41
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-05-23 11:06:49
 * @FilePath: /genome/changelog.md
 * @Description:
-->

changelog for genome
====================

---

## changelog

- 0.2.4:
  - fix:
    - update `bin_statistic.contig2bin`,
      previously it will open a lot of file and may panic if there are more than 1024 bins.
- 0.2.3:
  - feat!:
    - rename API for snakemake output:
      - `-`: indicate the software
      - `_`: indicate the param select in given software
      - `.`: indicate a format, and is preferred than `-`
    - `_version` for report version when installing this reposity
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
      - the monkey fixes will be remove in next few versions (tests required)
    - change contig2bin and format_bin_input
      - contig2bin return `out_dir (str), bin_id (list[str], without suffix)`
      - format_bin_input return
          `out_dir (str), bin_id (list[str], without suffix), suffix`
- 0.1.4:
  - remove `genome.pyrule.gene`.
  - feat
    - update `genome.binning`. Now it no longer accept old api, and related workflow can be refered directly in snakemake module style!
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
