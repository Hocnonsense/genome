<!--
 * @Date: 2023-08-07 15:18:41
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-09-04 22:37:43
 * @FilePath: /genome/changelog.md
 * @Description:
-->
changelog for genome
===

---
## changelog
- 0.1.3:
    - next changes
        - remove `genome.pyrule.gene` in next minor release.
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



# [***$\not$<!-- @Hwrn -->*~~`\`~~**](README.md)