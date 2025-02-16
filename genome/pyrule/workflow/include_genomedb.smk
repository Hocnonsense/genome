"""
 * @Date: 2025-01-12 16:40:11
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2025-02-16 11:24:32
 * @FilePath: /genome/genome/pyrule/workflow/include_genomedb.smk
 * @Description:

 include rules in genome package, intended to be used easily by outer snakemake workflow

 Usage:

'''
import genome.pyrule


include: genome.pyrule.rules_dir / "include_genomedb.smk"


bin_ls = "files/stage2/02_magdb/g__Comamonas/subdb/Comamonas-bins.ls"


rule expected:
    input:
        "files/stage2/02_magdb/g__Comamonas/subdb/Comamonas-bins/pan.id2prefix.tsv",
        "files/stage2/02_magdb/g__Comamonas/subdb/Comamonas-bins/pan-prodigal_single.concat.statistic.tsv",
        "files/stage2/02_magdb/g__Comamonas/subdb/Comamonas-bins/pan-prodigal_single-ge33.concat.genomeko.csv",
'''"""

try:
    import genome
except ImportError:
    import sys as __sys
    from pathlib import Path as __Path

    __current_file = __Path(workflow.snakefile)
    # workflow
    __search_path = __current_file.parents[3]
    __sys.path.append(str(__search_path))
    del __sys, __Path
    import genome


if config.get("checkm2_db_path"):

    include: "checkm2.smk"


include: "gene_clust.smk"
include: "genome.smk"


if config.get("gunc_db_path"):

    include: "gunc.smk"


if config.get("mantis_config"):

    include: "mantis.smk"


include: "pan_concat.smk"
include: "genomedb.smk"
include: "drep.smk"
