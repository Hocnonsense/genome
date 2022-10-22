
MIN_BIN_CONTIG_LEN  = config.get("MIN_BIN_CONTIG_LEN", "1500")
contig_raw          = config.get("contig", "{any}.fa")
bams_ls             = config.get("bams_ls", "{any}-bams.list")
jgi                 = config.get("jgi", "{any}-jgi.depth")
bin_single          = config.get("bin_single", "{any}-bins/single")
bin_methods         = config.get("bin_methods", [])

contig              = "/".join([bin_single, f"contig.{MIN_BIN_CONTIG_LEN}.fa"])


include: "./single.smk"
include: "./union.smk"
