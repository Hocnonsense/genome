"""
 * @Date: 2022-10-04 21:15:46
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2025-01-23 14:09:21
 * @FilePath: /genome/genome/pyrule/workflow/drep.smk
 * @Description:
"""


rule drep_create:
    input:
        binls="{any}-drep/bins.ls",
        checkms="{any}-drep/bins.csv",
    output:
        drep_out=[
            ("{any}-drep/data_tables/" f"{table}.csv")
            for table in (
                *("Bdb", "Cdb", "genomeInfo", "genomeInformation"),
                *("Mdb", "Ndb", "Sdb", "Wdb", "Widb"),
            )
        ],
    params:
        drep_out="{any}-drep",
    conda:
        "../envs/drep.yaml"
    threads: 64
    shell:
        """
        mkdir -p {params.drep_out}

        #   --primary_chunksize `cat {input.binls} |wc -l` \
        # by default, dRep will split genome into chunks of 5000,
        #   and then merge them with `mash paste` and finally compare them directly
        dRep dereplicate \
            {params.drep_out} \
            -p {threads} \
            -comp 50 -con 10 \
            -pa 0.9 -sa 0.95 \
            -g {input.binls} \
            --genomeInfo {input.checkms} \
            --debug
        """


rule drep_create_uncheck:
    input:
        binls="{any}-drep/bins_uncheck.ls",
    output:
        drep_out=[
            ("{any}-drep/data_tables/" f"{table}.csv")
            for table in (
                *("Bdb", "Cdb", "genomeInfo", "genomeInformation"),
                *("Mdb", "Ndb", "Sdb", "Wdb", "Widb"),
            )
        ],
    params:
        drep_out="{any}-drep",
    conda:
        "../envs/drep.yaml"
    threads: 64
    shell:
        """
        mkdir -p {params.drep_out}

        dRep dereplicate \
            {params.drep_out} \
            -p {threads} \
            -comp 50 -con 10 \
            -pa 0.9 -sa 0.95 \
            -g {input.binls} \
            --debug

        rm -r {params.drep_out}/data
        """


rule drep_create_nocheck:
    input:
        binls="{any}-drep/bins_nocheck.ls",
    output:
        drep_out=[
            ("{any}-drep/data_tables/" f"{table}.csv") for table in ("Bdb", "Cdb")
        ],
    params:
        drep_out="02_magdb/meer_rv2_vs_pub_{trench}/subdb/refgenome-drep",
    conda:
        "../envs/drep.yaml"
    threads: 64
    shell:
        """
        mkdir -p {params.drep_out}

        dRep dereplicate \
            {params.drep_out} \
            -p {threads} \
            --ignoreGenomeQuality \
            -pa 0.9 -sa 0.95 \
            -g {params.drep_out}/bins_nocheck.ls \
            --skip_plots \
            --debug

        rm -r {params.drep_out}/data
        """


ruleorder: drep_create > drep_create_uncheck > drep_create_nocheck


rule drep2redundant:
    input:
        Bdb="{any}-drep/data_tables/Bdb.csv",
        Cdb="{any}-drep/data_tables/Cdb.csv",
    output:
        binls="{any}-redundant_bins.ls",
    params:
        drep_wd="{any}-drep",
    run:
        import pandas as pd
        from pathlib import Path

        wd = Path(params.drep_wd)

        Cdb = pd.read_csv(wd / "data_tables" / "Cdb.csv")
        Bdb = pd.read_csv(wd / "data_tables" / "Bdb.csv")
        Bdb.merge(Cdb)[["location"]].to_csv(output.binls, index=False, header=False)


rule drep2gtdbtk:
    input:
        Wdb="{any}-drep/data_tables/Widb.csv",
        Bdb="{any}-drep/data_tables/Bdb.csv",
    output:
        binls="{any}-dereplicated_bins.ls",
    params:
        drep_wd="{any}-drep",
    run:
        import pandas as pd
        from pathlib import Path

        wd = Path(params.drep_wd)

        Wdb = pd.read_csv(wd / "data_tables" / "Wdb.csv")
        Bdb = pd.read_csv(wd / "data_tables" / "Bdb.csv")
        Wdb.merge(Bdb)[["location"]].to_csv(output.binls, index=False, header=False)


def threads_fastani(wc_l: str | int, maxthreds=64):
    """
    >>> threads_fastani(1)
    1
    >>> threads_fastani(4)
    1
    >>> threads_fastani(8)
    3
    >>> threads_fastani(16)
    12
    >>> threads_fastani(20)
    20
    >>> threads_fastani(64)
    64
    """
    compares = int(wc_l) ** 2
    return min(maxthreds, max(compares // 20, 1))


rule fastani:
    input:
        binls="{any}bins.ls",
    output:
        fastani="{any}bins-fastani.tsv",
    threads: lambda _, input: threads_fastani(shell(f"wc -l {input.binls}", read=True).split()[0])
    conda:
        "../envs/drep.yaml"
    shell:
        """
        fastANI -t {threads} \\
            --ql {input.binls} \\
            --rl {input.binls} \\
            --minFraction 0 \\
            -o {output.fastani}
        """
