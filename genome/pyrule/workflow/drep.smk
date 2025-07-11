"""
 * @Date: 2022-10-04 21:15:46
* @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
* @LastEditTime: 2025-07-09 20:56:51
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
        drep_out=subpath(input.binls, ancestor=1),
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
        drep_out=subpath(input.binls, ancestor=1),
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
        drep_out=subpath(input.binls, ancestor=1),
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
            -g {input.binls} \
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
    run:
        import pandas as pd

        Cdb = pd.read_csv(input.Cdb)
        Bdb = pd.read_csv(input.Bdb)
        Bdb.merge(Cdb)[["location"]].to_csv(output.binls, index=False, header=False)


rule drep2gtdbtk:
    input:
        Wdb="{any}-drep/data_tables/Widb.csv",
        Bdb="{any}-drep/data_tables/Bdb.csv",
    output:
        binls="{any}-dereplicated_bins.ls",
    run:
        import pandas as pd

        Wdb = pd.read_csv(input.Wdb)
        Bdb = pd.read_csv(input.Bdb)
        Wdb.merge(Bdb)[["location"]].to_csv(output.binls, index=False, header=False)


def threads_fastani(wc_l: str | int, maxthreds=64, per_thread_factor=20):
    """
    Calculate threads for FastANI based on number of genomes.
    Divides total comparisons by 20 as each thread can efficiently handle ~20 comparisons.

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
    return min(maxthreds, max(compares // per_thread_factor, 1))


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


rule fastani2af:
    input:
        fastani="{any}bins-fastani.tsv",
        id2prefix="{any}bins/pan.id2prefix.tsv",
    output:
        fastani_af="{any}bins-fastani_af.tsv",
    run:
        from genome.clust import read_fastani

        with open(input.id2prefix) as fi:
            prefix2id = dict(reversed(line.strip().split("\t")[:2]) for line in fi)
            fa2id = lambda v: prefix2id.get(os.path.splitext(v)[0], v)

        ani_long = read_fastani(input.fastani).assign(
            reference=lambda df: df["reference"].apply(fa2id),
            query=lambda df: df["query"].apply(fa2id),
        )
        ani_long.to_csv(output.fastani_af, sep="\t", index=False)


rule fastani2cut:
    input:
        fastani_af="{any}bins-fastani_af.tsv",
    output:
        fastani_cut="{any}bins-fastani_af-{cuts}.tsv",
    params:
        cuts="{cuts}",
    run:
        import pandas as pd
        from genome.clust import run_pairwise_ani

        ani_long = pd.read_csv(input.fastani_af, sep="\t")
        af, ani = 0.5, 0.95
        output_df = None
        for cut in params.cuts.split("-"):
            if "a" in cut:
                ani = float("0." + cut.split("a")[1])
                cut = cut.split("a")[0]
            if "c" in cut:
                af = float("0." + cut.split("c")[1])
            name = f"{af}ani{ani}".replace("ani0.", "fani")
            cx = run_pairwise_ani(ani_long, "average", False, af, ani)[0].rename(
                columns={"cluster": name}
            )
            if output_df is None:
                output_df = cx
            else:
                output_df = output_df.merge(cx, on="genome")
        output_df.set_index("genome").sort_index().to_csv(output.fastani_cut, sep="\t")
