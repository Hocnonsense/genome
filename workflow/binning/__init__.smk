MIN_BIN_CONTIG_LEN = config.get("MIN_BIN_CONTIG_LEN", "1500")
bin_methods = config.get("bin_methods", [])


include: "./single.smk"
include: "./union.smk"


rule clean_input_references:
    input:
        config="{any}-bins.yaml",
    output:
        contig="{any}-bins/input/" f"filter_lt.{MIN_BIN_CONTIG_LEN}.fa",
        lsbams="{any}-bins/input/bams.ls",
        jgi="{any}-bins/input/jgi.tsv",
    run:
        import yaml
        from pathlib import Path
        from Bio import SeqIO

        config_path = Path(input.config).parent
        with open(input.config) as yi:
            input_ = {k: config_path / v for k, v in yaml.safe_load(yi).items()}


        SeqIO.write(
            (
                i
                for i in SeqIO.parse(input_["contig"], "fasta")
                if MIN_BIN_CONTIG_LEN <= len(i.seq)
            ),
            output.contig,
            format="fasta",
        )

        shell(
            f"""
            ln -s {input_["lsbams"]} {output.lsbams}
            """
        )

        with (
            open(output.jgi, "w") as fo,
            open(input_["jgi"]) as fi,
            open(output.contig) as fa,
        ):
            fo.write(next(fi))
            jgi: dict[str, str] = {i.split()[0]: i for i in fi}
            for line in fa:
                if line.startswith(">"):
                    fo.write(jgi[line[1:].split()[0]])
            fo.flush()
