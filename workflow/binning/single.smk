
rule filtered_contig:
    input:
        contig  = contig_raw,
    output:
        contig  = temp(contig),
    message:
        "concoct cannot filter short contigs itself"
    run:
        from Bio import SeqIO
        SeqIO.write(
            (i for i in SeqIO.parse(input.contig, "fasta") if len(i.seq) > MIN_BIN_CONTIG_LEN),
            output.contig,
            format="fasta",
        )


rule metabat2:
    input:
        contig  = contig,
        jgi     = jgi,
    output:
        ctg2mag = "/".join([bin_single, "metabat2_{maxP}_{minS}.tsv"]),
    params:
        dout = "metabat2_{maxP}_{minS}",
        extension = "fa",
        maxP = "{maxP}",
        minS = "{minS}",
    threads: 1
    conda: "../../envs/binning.yaml"
    #shadow: "shallow"
    shell:
        """
        rm -f {params.dout}
        # mkdir -p $dout  # DONOT create the directory

        metabat2 \
            -i {input.contig} \
            -a {input.jgi} \
            -o {params.dout}/{params.dout} \
            -t {threads} \
            -v \
            --minContig {MIN_BIN_CONTIG_LEN} \
            --minS {params.minS} --maxP {params.maxP}

        for i in {params.dout}/*.{params.extension}
        do
            binname=$(echo $(basename $i) | sed "s/\\.{params.extension}//g")
            grep ">" $i | perl -pe "s/\n/\t$binname\n/g" | perl -pe "s/>//g"
        done \
        > {output.ctg2mag}
        """


rule maxbin2:
    input:
        contig  = contig,
        jgi     = jgi,
    output:
        ctg2mag = "/".join([bin_single, "maxbin2_{markerset}.tsv"]),
    params:
        dout = "maxbin2_{markerset}",
        extension = "fasta",
        markerset = "{markerset}",
    wildcard_constraints:
        markerset = "107|40",
    threads: 1
    conda: "../../envs/binning.yaml"
    priority: 1
    #shadow: "shallow"
    shell:
        """
        rm -f {params.dout}
        mkdir -p {params.dout}

        run_MaxBin.pl \
            -min_contig_length {MIN_BIN_CONTIG_LEN} \
            -contig {input.contig} \
            -abund {input.jgi} \
            -out {params.dout}/{params.dout} \
            -markerset {wildcards.markerset} -thread {threads}

        for i in {params.dout}/*.{params.extension}
        do
            binname=$(echo $(basename $i) | sed "s/\\.{params.extension}//g")
            grep ">" $i | perl -pe "s/\n/\t$binname\n/g" | perl -pe "s/>//g"
        done \
        > {output.ctg2mag}
        """


rule concoct:
    input:
        contig  = contig,
        bams_ls = bams_ls,
    output:
        ctg2mag = "/".join([bin_single, "concoct.tsv"]),
    params:
        dout = "concoct",
    threads: 8
    conda: "../../envs/binning.yaml"
    #shadow: "shallow"
    shell:
        """
        rm -f {params.dout}
        mkdir -p {params.dout}

        cut_up_fasta.py \
            {input.contig} -c 10000 -o 0 \
            --merge_last -b {params.dout}.bed \
            > {params.dout}.fasta
        concoct_coverage_table.py \
            {params.dout}.bed `cat {input.bams_ls}` \
            > {params.dout}.tsv
        concoct --coverage_file {params.dout}.tsv \
            --composition_file {params.dout}.fasta \
            -t {threads} -r 150 -b {params.dout}/ -s 599 --no_original_data
        merge_cutup_clustering.py \
            {params.dout}/clustering_gt1000.csv \
            > {params.dout}/clustering_gt1000_merge.csv

        awk -v FS="," -v OFS="\t" \
            '{{if (NR==1) {{next}}; {{print $1,"{params.dout}_"$2}}}}' \
            {params.dout}/clustering_gt1000_merge.csv \
        | sort | awk '{{print $2,$1}}' \
        | sort | awk -v OFS="\t" '{{print $2,$1}}' \
        > {output.ctg2mag}
        """


rule metadecoder:
    input:
        contig  = contig,
        bams_ls = bams_ls,
    output:
        ctg2mag = "/".join([bin_single, "metadecoder.tsv"]),
    params:
        dout = "metadecoder",
        extension = "fasta",
    threads: 8  # for maxbin
    conda: "../../envs/metadecoder.yaml"
    #shadow: "shallow"
    shell:
        """
        rm -f {params.dout}
        mkdir -p {params.dout}

        sams=""
        for i in `cat {input.bams_ls}`
        do
            sam={params.dout}/`basename $i`.sam
            samtools view -h $i -@ 1 -o $sam
            sams="$sams $sam"
        metadecoder coverage \
            -s $sams \
            -o {params.dout}/metadecoder.coverage.tsv

        metadecoder seed \
            --threads {threads} \
            -f {input.contig} \
            -o {params.dout}/metadecoder.seed

        metadecoder cluster \
            -f {input.contig} \
            -c {params.dout}/metadecoder.coverage.tsv \
            -s {params.dout}/metadecoder.seed \
            -o {params.dout}

        for i in {params.dout}/*.{params.extension}
        do
            binname=$(echo $(basename $i) | sed "s/\\.{params.extension}//g")
            grep ">" $i | perl -pe "s/\n/\t$binname\n/g" | perl -pe "s/>//g"
        done \
        > {output.ctg2mag}
        """


rule vamb:
    input:
        contig  = contig,
        jgi     = jgi,
    output:
        ctg2mag = "/".join([bin_single, "vamb.tsv"]),
    params:
        dout = "vamb",
        extension = "fna",
    threads: 1
    conda: "../../envs/metadecoder.yaml"
    #shadow: "shallow"
    shell:
        """
        rm -f {params.dout}
        mkdir -p {params.dout}

        vamb \
            --outdir {params.dout} \
            --fasta {input.contig} \
            --jgi {input.jgi} \
            -o C \
            --minfasta 200000

        rename {params.dout}/bins/ {params.dout}/{params.dout}- {params.dout}/bins/*.{params.extension}

        for i in {params.dout}/*.{params.extension}
        do
            binname=$(echo $(basename $i) | sed "s/\\.{params.extension}//g")
            grep ">" $i | perl -pe "s/\n/\t$binname\n/g" | perl -pe "s/>//g"
        done \
        > {output.ctg2mag}
        """