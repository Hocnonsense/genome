"""
 * @Date: 2022-10-27 19:16:12
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-05-30 10:27:41
 * @FilePath: /genome/genome/pyrule/workflow/binning/single.smk
 * @Description:
"""


rule metabat2:
    input:
        contig="{any}-bins/input/" f"filter_GE{MIN_BIN_CONTIG_LEN}.fa",
        jgi="{any}-bins/input/jgi.tsv",
    output:
        ctg2mag="{any}-bins/single/metabat2_{maxP}_{minS}.tsv",
    params:
        dout="metabat2_{maxP}_{minS}",
        extension="fa",
        maxP="{maxP}",
        minS="{minS}",
        MIN_BIN_CONTIG_LEN=max(int(MIN_BIN_CONTIG_LEN), 1500),
    threads: 1
    conda:
        "../../envs/binning.yaml"
    shadow:
        "shallow"
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
            --minContig {params.MIN_BIN_CONTIG_LEN} \
            --minS {params.minS} --maxP {params.maxP}

        if [ -f {params.dout}/*.{params.extension} ]
        then
            for i in {params.dout}/*.{params.extension}
            do
                binname=$(echo $(basename $i) | sed "s/\\\\.{params.extension}//g")
                grep ">" $i | perl -pe "s/\\n/\\t$binname\\n/g" | perl -pe "s/>//g"
            done \
            > {output.ctg2mag}
        else
            touch {output.ctg2mag}.fail
            touch {output.ctg2mag}
        fi
        """


rule maxbin2:
    input:
        contig="{any}-bins/input/" f"filter_GE{MIN_BIN_CONTIG_LEN}.fa",
        jgi="{any}-bins/input/jgi.tsv",
    output:
        ctg2mag="{any}-bins/single/maxbin2_{markerset}.tsv",
    params:
        dout="maxbin2_{markerset}",
        extension="fasta",
        markerset="{markerset}",
        MIN_BIN_CONTIG_LEN=MIN_BIN_CONTIG_LEN,
    wildcard_constraints:
        markerset="107|40",
    threads: 1
    conda:
        "../../envs/binning.yaml"
    priority: 1
    shadow:
        "shallow"
    shell:
        """
        rm -f {params.dout}
        mkdir -p {params.dout}

        run_MaxBin.pl \
            -min_contig_length {params.MIN_BIN_CONTIG_LEN} \
            -contig {input.contig} \
            -abund {input.jgi} \
            -out {params.dout}/{params.dout} \
            -markerset {wildcards.markerset} -thread {threads}

        if [ -f {params.dout}/*.{params.extension} ]
        then
            for i in {params.dout}/*.{params.extension}
            do
                binname=$(echo $(basename $i) | sed "s/\\\\.{params.extension}//g")
                grep ">" $i | perl -pe "s/\\n/\\t$binname\\n/g" | perl -pe "s/>//g"
            done \
            > {output.ctg2mag}
        else
            touch {output.ctg2mag}.fail
            touch {output.ctg2mag}
        fi
        """


rule concoct:
    input:
        contig="{any}-bins/input/" f"filter_GE{MIN_BIN_CONTIG_LEN}.fa",
        lsbams="{any}-bins/input/bams.ls",
    output:
        ctg2mag="{any}-bins/single/concoct.tsv",
    params:
        dout="concoct",
    threads: 8
    conda:
        "../../envs/concoct.yaml"
    shadow:
        "shallow"
    shell:
        """
        rm -f {params.dout}
        mkdir -p {params.dout}

        cut_up_fasta.py \
            {input.contig} -c 10000 -o 0 \
            --merge_last -b {params.dout}.bed \
            > {params.dout}.fasta
        concoct_coverage_table.py \
            {params.dout}.bed `cat {input.lsbams}` \
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
        contig="{any}-bins/input/" f"filter_GE{MIN_BIN_CONTIG_LEN}.fa",
        lsbams="{any}-bins/input/bams.ls",
    output:
        ctg2mag="{any}-bins/single/metadecoder.tsv",
    params:
        dout="metadecoder",
        extension="fasta",
    threads: 8  # for maxbin
    conda:
        "../../envs/metadecoder.yaml"
    shadow:
        "shallow"
    shell:
        """
        rm -f {params.dout}
        mkdir -p {params.dout}

        sams=""
        for i in `cat {input.lsbams}`
        do
            sam={params.dout}/`basename $i`.sam
            samtools view -h $i -@ 1 -o $sam
            sams="$sams $sam"
        done

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
            -o {params.dout}/{params.dout}

        if [ -f {params.dout}/*.{params.extension} ]
        then
            for i in {params.dout}/*.{params.extension}
            do
                binname=$(echo $(basename $i) | sed "s/\\\\.{params.extension}//g")
                grep ">" $i | perl -pe "s/\\n/\\t$binname\\n/g" | perl -pe "s/>//g"
            done \
            > {output.ctg2mag}
        else
            touch {output.ctg2mag}.fail
            touch {output.ctg2mag}
        fi
        """


rule vamb:
    input:
        contig="{any}-bins/input/" f"filter_GE{MIN_BIN_CONTIG_LEN}.fa",
        jgi="{any}-bins/input/jgi.tsv",
    output:
        ctg2mag="{any}-bins/single/vamb.tsv",
    params:
        dout="vamb",
        extension="fna",
    threads: 1
    conda:
        "../../envs/vamb.yaml"
    shadow:
        "shallow"
    shell:
        """
        rm -f {params.dout} {output.ctg2mag}.fail
        # mkdir -p {params.dout}

        vamb \
            --outdir {params.dout} \
            --fasta {input.contig} \
            --jgi {input.jgi} \
            --minfasta 200000 \
        || touch {output.ctg2mag}.fail

        if [ ! -f {output.ctg2mag}.fail ]
        then
            rename \
                {params.dout}/bins/ {params.dout}/{params.dout}- \
                {params.dout}/bins/*.{params.extension}

            for i in {params.dout}/*.{params.extension}
            do
                binname=$(echo $(basename $i) | sed "s/\\\\.{params.extension}//g")
                grep ">" $i | perl -pe "s/\\n/\\t$binname\\n/g" | perl -pe "s/>//g"
            done \
            > {output.ctg2mag}
        else
            touch {output.ctg2mag}
        fi
        """
