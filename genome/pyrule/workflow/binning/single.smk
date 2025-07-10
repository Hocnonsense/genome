"""
 * @Date: 2022-10-27 19:16:12
* @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
* @LastEditTime: 2025-07-10 17:40:57
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
        folder="metabat2_{maxP}_{minS}",
        extension="fa",
        maxP="{maxP}",
        minS="{minS}",
        MIN_BIN_CONTIG_LEN=max(int(MIN_BIN_CONTIG_LEN), 1500),
    threads: 1
    conda:
        "../../envs/binning.yaml"
    shadow:
        "minimal"
    shell:
        """
        rm -f {params.folder}
        # mkdir -p $folder  # DO NOT create the directory

        metabat2 \
            -i {input.contig} \
            -a {input.jgi} \
            -o {params.folder}/{params.folder} \
            -t {threads} \
            -v \
            --minContig {params.MIN_BIN_CONTIG_LEN} \
            --minS {params.minS} --maxP {params.maxP}

        if [ -f {params.folder}/*.{params.extension} ]
        then
            for i in {params.folder}/*.{params.extension}
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
        folder="maxbin2_{markerset}",
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
        "minimal"
    shell:
        """
        rm -f {params.folder}
        mkdir -p {params.folder}

        run_MaxBin.pl \
            -min_contig_length {params.MIN_BIN_CONTIG_LEN} \
            -contig {input.contig} \
            -abund {input.jgi} \
            -out {params.folder}/{params.folder} \
            -markerset {wildcards.markerset} -thread {threads}

        if [ -f {params.folder}/*.{params.extension} ]
        then
            for i in {params.folder}/*.{params.extension}
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
        folder="concoct",
    threads: 8
    conda:
        "../../envs/concoct.yaml"
    shadow:
        "minimal"
    shell:
        """
        rm -f {params.folder}
        mkdir -p {params.folder}

        cut_up_fasta.py \
            {input.contig} -c 10000 -o 0 \
            --merge_last -b {params.folder}.bed \
            > {params.folder}.fasta
        concoct_coverage_table.py \
            {params.folder}.bed `cat {input.lsbams}` \
            > {params.folder}.tsv
        concoct --coverage_file {params.folder}.tsv \
            --composition_file {params.folder}.fasta \
            -t {threads} -r 150 -b {params.folder}/ -s 599 --no_original_data
        merge_cutup_clustering.py \
            {params.folder}/clustering_gt1000.csv \
            > {params.folder}/clustering_gt1000_merge.csv

        awk -v FS="," -v OFS="\t" \
            '{{if (NR==1) {{next}}; {{print $1,"{params.folder}_"$2}}}}' \
            {params.folder}/clustering_gt1000_merge.csv \
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
        folder="metadecoder",
        extension="fasta",
    threads: 8  # for maxbin
    conda:
        "../../envs/metadecoder.yaml"
    shadow:
        "minimal"
    shell:
        """
        rm -f {params.folder}
        mkdir -p {params.folder}

        sams=""
        for i in `cat {input.lsbams}`
        do
            sam={params.folder}/`basename $i`.sam
            samtools view -h $i -@ 1 -o $sam
            sams="$sams $sam"
        done

        metadecoder coverage \
            -s $sams \
            -o {params.folder}/metadecoder.coverage.tsv

        metadecoder seed \
            --threads {threads} \
            -f {input.contig} \
            -o {params.folder}/metadecoder.seed

        metadecoder cluster \
            -f {input.contig} \
            -c {params.folder}/metadecoder.coverage.tsv \
            -s {params.folder}/metadecoder.seed \
            -o {params.folder}/{params.folder}

        if [ -f {params.folder}/*.{params.extension} ]
        then
            for i in {params.folder}/*.{params.extension}
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
        folder="vamb",
        extension="fna",
    threads: 1
    conda:
        "../../envs/vamb.yaml"
    shadow:
        "minimal"
    shell:
        """
        rm -f {params.folder} {output.ctg2mag}.fail
        # mkdir -p {params.folder}

        vamb \
            --outdir {params.folder} \
            --fasta {input.contig} \
            --jgi {input.jgi} \
            --minfasta 200000 \
        || touch {output.ctg2mag}.fail

        if [ ! -f {output.ctg2mag}.fail ]
        then
            rename \
                {params.folder}/bins/ {params.folder}/{params.folder}- \
                {params.folder}/bins/*.{params.extension}

            for i in {params.folder}/*.{params.extension}
            do
                binname=$(echo $(basename $i) | sed "s/\\\\.{params.extension}//g")
                grep ">" $i | perl -pe "s/\\n/\\t$binname\\n/g" | perl -pe "s/>//g"
            done \
            > {output.ctg2mag}
        else
            touch {output.ctg2mag}
        fi
        """


rule rosella:
    input:
        contig="{any}-bins/input/" f"filter_GE{MIN_BIN_CONTIG_LEN}.fa",
        jgi="{any}-bins/input/jgi.tsv",
    output:
        ctg2mag="{any}-bins/single/rosella.tsv",
    params:
        folder="rosella",
        extension="fna",
    threads: 1
    conda:
        "../../envs/rosella.yaml"
    shadow:
        "minimal"
    shell:
        """
        rm -f {params.folder} {output.ctg2mag}.fail
        # mkdir -p {params.folder}

        rosella recover \
            -C {input.jgi} \
            -r {input.contig} \
            --output-directory {params.folder} \
            -t {threads} \
        || touch {output.ctg2mag}.fail

        rm -f {params.folder}/**_unbinned.{params.extension}
        if [ -f {params.folder}/*.{params.extension} ] && [ ! -f {output.ctg2mag}.fail ]
        then
            for i in {params.folder}/*.{params.extension}
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
