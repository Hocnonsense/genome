#!/usr/bin/env python
import argparse
import csv
from collections import defaultdict
import os
import errno
from multiprocessing import Pool


def parse_list(listfile, tag_column):
    genome_type = {}
    tag_to_genome = {}
    genomeid_list = []
    with open(listfile) as f:
        next(f)  # skip the header
        for row in csv.reader(f, dialect='excel-tab'):
            #  row[0] is genome id, row[tag_column] is tag id
            #  row[1] is chromosome or plasmid
            genome_type[row[0]] = row[1]
            tag_to_genome[row[tag_column]] = row[0]
            genomeid_list.append(row[0])
    return (genome_type, tag_to_genome, genomeid_list)


def cluster_stat(row, genome_type, tag_to_genome, genomeid_list):
    Rep_genome, Rep_gene = row[0].split('@')  # genomeid@geneid
    Rep_genome = tag_to_genome[Rep_genome]
    Gene_list_chromosome = []
    Gene_list_plasmid = []
    Genome_list = []
    genome_gene_list = defaultdict(list)
    for element in row:
        tag, gene = element.split('@')
        # in case list contains only a subset of genomes
        if tag not in tag_to_genome:
            continue
        genome = tag_to_genome[tag]
        genome_gene_list[genome].append(gene)
        Genome_list.append(genome)
        if genome_type[genome] == 'chromosome':
            Gene_list_chromosome.append(gene)
        else:
            Gene_list_plasmid.append(gene)
    # genome list should be non-redundant, use set() to do that
    Genome_list = set(Genome_list)
    stat = [Rep_gene, '-', Rep_genome, len(Gene_list_chromosome),
            ','.join(Gene_list_chromosome)
            if len(Gene_list_chromosome) > 0
            else '-',
            len(Gene_list_plasmid),
            ','.join(Gene_list_plasmid)
            if len(Gene_list_plasmid) > 0 else '-', len(Genome_list),
            ','.join(Genome_list)]
    gene_count = [len(genome_gene_list[x]) for x in genomeid_list]
    gene_list = [','.join(genome_gene_list[x]) if len(
        genome_gene_list[x]) > 0 else '-' for x in genomeid_list]
    return [stat+gene_count, stat+gene_list]


def parallel_workhorse(tasks):
    func, p = tasks
    return func(*p)


def parse_mcl(listfile, mclfile, tag_column, nproc):
    genome_type, tag_to_genome, genomeid_list = parse_list(
        listfile, tag_column)
    app1 = [["Rep_gene", "Rep_len", "Rep_genome", "Gene_cnt_chromosome",
             "Gene_list_chromosome", "Gene_cnt_plasmid", "Gene_list_plasmid",
             "Genome_cnt", "Genome_list"]+genomeid_list]
    app2 = [["Rep_gene", "Rep_len", "Rep_genome", "Gene_cnt_chromosome",
             "Gene_list_chromosome", "Gene_cnt_plasmid", "Gene_list_plasmid",
             "Genome_cnt", "Genome_list"]+genomeid_list]
    with open(mclfile) as f:
        mcl = csv.reader(f, dialect='excel-tab')
        tasks = [(cluster_stat,
                  (row, genome_type, tag_to_genome, genomeid_list))
                 for row in mcl]
        # Initialize parallel pool of workers
        chunksize = len(tasks) // nproc
        pool = Pool(nproc, maxtasksperchild=chunksize)
        ####
        result = pool.map(parallel_workhorse, tasks, chunksize=chunksize)
        # Finish parallel pool
        pool.close()
        pool.join()
        ####
    app1 += [x[0] for x in result]
    app2 += [x[1] for x in result]
    return(app1, app2)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--list',  help='list file', required=True)
    parser.add_argument(
        '--tag-column',
        help='which column in list file is used as the genome tag (default=1)',
        default=1, type=int)
    parser.add_argument('--mcl',  help='mcl dump file', required=True)
    parser.add_argument('--out',  help='output directory path', required=True)
    parser.add_argument(
        '--nproc', help='number of processes (default=1)', default=1, type=int)

    args = parser.parse_args()

    try:  # make up a temp folder
        os.makedirs(args.out)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    app = parse_mcl(args.list, args.mcl, args.tag_column-1,
                    args.nproc)  # column counts from 0 in python

    listfile = os.path.basename(args.list)
    mclfile = os.path.basename(args.mcl)
    with open(
            os.path.join(
                args.out,
                '_'.join([listfile, mclfile, 'appearance'])),
            'w') as o:
        wr = csv.writer(o, delimiter='\t')
        wr.writerows(app[0])
    with open(
            os.path.join(
                args.out,
                '_'.join([listfile, mclfile, 'appearance2'])),
            'w') as o:
        wr = csv.writer(o, delimiter='\t')
        wr.writerows(app[1])
