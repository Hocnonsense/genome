#!/usr/bin/env python
"""Prepare diamond or blast alignment file into MCL-required format.

Also known as the genemap process. Only self connections and bi-directional
best hits are retained.
"""
import pandas as pd
import numpy as np
import argparse
import os
import sys
import tempfile
import subprocess
import shutil
import logging
from contextlib import contextmanager


def parse_ptt(ptt_dir, ncid):
    """Parse ptt file."""
    ptt = pd.read_csv(ptt_dir + '/' + ncid + '.ptt',
                      header=None,
                      usecols=[2, 3],
                      names=['len', 'geneid'],
                      sep='\t')
    return ptt


def parse_alignment(blast_dir, ncid1, ncid2, ptt1, ptt2, ident_cutoff,
                    evalue_cutoff, cov_cutoff):
    """Parse alignment result."""
    align = pd.read_csv(blast_dir + '/' + ncid1 + '_' + ncid2 + '.blast',
                        header=None,
                        names=['qseqid', 'sseqid', 'pident', 'length',
                               'mismatch', 'gapopen', 'qstart', 'qend',
                               'sstart', 'send', 'evalue', 'bitscore'],
                        sep='\t')
    align = pd.merge(pd.merge(align,
                              ptt1,
                              left_on='qseqid', right_on='geneid'),
                     ptt2,
                     left_on='sseqid', right_on='geneid')
    # If there's no best hit at all, just skip the whole session
    if align.empty:
        return align
    # calculate the alignment coverage
    align.eval("""
    qcov = (qend-qstart+1)/len_x
    scov = (send-sstart+1)/len_y
    """, inplace=True)
    # filtering base on threshold
    align.query(('pident>@ident_cutoff & evalue<@evalue_cutoff & '
                 'qcov>@cov_cutoff & scov>@cov_cutoff'), inplace=True)
    align.reset_index(drop=True, inplace=True)
    # find the top1 best match for each query gene
    align_best_values = align.sort_values(by='bitscore',
                                          ascending=False
                                          ).groupby('qseqid',
                                                    sort=False).head(1)
    # show the best bitscore for each query gene
    align_merge = pd.merge(align, align_best_values, how='left', on='qseqid')
    # find the rows that have the best values
    # can be miltiple rows for the same query gene
    align_best_index = align_merge.query('bitscore_x == bitscore_y').index
    # get the final best hit list
    align_best = align.loc[align_best_index]
    if ncid1 == ncid2:  # self alignment
        # E-value of self-alignment should be 0
        self_align_index = align_best.query('qseqid == sseqid').index
        align_best.loc[self_align_index, 'evalue'] = 0
        # find the query genes that do not have self alignment
        gene_unaligned_index = -ptt1['geneid'].isin(
            align_best.loc[self_align_index, 'qseqid'])
        gene_unaligned = ptt1['geneid'][gene_unaligned_index]
        # build up fake alignment
        gene_unaligned = pd.DataFrame({'qseqid': gene_unaligned,
                                       'sseqid': gene_unaligned,
                                       'evalue': 0
                                       })
        # concat alignments
        align_best = pd.concat([align_best, gene_unaligned], sort=False)

    return align_best


def output_alignment(align1: pd.DataFrame,
                     align2: pd.DataFrame,
                     tag1: str, tag2: str, out: str):
    """Write down the filtered alignments into temp file."""
    # If no best hit was found, skip the whole session
    if align1.empty or align2.empty:
        return
    # find the bi-directional best hit
    align_merge = pd.merge(align1, align2, left_on=['qseqid', 'sseqid'],
                           right_on=['sseqid', 'qseqid'], how='inner')
    output = pd.DataFrame({
        'qseqid': tag1 + '@' + align_merge['qseqid_x'],
        'sseqid': tag2 + '@' + align_merge['sseqid_x'],
    })
    # use the lower E-value
    evalue = align_merge[['evalue_x', 'evalue_y']].min(axis=1)
    # negative log10 transform for mcl
    # make the maximum E-value to 200
    evalue = np.negative(np.log10(evalue + 1e-200))
    output['evalue'] = evalue
    if tag1 == tag2:
        # add a new list to check duplicated rows
        # with reversed qseqid and sseqid
        output['checkname'] = list(
            map(lambda x, y: ''.join(sorted([x, y])),
                output['qseqid'].values,
                output['sseqid'].values))
        # drop duplicates
        output.drop_duplicates('checkname', inplace=True)
        del output['checkname']

    output.to_csv(out, header=False, index=False,
                  sep='\t', float_format='%.2g')


def filter_alignment(ncid1, ncid2, tag1, tag2, blast_dir, ptt_dir, out,
                     ident_cutoff, evalue_cutoff, cov_cutoff):
    """Parse, filter and output alignment results."""
    ptt1 = parse_ptt(ptt_dir, ncid1)
    ptt2 = parse_ptt(ptt_dir, ncid2)
    align1 = parse_alignment(blast_dir, ncid1, ncid2, ptt1, ptt2, ident_cutoff,
                             evalue_cutoff, cov_cutoff)
    align2 = parse_alignment(blast_dir, ncid2, ncid1, ptt2, ptt1, ident_cutoff,
                             evalue_cutoff, cov_cutoff)
    output_alignment(align1, align2, tag1, tag2, out)


def get_job_id(job_id, total_job, genome_num):
    def get_id(job, genome_num):
        """Get the corresponding alignment index."""
        for i in range(genome_num):
            count = genome_num - i
            if job > count:
                job = job - count
            else:
                j = i + job - 1
                return (i, j)

    def get_job(job_id, total_job, genome_num):
        """Get the list of jobs to run."""
        i = 0
        align_num = genome_num * (genome_num - 1) / 2 + genome_num
        while i * total_job + job_id <= align_num:
            yield i * total_job + job_id
            i += 1

    for job in get_job(job_id, total_job, genome_num):
        yield get_id(job, genome_num)


def diamond_check_db(fastadir, fastaid, suffix, outdir):
    """Make database with DIAMOND."""
    diamond_db = os.path.join(outdir, fastaid)
    if not os.access(diamond_db, os.F_OK):
        subprocess.check_call(
            (
                "diamond makedb "
                "--in " + os.path.join(fastadir, fastaid + suffix) + " "
                "-d " + diamond_db + " "
                "-p 1 --quiet"
            ),
            shell=True,
        )
    return diamond_db


def diamond_blastp(fastadir, query, db, suffix, level, outdir):
    """DIAMOND blastp alignment."""
    level = "" if level == "fast" else f"--{level}"

    # make database with DIAMOND
    diamond_db = diamond_check_db(fastadir, db, suffix, outdir)

    blast_out = os.path.join(outdir, f"{query}_{db}.blast")
    if not os.access(blast_out, os.F_OK):
        subprocess.check_call(
            (
                "diamond blastp "
                "-d " + diamond_db + " "  # database
                "-q " + os.path.join(fastadir, query + suffix) + " "  # query like "<dir>/faa/<query>.faa"
                "-o " + blast_out + " "
                "-k 10 -e 0.001 --id 30 --query-cover 70 --subject-cover 70 "
                "-b 10 --dbsize 1000000000 -p 1 "
                "--quiet " + level
            ),
            shell=True,
        )


def cross_diamond(fastadir, ncid1, ncid2, suffix, level, outdir):
    """Run pairwise DIAMOND."""
    diamond_blastp(fastadir, ncid1, ncid2, suffix, level, outdir)
    diamond_blastp(fastadir, ncid2, ncid1, suffix, level, outdir)


@contextmanager
def TemporaryDirectory(**kwargs):
    """Context manager so it's usable with "with" statement."""
    name = tempfile.mkdtemp(**kwargs)
    try:
        yield name
    finally:
        shutil.rmtree(name)


def set_args(parser: argparse.ArgumentParser):
    parser.add_argument('--list', help='List file', required=True)
    parser.add_argument('--tag-column',
                        help=('which column in list file is used '
                              'as the genome tag (default=%(default)i)'),
                        default=1, type=int)
    parser.add_argument('--ptt', help='PTT dir', required=True)
    parser.add_argument('--fasta', help='Protein sequence dir', required=True)
    parser.add_argument('--suffix',
                        help='Suffix for sequence files (default=%(default)s)',
                        default='.faa')
    parser.add_argument('--out', help='Output file path', required=True)
    parser.add_argument('--job-id', help='Job id used for parallelization',
                        type=int, default=1)
    parser.add_argument('--total-job',
                        help='Total number of jobs used in parallelization',
                        type=int, default=1)
    parser.add_argument('--ident-cutoff',
                        help='Cutoff for identity (default=%(default)f)',
                        default=30, type=float)
    parser.add_argument('--evalue-cutoff',
                        help='Cutoff for E-value (default=%(default)f)',
                        default=1e-3, type=float)
    parser.add_argument('--cov-cutoff',
                        help=('Cutoff for alignment coverage '
                              '(default=%(default)f)'),
                        default=0.7, type=float)
    parser.add_argument('--sensitive-level',
                        help=('Choose which mode to run diamond blastp '
                              '(default=%(default)s)'),
                        default='sensitive', type=str,
                        choices=['fast', 'mid-sensitive', 'sensitive',
                                 'more-sensitive', 'very-sensitive',
                                 'ultra-sensitive'])
    parser.add_argument('--temp-dir',
                        help=('Directory for storing temporary '
                              'files (default=%(default)s)'),
                        default=None, type=str)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    set_args(parser)

    args = parser.parse_args()

    tag_column = args.tag_column - 1  # count starts from 0 in python
    id_list = pd.read_csv(args.list, header=None, skiprows=1, sep='\t',
                          usecols=[0, tag_column], dtype=np.unicode_)

    # if job_id==total_job==1, run process in every iteration,
    # otherwise only run once when job_id==current_id
    # make up a output folder
    os.makedirs(args.out + '_tmp', exist_ok=True)
    # crate the corresponding file for each job, fail if file exists
    with open(args.out + '_tmp/' + str(args.job_id), 'x') as out:
        pass

    # initiate the root logger
    logging.basicConfig(
        format='%(asctime)s\t%(levelname)s:%(message)s',
        level=logging.INFO)
    # temporary output file on local node
    with TemporaryDirectory(dir=args.temp_dir) as tempdirname:
        with open(os.path.join(tempdirname, 'temp'), 'w+') as tmp:
            progress = 0
            logging.info('{:.0f} in total...'.format(
                len(id_list) * (len(id_list) + 1) / 2 // args.total_job))
            sys.stdout.flush()

            for i, j in get_job_id(args.job_id, args.total_job, len(id_list)):
                ncid1, tag1 = id_list.loc[i, [0, tag_column]]
                ncid2, tag2 = id_list.loc[j, [0, tag_column]]
                # run pairwise alignment with DIAMOND
                cross_diamond(args.fasta, ncid1, ncid2,
                              args.suffix, args.sensitive_level, tempdirname)
                filter_alignment(ncid1, ncid2, tag1, tag2, tempdirname,
                                 args.ptt, tmp, args.ident_cutoff,
                                 args.evalue_cutoff, args.cov_cutoff)
                if progress % 100 == 0 and progress != 0:  # report progress
                    logging.info('{} are done...'.format(progress))
                    sys.stdout.flush()
                    tmp.seek(0)  # rewind tmp file to the begining
                    with open(
                            args.out + '_tmp/' + str(args.job_id), 'a') as out:
                        while True:
                            inblock = tmp.read(1024 * 1024)
                            if not inblock:
                                break
                            out.write(inblock)
                    tmp.seek(0)
                    tmp.truncate(0)
                progress += 1
            tmp.seek(0)  # rewind tmp file to the begining
            with open(args.out + '_tmp/' + str(args.job_id), 'a') as out:
                while True:
                    inblock = tmp.read(1024 * 1024)
                    if not inblock:
                        break
                    out.write(inblock)
