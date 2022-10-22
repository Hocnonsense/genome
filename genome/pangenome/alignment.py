# -*- coding: utf-8 -*-
"""
 * @Date: 2021-12-04 15:13:17
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-10-13 09:47:46
 * @FilePath: /genome/genome/pangenome/alignment.py
 * @Description:
    To recognize similar genes of given genomes.
    It would include these steps:
    - [x] structure annotation of genome if given genome havn"t been
          annotated.
    - [x] diamond blast of all genes
"""

import argparse
import collections
import itertools
import logging
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime
from typing import IO, Dict, List, Sequence, Tuple

from Bio import SeqIO, AlignIO
from joblib import Parallel, delayed

from . import run_diamond_mcl_one2one as genemap
from . import mcl_to_appearance as appearance

logger = logging.getLogger(__name__)


class utils:
    @staticmethod
    def runsh(
        tokens: Sequence[str],
        env: Dict[str, str] = {},
        trim: bool = True,
        check_call=False,
    ) -> Tuple[str, str]:
        """@description: 在管道中运行 shell 命令
            Run a shell command-line (in token list form) and return its output.
            This does not expand filename patterns or environment variables or do other
            shell processing steps.

            This sets environment variables `PATH` and (if available) `LD_LIBRARY_PATH`.
            Sherlock needs the latter to find libcrypto.so to run `git`.

        @param
            tokens: The command line as a list of string tokens.
            env: add useful environment variables.
            trim: Whether to trim off trailing whitespace. This is useful
                because the subprocess output usually ends with a newline.
            check_call: raise subprocess.CalledProcessError if
                return code is not 0.
        @returns 命令行输出
            stdout, stderr
        """
        if ">" in "".join(tokens):
            raise NotImplementedError(
                "runsh is using subprecess.Popen, cannot recognize popens!"
            )

        environ = {
            "PATH": os.environ["PATH"],
            # 临时设置动态链接库的地址
            "LD_LIBRARY_PATH": os.environ.get("LD_LIBRARY_PATH", ""),
            "PYTHONHASHSEED": "1",
            **env,
        }
        out = subprocess.Popen(
            tokens,
            env=environ,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            shell=True,
        )

        stdout, stderr = (std.decode("utf-8") for std in out.communicate())
        if trim:
            stdout, stderr = stdout.strip("\n"), stderr.strip("\n")
        if out.returncode and check_call:
            raise subprocess.CalledProcessError(out.returncode, tokens, stdout, stderr)
        return stdout, stderr

    def basicConfig(logger_level=None, _logger_level=["INFO"]):
        """This function is used in funtion called by joblib.delayed"""
        if logger_level:
            _logger_level[0] = logger_level
        else:
            logger_level = _logger_level[0]
        logging.basicConfig(level=logger_level)
        return logger

    def prodigal(genome, out_dirname="./faa", out_basename=""):
        out_prefix = os.path.join(out_dirname, out_basename or os.path.basename(genome))
        logger = utils.basicConfig()
        logger.info(f"prodigal: {genome}")
        stdout, stderr = utils.runsh(
            f"prodigal -i {genome}"
            f"         -d {out_prefix}.fna"
            f"         -a {out_prefix}.faa"
            f"         -o {out_prefix}.gff -f gff",
        )
        return out_prefix

    @staticmethod
    def diamond_check_db(fastadir, fastaid, suffix, outdir):
        """Make database with DIAMOND."""
        diamond_db = os.path.join(outdir, fastaid)
        if not os.access(diamond_db, os.F_OK):
            stdout, stderr = utils.runsh(
                "diamond makedb "
                "--in " + os.path.join(fastadir, fastaid + suffix) + " "
                "-d " + diamond_db + " "
                "-p 1 --quiet"
            )
        return diamond_db

    @staticmethod
    def diamond_blastp(fastadir, query, db, suffix, level, outdir):
        """DIAMOND blastp alignment."""
        level = "" if level == "fast" else f"--{level}"

        # make database with DIAMOND
        diamond_db = utils.diamond_check_db(fastadir, db, suffix, outdir)

        blast_out = os.path.join(outdir, f"{query}_{db}.blast")
        if not os.access(blast_out, os.F_OK):
            stdout, stderr = utils.runsh(
                "diamond blastp" "    -d " + diamond_db + " "  # database
                # query like "<dir>/faa/<query>.faa"
                "    -q " + os.path.join(fastadir, query + suffix) + " "
                "    -o " + blast_out + " "
                "    -k 10 -e 0.001 --id 30 --query-cover 70 --subject-cover 70 "
                "    -b 10 --dbsize 1000000000 -p 1 "
                "    " + level,
                check_call=True,
            )

    @staticmethod
    def cross_diamond(fastadir, ncid1, ncid2, suffix, level, outdir):
        """Run pairwise DIAMOND."""
        utils.diamond_blastp(fastadir, ncid1, ncid2, suffix, level, outdir)
        utils.diamond_blastp(fastadir, ncid2, ncid1, suffix, level, outdir)


def faa_2_ptt(genome: str, output_dir: str):
    in_file = os.path.join(output_dir, "faa", genome + ".faa")
    out_file = os.path.join(output_dir, "ptt", genome + ".ptt")
    with open(out_file, "w") as fout:
        for record in SeqIO.parse(in_file, "fasta"):
            print("NA", "NA", len(record.seq), record.id, sep="\t", file=fout)
    return out_file


def diamond_mcl(genomes, output, threads, diamond_args):
    (ident_cutoff, evalue_cutoff, cov_cutoff, sensitive_level) = diamond_args

    def diamond_mcl_1to1(genome1, genome2, tid):
        tout: IO = touts[tid]
        tempdir_t = os.path.join(tempdir, str(tid))
        os.makedirs(tempdir_t, exist_ok=True)

        logger = utils.basicConfig()
        logger.debug(f"diamand: {tid}({tout.name}) : {genome1}, {genome2}")
        for i in range(1):
            try:
                utils.cross_diamond(
                    genome_dir, genome1, genome2, ".faa", sensitive_level, tempdir_t
                )
                logger.debug(f"filter diamand finished: {tid} : {genome1}, {genome2}")
                break
            except subprocess.CalledProcessError as e:
                logger.warning(f"{tid} : {genome1}, {genome2} failed for {i}")
                print(e.returncode, e.cmd, e.output, e.stderr, sep="\n\n")
                e1 = e
        else:
            raise e1
        genemap.filter_alignment(
            genome1,
            genome2,
            genome1,
            genome2,
            tempdir_t,
            os.path.join(output, "ptt"),
            tout,
            ident_cutoff,
            evalue_cutoff,
            cov_cutoff,
        )
        return

    genome_dir = os.path.join(output, "faa")
    with tempfile.TemporaryDirectory(dir=output) as tempdir:
        logger.info(f"diamond_mcl_1to1 at {tempdir}")

        touts = [
            tempfile.SpooledTemporaryFile(mode="w+", prefix=f"genemap_{i}_", dir=output)
            for i in range(threads)
        ]
        [tout.rollover() for tout in touts]  # write to directory

        Parallel(threads, backend="threading")(  #
            delayed(diamond_mcl_1to1)(genome1, genome2, parallel_id)
            for ((genome1, genome2), parallel_id) in zip(
                (
                    (genome1, genome2)
                    for i, genome2 in enumerate(genomes)
                    for genome1 in genomes[0:i]
                ),
                itertools.cycle(range(threads)),
            )
        )

    genemap_path = os.path.join(output, "genemap")
    with open(genemap_path, "w") as fout:
        for tout in touts:
            tout.rollover()
            tout.seek(0)
            while True:
                inblock = tout.read(1024 * 1024)
                if not inblock:
                    break
                fout.write(inblock)
    logger.info(f"genemap wrote to {genemap_path}")
    return genemap_path


def mcl_appearance(genomes, output, threads, mcl_args):
    (slmcl_path,) = mcl_args
    genemap_path = os.path.join(output, "genemap")
    slmcl_out = os.path.join(output, "SL_MCL")
    utils.runsh(
        f"{slmcl_path} "
        f"    -i {genemap_path}"
        f"    -o {slmcl_out}"
        f"    -p {threads}"
    )

    (genome_type, tag_to_genome, genomeid_list) = appearance.parse_list(
        os.path.join(output, "list.tsv"), 0
    )

    def cluster_stat(line: str):
        values = line.strip().split("\t")
        Rep_genome, Rep_gene = values[0].split("@")  # genomeid@geneid
        Rep_genome = tag_to_genome[Rep_genome]
        Gene_list_chromosome: List = []
        Gene_list_plasmid: List = []
        Genome_list = set()
        genome_gene_list = collections.defaultdict(list)
        for gene in values:
            tag, gene = gene.split("@")
            # in case list contains only a subset of genomes
            if tag not in tag_to_genome:
                continue
            genome = tag_to_genome[tag]
            genome_gene_list[genome].append(gene)
            Genome_list.add(genome)
            (
                Gene_list_chromosome
                if genome_type[genome] == "chromosome"
                else Gene_list_plasmid
            ).append(gene)
        # genome list should be non-redundant, use set() to do that
        stat = [
            Rep_gene,
            "-",
            Rep_genome,
            len(Gene_list_chromosome),
            ",".join(Gene_list_chromosome) or "-",
            len(Gene_list_plasmid),
            ",".join(Gene_list_plasmid) or "-",
            len(Genome_list),
            ",".join(Genome_list),
        ]
        gene_count = [len(genome_gene_list[x]) for x in genomeid_list]
        gene_list = [",".join(genome_gene_list[x]) or "-" for x in genomeid_list]
        return stat, gene_count, gene_list

    mcl_file = os.path.join(slmcl_out, "genemap_t3_E2_I1.5.mcl")
    logger.info(f"analyzing {genemap_path}")
    with open(mcl_file) as f:
        stats = Parallel(threads)(delayed(cluster_stat)(line) for line in f)

    appearance_out = os.path.join(output, "appearance")
    os.makedirs(appearance_out)
    appearance_file: str = os.path.join(
        appearance_out, "_".join(("list.tsv", os.path.basename(mcl_file), "{}", ".tsv"))
    )
    with open(appearance_file.format("count"), "w") as cout, open(
        appearance_file.format("list"), "w"
    ) as lout, open(appearance_file.format("list_SC"), "w") as scout:
        title = [
            "Rep_gene",
            "Rep_len",
            "Rep_genome",
            "Gene_cnt_chromosome",
            "Gene_list_chromosome",
            "Gene_cnt_plasmid",
            "Gene_list_plasmid",
            "Genome_cnt",
            "Genome_list",
        ] + genomeid_list
        print(*title, sep="\t", file=cout)
        print(*title, sep="\t", file=lout)
        print(*title, sep="\t", file=scout)
        for stat, gene_count, gene_list in stats:
            print(*stat, *gene_count, sep="\t", file=cout)
            print(*stat, *gene_list, sep="\t", file=lout)
            if stat[3] == len(genomes) and stat[7] == len(genomes):
                print(*stat, *gene_list, sep="\t", file=scout)
    logger.info(f"genemap wrote to {genemap_path}")


def build_sc_tree(genomes, output, threads):
    appearance_file = os.path.join(
        output,
        "appearance",
        "_".join(("list.tsv", "genemap_t3_E2_I1.5.mcl", "list_SC", ".tsv")),
    )
    clusts = []
    with open(appearance_file) as fin:
        title = fin.readline()
        for line in fin:
            clusts.append(line.rstrip().split("\t")[9:])
    genome_seqs = {
        genome: SeqIO.to_dict(
            SeqIO.parse(os.path.join(output, "faa", f"{genome}.faa"), format="fasta")
        )
        for genome in genomes
    }

    msa_out = os.path.join(output, "msa_alignments_full")
    msa_protein_out = os.path.join(msa_out, "scprotein")
    os.makedirs(msa_protein_out)
    logger.info(f"align {len(clusts)} single-copy genes in {msa_protein_out}")

    def muscle_clust_genes(clustid, genes):
        clust_prefix = os.path.join(msa_protein_out, f"{clustid}")
        with open(clust_prefix + ".faa", "w") as fout:
            for genome, gene in zip(genomes, genes):
                record: SeqIO.SeqRecord = genome_seqs[genome][gene]
                record.id = genome  # f"{genome}@{record.id}"
                SeqIO.write(record, fout, "fasta-2line")
        utils.runsh(
            f"muscle "
            f"    -quiet "
            f"    -in {clust_prefix}.faa "
            f"    -out {clust_prefix}.afa"
        )
        return AlignIO.read(f"{clust_prefix}.afa", "fasta")

    cat_seqs = {genome: "" for genome in genomes}
    for align in Parallel(threads, backend="threading")(
        delayed(muscle_clust_genes)(clustid, genes)
        for clustid, genes in enumerate(clusts)
    ):
        for record in align:
            cat_seqs[record.id] += record.seq
    with open(f"{msa_protein_out}.afa", "w") as fout:
        for genome in genomes:
            print(f">{genome}\n{cat_seqs[genome]}", file=fout)

    utils.runsh(f"Gblocks {msa_protein_out}.afa -b5=h")
    logger.info(f"building tree")
    stdout, _ = utils.runsh(f"FastTree -gamma {msa_protein_out}.afa-gb")
    with open(f"{msa_protein_out}.tree", "w") as tout:
        print(stdout, file=tout)
    return f"{msa_protein_out}.tree"


def main(genomes, output, threads, diamond_args, mcl_args):
    os.makedirs(os.path.join(output, "ptt"))
    Parallel(threads)(delayed(faa_2_ptt)(genome, output) for genome in genomes)
    diamond_mcl(genomes, output, threads, diamond_args)
    mcl_appearance(genomes, output, threads, mcl_args)
    build_sc_tree(genomes, output, threads)
    return 0


def get_args() -> Tuple:
    def infer_genomes(input_files):
        try:
            for protein in input_files:
                assert os.path.isfile(protein), protein
        except AssertionError:
            path, genomes = input_files[0], input_files[1:]
            assert os.path.isdir(path), path
            input_files = []
            for genome in genomes:
                protein = os.path.join(path, genome + ".faa")
                assert os.path.isfile(protein)
                input_files.append(protein)
        return input_files

    parser = argparse.ArgumentParser(description=__doc__)
    set_args(parser)
    args = parser.parse_args()
    utils.basicConfig(args.loglevel.upper())

    threads = args.threads
    output = args.output_dir
    if os.path.exists(output):
        logger.warning(f"{output} already exists, cleaning...")
        shutil.rmtree(output)
    genome_dir = os.path.join(output, "faa")
    os.makedirs(genome_dir)

    genomes = infer_genomes(args.input_files)
    if args.gene:
        logger.info(f"input genome are already annotated to genes, just move them.")
        genomes = [shutil.copy(genome, genome_dir) for genome in genomes]
    else:
        genomes = Parallel(threads)(
            delayed(utils.prodigal)(genome, genome_dir) for genome in genomes
        )
        genomes = [genome + ".faa" for genome in genomes]
    genomes = [os.path.basename(genome).rsplit(".faa", 1)[0] for genome in genomes]

    if args.list_chromosome_plasmid:
        shutil.copy(args.list_chromosome_plasmid, os.path.join(output, "list.tsv"))
    else:
        logging.info("Assuming that all sequences are from " "chromosome")
        with open(os.path.join(output, "list.tsv"), "w") as fout:
            print(
                "...",
                *(f"{genome}\tchromosome" for genome in genomes),
                sep="\n",
                file=fout,
            )

    diamond_args = (
        args.ident_cutoff,
        args.evalue_cutoff,
        args.cov_cutoff,
        args.sensitive_level,
    )
    mcl_args = (args.slmcl_path,)

    return genomes, output, threads, diamond_args, mcl_args


def set_args_mcl(parser: argparse.ArgumentParser):
    parser.add_argument(
        "--slmcl-path",
        default=os.path.join(sys.path[0], "..", "slmcl_v0.2.2_linux64"),
        help="path to slmcl",
    )


def set_args_diamond(parser: argparse.ArgumentParser):
    parser.add_argument(
        "--ident-cutoff",
        help="Cutoff for identity (default=%(default)f)",
        default=30,
        type=float,
    )
    parser.add_argument(
        "--evalue-cutoff",
        help="Cutoff for E-value (default=%(default)f)",
        default=1e-3,
        type=float,
    )
    parser.add_argument(
        "--cov-cutoff",
        help=("Cutoff for alignment coverage " "(default=%(default)f)"),
        default=0.7,
        type=float,
    )
    parser.add_argument(
        "--sensitive-level",
        help=("Choose which mode to run diamond blastp " "(default=%(default)s)"),
        default="sensitive",
        type=str,
        choices=[
            "fast",
            "mid-sensitive",
            "sensitive",
            "more-sensitive",
            "very-sensitive",
            "ultra-sensitive",
        ],
    )


def set_args(parser: argparse.ArgumentParser):
    parser.add_argument(
        "--loglevel", default="INFO", type=str, help="set level of logger"
    )
    parser.add_argument(
        "-t", "--threads", default=1, type=int, help="maximum number of threads"
    )
    parser.add_argument(
        "--gene",
        action="store_true",
        help="input file is proteins of "
        "already-annotated-genes "
        "instead of raw genome.fna. "
        "prodigal will be skipped.",
    )
    parser.add_argument(
        "-i",
        "--input-files",
        nargs="+",
        type=str,
        required=True,
        help="path to fasta files. "
        "If only a subset of genome are "
        "selected, <path-to-genome>/ "
        "<basemane1> <basemane2> "
        "is recommended",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        default="./alignment",
        help="Output dir. The dir will be " "overwritten if it already exists.",
    )

    parser.add_argument(
        "--list-chromosome-plasmid",
        default="",
        type=str,
        help="Output dir. The dir will be " "overwritten if it already exists.",
    )

    set_args_diamond(parser)
    set_args_mcl(parser)


def run():
    args = get_args()

    now = datetime.now()
    logger.warning(">>> job start at " + now.strftime("%Y-%m-%d %H:%M:%S"))
    state = main(*args)
    logger.warning(">>> job run time: " + str(datetime.now() - now))
    if state == 0:
        logger.info("success!")


if __name__ == "__main__":
    run()
