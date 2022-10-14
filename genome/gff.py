# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-12 19:32:50
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-10-14 22:18:09
 * @FilePath: /genome/genome/gff.py
 * @Description:
"""

import gzip
import math
from pathlib import Path
from pickle import dump, load
from typing import (
    Final,
    Generator,
    Iterable,
    Literal,
    NamedTuple,
    TextIO,
    overload,
    Union,
)

import gffutils
from BCBio import GFF as _GFF
from Bio import SeqFeature, SeqIO, SeqRecord
from gffutils.biopython_integration import to_seqfeature
from gffutils.iterators import (
    feature_from_line,
    six,
    _FileIterator as _GffutilsFileIterator,
)
from gffutils.exceptions import EmptyInputError
from numpy import mean


PathLike = Union[str, Path]
GeneralInput = Union[PathLike, TextIO]
GffOutFormat = Literal["faa", "fna"]


def as_test_io(data: GeneralInput) -> TextIO:
    """Allowed input:
    - TextIO
        - anything apply "read" method
    - filename
        - string
        - Path
    """
    if hasattr(data, "read"):
        return data  # type: ignore  # I'm sure this will return a TextIO
    assert not isinstance(data, TextIO)
    if str(data).endswith(".gz"):
        return gzip.open(data, "r")  # type: ignore  # I'm sure this will return a TextIO
    return open(data, "r")


def write(
    recs: Iterable[SeqRecord.SeqRecord],
    out_handle: Union[PathLike, TextIO],
    include_fasta=False,
):
    """
    High level interface to write GFF files into SeqRecords and SeqFeatures.
    Add type hints to this function.
    """
    if str(out_handle).endswith(".gz"):
        out_handle = gzip.open(out_handle, "r")  # type: ignore  # I'm sure this will return a TextIO
    return _GFF.write(recs, out_handle, include_fasta)


class _FastaGffFileIterator(_GffutilsFileIterator):
    def open_function(self, data: GeneralInput):
        return as_test_io(data)

    def _custom_iter(self):
        self.fasta_start_pointer = -1
        self.directives = []
        valid_lines = 0

        # with self.open_function(self.data) as fh:
        fh = self.open_function(self.data)
        if not fh.closed:
            i = 0
            while True:
                line = fh.readline()
                if not line:
                    return
                i += 1

                if isinstance(line, six.binary_type):
                    line = line.decode("utf-8")
                self.current_item = line
                self.current_item_number = i

                if line.startswith("##FASTA"):
                    self.fasta_start_pointer: int = fh.tell()
                    return
                if line.startswith(">"):
                    self.fasta_start_pointer: int = fh.tell() - len(line)
                    return

                line = line.rstrip("\n\r")

                if line.startswith("##"):
                    self._directive_handler(line)
                    continue

                if line.startswith(("#")) or len(line) == 0:
                    continue

                # (If we got here it should be a valid line)
                valid_lines += 1
                yield feature_from_line(line, dialect=self.dialect)
        else:
            raise IOError("data closed unexpectedly")

    def parse_seq(self):
        with self.open_function(self.data) as fh:
            fh.seek(self.fasta_start_pointer)
            rec: SeqRecord.SeqRecord
            for rec in SeqIO.parse(fh, "fasta"):
                yield rec


class Parse:
    def __init__(self, gff_file: PathLike, create_now=True):
        self.gff_file = gff_file
        self.db: gffutils.FeatureDB = None
        if create_now:
            self.create()

    def create(self, dbfn=":memory:", verbose=True, **kwargs):
        self.fasta_gff = _FastaGffFileIterator(self.gff_file)
        try:
            self.db = gffutils.create_db(
                self.fasta_gff, dbfn=dbfn, verbose=verbose, **kwargs
            )
            if self.fasta_gff.fasta_start_pointer == -1:
                raise RuntimeWarning("No sequences found in file. Please check.")
        except EmptyInputError:
            self.fasta_gff.fasta_start_pointer = 0

    def __call__(
        self, limit_info: str = None
    ) -> Generator[SeqRecord.SeqRecord, None, None]:
        """
        High level interface to parse GFF files into SeqRecords and SeqFeatures.
        Add type hints to this function.
        """
        for rec in self.fasta_gff.parse_seq():
            if self.db:
                rec.features.extend(
                    (to_seqfeature(fet) for fet in self.db.region(seqid=rec.id))
                )
            yield rec

    def extract(self, translate=True, min_aa_length=33):
        """
        min_aa_length acturally refer to at least 32 aa complete protein,
        for a complete terminal codon.
        """
        min_gene_length = min_aa_length * 3

        rec: SeqRecord.SeqRecord = None
        for rec in self():
            fet: SeqFeature.SeqFeature = None
            for fet in rec.features:
                if fet.type != "CDS":
                    continue
                seq: SeqRecord.SeqRecord = fet.extract(rec)
                if len(seq) < min_gene_length:
                    continue
                if translate:
                    seq = seq.translate()
                seq.id = rec.id + "_" + str(int(fet.id.rsplit("_", 1)[1]))
                seq.description = " # ".join(
                    (
                        str(i)
                        for i in (
                            "",
                            fet.location.start,
                            fet.location.end,
                            fet.location.strand,
                            ";".join(
                                f"{k}={','.join(v)}" for k, v in fet.qualifiers.items()
                            ),
                        )
                    )
                ).strip()
                yield seq


def calculateN50(seqLens: list[int]):
    thresholdN50 = sum(seqLens) / 2.0

    seqLens.sort(reverse=True)

    testSum = 0
    for seqLen in seqLens:
        testSum += seqLen
        if testSum >= thresholdN50:
            N50 = seqLen
            break

    return N50


class BinStatistic:
    """
    modified from checkm

    - parse cannot be stored in BinStatistic
    """

    class _SeqStat(NamedTuple):
        len: int = 0
        gc_pct: float = 0.0
        gc: int = 0
        at: int = 0
        n: int = 0
        n_aa: int = 0
        len_aa: int = 0

    @overload
    def __init__(self, gff_file: PathLike, min_contig_len=0):
        ...

    @overload
    def __init__(self, last_seq_stats: "BinStatistic", min_contig_len=0):
        ...

    def __init__(self, gff_file, min_contig_len=0):

        if isinstance(gff_file, BinStatistic):
            self.gff_file = gff_file.gff_file
            self.parse = lambda: ()
            self._seq_stats = gff_file._seq_stats
        else:
            self.gff_file = Path(gff_file)
            self.parse = Parse(self.gff_file)
            self._seq_stats: dict[str, "BinStatistic._SeqStat"] = {}

        self.min_contig_len = min_contig_len

        gc, gc_std = self.calculate_gc_std()
        gss = self.calculate_seq_stats()
        coding_len, num_orfs = self.calculate_prot_coding_length()

        self.gc: Final[float] = gc
        self.gc_std: Final[float] = gc_std
        self.bp_size: Final[int] = gss.sum
        self.max_contig_len: Final[int] = gss.max
        self.contigs_num: Final[int] = gss.num
        self.contig_n50: Final[int] = gss.n50
        self.ambiguous_bases_num: Final[int] = gss.numN
        self.coding_density: Final[float] = (
            float(coding_len) / gss.sum if coding_len > 0 and gss.sum else -1
        )
        self.genes_size: Final[int] = num_orfs

    def __repr__(self):
        return (
            f"Statistics of bin <{self.gff_file}>: \n"
            f"GC:              {self.gc:4f}\n"
            f"GC std:          {self.gc_std:4f}\n"
            f"Genome size:     {self.bp_size}\n"
            f"Longest contig:  {self.max_contig_len}\n"
            f"Contig number:   {self.contigs_num}\n"
            f"Contig N50:      {self.contig_n50}\n"
            f"Ambiguous bases: {self.ambiguous_bases_num}\n"
            f"Coding density:  {self.coding_density:4f}\n"
            f"predicted genes: {self.genes_size}"
        )

    def contigs(self) -> Generator[SeqRecord.SeqRecord, None, None]:
        # read contigs
        return (rec for rec in self.parse())

    @property
    def seq_stats(self):
        # calculate once, use every time
        if not self._seq_stats:
            for seq in self.contigs():
                if len(seq.seq) < self.min_contig_len:
                    continue
                a, c, g, t, u, n = (seq.seq.lower().count(base) for base in "acgtun")

                at = a + u + t
                gc = g + c

                if (gc + at) > 0:
                    gcContent = float(gc) / (gc + at)
                else:
                    gcContent = 0.0

                cc: list[int] = [
                    # fet: SeqFeature.SeqFeature
                    int(fet.location.end - fet.location.start)
                    for fet in seq.features
                    if fet.type == "CDS"
                ]

                self._seq_stats[seq.id] = type(self)._SeqStat(
                    len(seq), gcContent, gc, at, n, len(cc), sum(cc)
                )

        return (
            (seq_id, seq_stat)
            for seq_id, seq_stat in self._seq_stats.items()
            if seq_stat.len > self.min_contig_len
        )

    def calculate_gc_std(
        self,
        min_seq_len_gc_std=1000,
    ):
        """
        Calculate fraction of nucleotides that are G or C.
        """
        total_gc = 0
        total_at = 0
        gc_pct_per_seq = []
        for _, seq_stat in self.seq_stats:
            total_gc += seq_stat.gc
            total_at += seq_stat.at

            if seq_stat.len > min_seq_len_gc_std:
                gc_pct_per_seq.append(seq_stat.gc_pct)

        if (total_gc + total_at) > 0:
            gc_pct = float(total_gc) / (total_gc + total_at)
        else:
            gc_pct = 0.0

        var_gc = 0
        if len(gc_pct_per_seq) > 1:
            var_gc = mean(list(map(lambda x: (x - gc_pct) ** 2, gc_pct_per_seq)))

        return gc_pct, math.sqrt(var_gc)

    class _GenomeSeqsStatistics(NamedTuple):
        sum: int
        max: int
        num: int
        n50: int
        numN: int

    def calculate_seq_stats(self):
        """Calculate scaffold length statistics (min length, max length, total length, N50, # contigs)."""
        contig_lens = []
        numAmbiguousBases = 0
        for _, seq_stat in self.seq_stats:
            contig_lens.append(seq_stat.len)

            numAmbiguousBases += seq_stat.n

        contig_N50 = calculateN50(contig_lens)

        return type(self)._GenomeSeqsStatistics(
            sum(contig_lens),
            max(contig_lens),
            len(contig_lens),
            contig_N50,
            numAmbiguousBases,
        )

    def calculate_prot_coding_length(self):
        """Calculate coding density of putative genome bin."""
        len_aa = 0
        n_aa = 0
        for _, seq_stat in self.seq_stats:
            len_aa += seq_stat.len_aa
            n_aa += seq_stat.n_aa

        return len_aa, n_aa

    def dump(self, filename: PathLike = None):
        self.parse = lambda: ()
        if filename is None:
            pickle_filename = Path(f"{self.gff_file}-stat.pkl")
        else:
            pickle_filename = Path(filename)
        pickle_filename.parent.mkdir(parents=True, exist_ok=True)
        with open(pickle_filename, "wb") as po:
            dump(self, po)

    @classmethod
    def load(cls, filename: PathLike):
        with open(filename, "rb") as pi:
            return load(pi)
