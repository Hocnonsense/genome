# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-12 19:32:50
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-10-15 17:05:56
 * @FilePath: /genome/genome/gff.py
 * @Description:
"""

import gzip
from pathlib import Path
from typing import Generator, Iterable, Literal, TextIO, Union

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
