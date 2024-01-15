# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-12 19:32:50
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-01-16 00:01:10
 * @FilePath: /genome/genome/gff.py
 * @Description:
"""

import gzip
import warnings
from pathlib import Path
from typing import Callable, Generator, Iterable, Optional, TextIO, Union

import gffutils
import gffutils.feature
from Bio import SeqFeature, SeqIO, SeqRecord
from gffutils.exceptions import EmptyInputError
from gffutils.iterators import _FileIterator as _GffutilsFileIterator
from gffutils.iterators import feature_from_line, six

from . import GFFOutput

PathLike = Union[str, Path]
GeneralInput = Union[PathLike, TextIO]


def as_text_io(data: GeneralInput) -> TextIO:
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
    return GFFOutput.write(recs, out_handle, include_fasta)


class _FastaGffFileIterator(_GffutilsFileIterator):
    def open_function(self, data: GeneralInput):
        return as_text_io(data)

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
                    self.fasta_start_pointer = fh.tell()
                    return
                if line.startswith(">"):
                    self.fasta_start_pointer = fh.tell() - len(line)
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

    # quickly init fasta_start_pointer
    def __set_fasta_start_pointer(self, pointer: int):
        self.__fasta_start_pointer = pointer

    def __get_fasta_start_pointer(self):
        if "__fasta_start_pointer" not in self.__dir__():
            for _ in self._custom_iter():
                pass
        return self.__fasta_start_pointer

    fasta_start_pointer = property(__get_fasta_start_pointer, __set_fasta_start_pointer)


def infer_prodigal_gene_id(rec_id: str, fet_id: str):
    return rec_id + "_" + str(int(fet_id.rsplit("_", 1)[1]))


def infer_refseq_gene_id(rec_id: str, fet_id: str):
    return fet_id.rsplit("-", 1)[1]


def infer_trnascan_rna_id(rec_id: str, fet_id: str):
    # assert fet_id.startswith(rec_id)
    return fet_id


_biopython_strand = {
    "+": 1,
    "-": -1,
    ".": 0,
}


def to_seqfeature(feature: gffutils.feature.Feature):
    """
    Converts a gffutils.Feature object to a Bio.SeqFeature object.

    The GFF fields `source`, `score`, `seqid`, and `frame` are stored as
    qualifiers.  GFF `attributes` are also stored as qualifiers.

    Parameters
    ----------
    feature : Feature object, or string
        If string, assume it is a GFF or GTF-format line; otherwise just use
        the provided feature directly.
    """
    qualifiers = {
        "source": [feature.source],
        "score": [feature.score],
        "seqid": [feature.seqid],
        "frame": [feature.frame],
    }
    qualifiers.update(feature.attributes)
    return SeqFeature.SeqFeature(
        # Convert from GFF 1-based to standard Python 0-based indexing used by
        # BioPython
        SeqFeature.SimpleLocation(
            feature.start - 1, feature.stop, strand=_biopython_strand[feature.strand]
        ),
        id=feature.id,
        type=feature.featuretype,
        qualifiers=qualifiers,
    )


class Parse:
    def __init__(self, gff_file: PathLike, create_now: Union[bool, PathLike] = True):
        self.fa = _FastaGffFileIterator(gff_file)
        self.gff_file = Path(gff_file)
        self.db: gffutils.FeatureDB = None
        if isinstance(create_now, (Path, str)):
            self.gff_file = Path(create_now)
            self.create()
            self.fa = _FastaGffFileIterator(create_now)
            self.fa.fasta_start_pointer = 0
        else:
            if create_now is not False:
                self.create()

    def reset_reference(self, refernce_file: PathLike, create_now=False):
        p = Parse(refernce_file, create_now=create_now)
        p.db = self.db
        return p

    def reset_db(self, db_file: PathLike, create_now=True):
        if self.db is None:
            pass
        p = Parse(db_file, create_now=create_now)
        p.fa = self.fa
        return p

    def create(self, dbfn=":memory:", verbose=True, **kwargs):
        try:
            self.db = gffutils.create_db(self.fa, dbfn=dbfn, verbose=verbose, **kwargs)
            if self.fa.fasta_start_pointer == -1:
                warnings.warn(
                    "No sequences found in file. Please check.", RuntimeWarning
                )
        except EmptyInputError:
            self.fa.fasta_start_pointer = 0

    def __call__(
        self, limit_info: Optional[str] = None
    ) -> Generator[SeqRecord.SeqRecord, None, None]:
        """
        High level interface to parse GFF files into SeqRecords and SeqFeatures.
        Add type hints to this function.
        """
        for rec in self.fa.parse_seq():
            if self.db:
                rec.features.extend(
                    (to_seqfeature(fet) for fet in self.db.region(seqid=rec.id))
                )
            yield rec

    def extract(
        self,
        fet_type="CDS",
        translate=True,
        min_aa_length=33,
        call_gene_id: Union[Callable[[str, str], str], str] = infer_prodigal_gene_id,
        auto_fix=True,
    ):
        """
        min_aa_length acturally refers to
          - at least 32 aa complete protein with a complete terminal codon.
          - or at least 33 aa complete protein without terminal codon.
        """
        min_gene_length = min_aa_length * 3
        if fet_type != "CDS":
            assert not translate

        if isinstance(call_gene_id, str):
            call_gene_id = globals()[call_gene_id]

        rec: SeqRecord.SeqRecord = None
        for rec in self():
            fet: SeqFeature.SeqFeature = None
            for fet in rec.features:
                if fet.type != fet_type:
                    continue
                seq: SeqRecord.SeqRecord = fet.extract(rec)
                if fet_type == "CDS" and len(seq) < min_gene_length:
                    continue
                if translate:
                    seq = seq.translate(
                        table=fet.qualifiers.get("transl_table", "Standard")[0]
                    )
                    if (
                        auto_fix
                        and fet.qualifiers.get("partial", "00")[0] == "0"
                        and seq.seq[0] != "M"
                    ):
                        seq.seq = "M" + seq.seq[1:]
                seq.id = call_gene_id(rec.id, fet.id)  # type: ignore
                seq.description = " # ".join(
                    (
                        str(i)
                        for i in (
                            "",
                            fet.location.start + 1,
                            fet.location.end,
                            fet.location.strand,
                            ";".join(
                                f"{k}={','.join(v)}" for k, v in fet.qualifiers.items()
                            ),
                        )
                    )
                ).strip()
                yield seq
