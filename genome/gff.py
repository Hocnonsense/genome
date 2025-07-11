# -*- coding: utf-8 -*-
"""
* @Date: 2022-10-12 19:32:50
* @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
* @LastEditTime: 2025-07-04 08:36:39
* @FilePath: /genome/genome/gff.py
* @Description:
"""

import gzip
import re
import warnings
from pathlib import Path
from typing import Callable, Generator, Iterable, TextIO, overload

import gffutils
import gffutils.feature
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
from gffutils.exceptions import EmptyInputError
from gffutils.iterators import _FileIterator as _GffutilsFileIterator
from gffutils.iterators import feature_from_line

from . import GFFOutput
from .translate import translate as _translate
from .translate import update_cds_annotations

PathLike = str | Path


def as_text_io(data: PathLike | TextIO) -> TextIO:
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
    recs: Iterable[SeqRecord],
    out_handle: PathLike | TextIO,
    include_fasta=False,
):
    """
    High level interface to write GFF files into SeqRecords and SeqFeatures.
    Add type hints to this function.
    """
    if not hasattr(out_handle, "write"):
        if str(out_handle).endswith(".gz"):
            out_handle = gzip.open(out_handle, "w")  # type: ignore[assignment]
        else:
            out_handle = open(out_handle, "w")  # type: ignore[assignment]
    return GFFOutput.write(recs, out_handle, include_fasta)


class _FastaGffFileIterator(_GffutilsFileIterator):
    open_function = staticmethod(as_text_io)  # type: ignore[reportAssignmentType]

    def _custom_iter(self):
        self.fasta_start_pointer = -1
        self.directives = []
        valid_lines = 0

        # with self.open_function(self.data) as fh:
        fh = self.open_function(self.data)
        if not fh.closed:
            i = 0
            while True:
                line: str | bytes = fh.readline()
                if not line:
                    return
                i += 1

                if isinstance(line, bytes):
                    line = line.decode("utf-8")
                self.current_item = line
                self.current_item_number = i

                assert isinstance(line, str)
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

                if line.startswith("#") or len(line) == 0:
                    continue

                # (If we got here it should be a valid line)
                valid_lines += 1
                yield feature_from_line(line, dialect=self.dialect)
        else:
            raise IOError("data closed unexpectedly")

    def parse_seq(self):
        with self.open_function(self.data) as fh:
            fh.seek(self.fasta_start_pointer)
            rec: SeqRecord
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


# region InferGeneId
class InferGeneId:
    FUNC_TYPE = Callable[[str | None, SeqFeature], str]
    funcs: dict[str, FUNC_TYPE] = {}

    @classmethod
    def rec(cls, fn: FUNC_TYPE):
        cls.funcs[fn.__name__] = fn
        return fn

    @classmethod
    def get(cls, fn) -> FUNC_TYPE:
        if isinstance(fn, str):
            return cls.funcs[fn]
        return fn


PRODIGAL_ID_PATERN = re.compile(r"^\d+_(\d+)$")


@InferGeneId.rec
def infer_gene_id(rec_id: str | None, fet: SeqFeature):
    if fet.qualifiers and (name := fet.qualifiers.get("Name")) and name[0]:
        return name[0]
    if (group := fet.qualifiers.get("gene_group")) and group[0]:
        return f"{rec_id}_{group[0]}"
    return fet.id


@InferGeneId.rec
def infer_prodigal_gene_id(rec_id: str | None, fet: SeqFeature):
    return f"{rec_id}_" + str(int(fet.id.rsplit("_", 1)[1]))


@InferGeneId.rec
def infer_refseq_gene_id(rec_id: str | None, fet: SeqFeature):
    return fet.id.rsplit("-", 1)[1]


@InferGeneId.rec
def feat_id(rec_id: str | None, fet: SeqFeature):
    # assert fet.id.startswith(rec_id)
    return fet.id


# endregion InferGeneId


_biopython_strand = {"+": 1, "-": -1, ".": 0}


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
    start, stop = sorted((feature.start or 0, feature.end or 0))
    return SeqFeature(
        # Convert from GFF 1-based to standard Python 0-based indexing used by
        # BioPython
        SimpleLocation(start - 1, stop, strand=_biopython_strand[feature.strand]),
        id=feature.id or "",
        type=feature.featuretype,
        qualifiers=qualifiers,
    )


def parse(gff: PathLike, fa: PathLike | None = None):
    if fa is None:
        return Parse(gff)
    return Parse(fa, gff)


class Parse:
    @overload
    def __init__(self, gff: PathLike, /): ...
    @overload
    def __init__(self, fa: PathLike, gff: PathLike, /): ...
    @overload
    def __init__(self, gff_or_fa: PathLike, create_now: bool): ...

    def __init__(self, gff_or_fa: PathLike, create_now: bool | PathLike = True):
        self.fa = _FastaGffFileIterator(gff_or_fa)
        self.db: gffutils.FeatureDB = None  # type: ignore[assignment]
        if isinstance(create_now, (Path, str)):
            gff_file = create_now
            self.create(verbose=False, expect_fa=False)  # init pointer of fa
            self.create(_FastaGffFileIterator(gff_file), verbose=False)
            self.gff_file = Path(gff_file)
        else:
            self.gff_file = Path(gff_or_fa)
            if create_now is not False:
                self.create(verbose=False)

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

    def create(
        self,
        gff: _FastaGffFileIterator | None = None,
        dbfn=":memory:",
        verbose=True,
        expect_fa=True,
        merge_strategy="create_unique",
        **kwargs,
    ):
        if gff is not None:
            self.db = gffutils.create_db(
                gff, dbfn=dbfn, verbose=verbose, merge_strategy=merge_strategy, **kwargs
            )
        else:
            try:
                self.db = gffutils.create_db(
                    gff or self.fa,
                    dbfn=dbfn,
                    verbose=verbose,
                    merge_strategy=merge_strategy,
                    **kwargs,
                )
                if expect_fa and self.fa.fasta_start_pointer == -1:
                    warnings.warn(
                        "No sequences found in file. Please check.", RuntimeWarning
                    )
            except EmptyInputError:
                self.fa.fasta_start_pointer = 0

    def __call__(
        self, limit_info: str | None = None
    ) -> Generator[SeqRecord, None, None]:
        """
        High level interface to parse GFF files into SeqRecords and SeqFeatures.
        Add type hints to this function.
        """
        if self.fa.fasta_start_pointer == -1:
            if self.db is None:
                return
            seqs = (SeqRecord(None, id=rec_id) for rec_id in self.db.seqids())
        else:
            seqs = self.fa.parse_seq()
        for rec in seqs:
            if self.db:
                rec.features.extend(
                    (to_seqfeature(fet) for fet in self.db.region(seqid=rec.id))
                )
            yield rec

    def extract(
        self,
        fet_type="CDS",
        call_gene_id: InferGeneId.FUNC_TYPE | str = infer_gene_id,
        translate=True,
        min_aa_length=33,
        auto_fix=True,
    ):
        return extract(
            self(),
            fet_type=fet_type,
            call_gene_id=call_gene_id,
            translate=translate,
            min_aa_length=min_aa_length,
            auto_fix=auto_fix,
        )


def to_dict(
    seqs: Iterable[SeqRecord],
):
    seqd: dict[str, list[SeqRecord]] = {}
    for seq in seqs:
        if not isinstance(seq.id, str):
            seq.id = str(seq.id)
        seql = seqd.setdefault(seq.id, [])
        if not any(i for i in seql if i.seq == seq.seq):
            seql.append(seq)
    seqd1 = {}
    for seql in seqd.values():
        seqd1[seql[0].id] = seql[0]
        for i, seq in enumerate(seql[1:], 1):
            seq.id = f"{seq.id}-{i}"
            seqd1[seq.id] = seq
    return seqd1


def extract(
    seq_record: Iterable[SeqRecord],
    fet_type="CDS",
    call_gene_id: InferGeneId.FUNC_TYPE | str = infer_gene_id,
    translate=True,
    min_aa_length=33,
    auto_fix=True,
):
    """
    `min_aa_length=33` actually refers to
      - at least 32 aa complete protein with a complete terminal codon.
      - or at least 33 aa complete protein without terminal codon.
    """
    min_gene_length = int(min_aa_length) * 3
    # if fet_type != "CDS":
    #    assert not translate
    call_gene_id = InferGeneId.get(call_gene_id)

    rec: SeqRecord
    for rec in seq_record:
        fet: SeqFeature
        for fet in rec.features:
            if fet.type != fet_type or fet.location is None:
                continue
            # region extract seq
            if fet_type == "CDS":
                if len(fet) < min_gene_length:
                    continue
            seq = extract1(rec, fet, call_gene_id)
            if fet_type == "CDS" and translate:
                seq = _translate(seq, fet, auto_fix=auto_fix)
            yield seq


def extract1(
    rec: SeqRecord,
    fet: SeqFeature,
    call_gene_id: InferGeneId.FUNC_TYPE = infer_gene_id,
):
    len_rec = len(rec)
    assert fet.location is not None
    start, end = int(fet.location.start), int(fet.location.end)  # type: ignore[reportArgumentType]
    if end >= len_rec:
        # circular genome without overlap
        seq: SeqRecord = fet.extract(rec + rec)  # type: ignore[reportArgumentType]
    else:
        seq = fet.extract(rec)  # type: ignore[reportArgumentType]
    seq.features.clear()  # Drop features from `SeqFeature.extract`
    seq.features.append(fet)
    # endregion extract seq
    seq.id = call_gene_id(rec.id, fet)
    seq.description = " # ".join(
        (
            str(i).strip()
            for i in (
                "",
                start + 1,
                end,
                fet.location.strand,
                ";".join(f"{k}={','.join(v)}" for k, v in fet.qualifiers.items()),
            )
        )
    )
    if fet.type == "CDS":
        seq.annotations |= update_cds_annotations(fet)
    return seq


def recover_qualifiers(description: str):
    _, qualifiers_s = description.rsplit(" # ", 1)
    qualifiers = {
        k: v.split(",") for k, v in (i.split("=", 1) for i in qualifiers_s.split(";"))
    }
    return qualifiers
