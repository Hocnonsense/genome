# -*- coding: utf-8 -*-
"""
* @Date: 2022-10-12 19:32:50
* @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
* @LastEditTime: 2025-04-29 10:22:01
* @FilePath: /genome/genome/gff.py
* @Description:
"""

from functools import wraps
import gzip
import warnings
from pathlib import Path
from typing import Callable, Generator, Iterable, TextIO
import re

import gffutils
import gffutils.feature
from Bio import SeqIO
from Bio.Data import IUPACData
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
from gffutils.exceptions import EmptyInputError
from gffutils.iterators import _FileIterator as _GffutilsFileIterator
from gffutils.iterators import feature_from_line
from matplotlib.pylab import overload

from . import GFFOutput

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
            out_handle = open(out_handle, "w")
    return GFFOutput.write(recs, out_handle, include_fasta)


class _FastaGffFileIterator(_GffutilsFileIterator):
    open_function = staticmethod(as_text_io)

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

                if isinstance(line, bytes):
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
    def rec(self, fn: FUNC_TYPE):
        self.funcs[fn.__name__] = fn
        return fn

    @classmethod
    def get(self, fn) -> FUNC_TYPE:
        if isinstance(fn, str):
            return self.funcs[fn]
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
def infer_trnascan_rna_id(rec_id: str | None, fet: SeqFeature):
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
    start, stop = sorted((feature.start, feature.end))
    return SeqFeature(
        # Convert from GFF 1-based to standard Python 0-based indexing used by
        # BioPython
        SimpleLocation(start - 1, stop, strand=_biopython_strand[feature.strand]),
        id=feature.id,
        type=feature.featuretype,
        qualifiers=qualifiers,
    )


# region TranslExcept
class TranslExcept:
    PATTERN = re.compile(
        r"\(pos:([0-9]+)\.\.([0-9]+),aa:([A-Za-z]+)\)|"
        r"\(pos:complement\(([0-9]+)\.\.([0-9]+)\),aa:([A-Za-z]+)\)"
    )

    def __new__(cls, text: "str|TranslExcept"):
        if isinstance(text, cls):
            return text
        return super().__new__(cls)

    def __init__(self, text: str):
        # GFF entry 是一行 GFF 文件的记录 (字符串)
        self.text = text
        _region, self.strand = self.parse_transl_except(text)
        self.start = int(_region[0])
        self.end = int(_region[1])
        self.aa: str = _region[2]

    def __repr__(self):
        if self.strand == 1:
            return f"(pos:{self.start}..{self.end},aa:{self.aa})"
        elif self.strand == -1:
            return f"(pos:complement({self.start}..{self.end}),aa:{self.aa})"
        raise NotImplementedError(f"Unknown strand {self.strand}")

    @classmethod
    def parse_transl_except(cls, text: str):
        """
        AP024703.1	DDBJ	CDS	3092001	3095066	.	+	0	ID=cds-BCX53216.1;Note=codon on position 197 is selenocysteine opal codon.;transl_except=(pos:3092589..3092591,aa:Sec)
        JACSQK010000001.1	Protein Homology	CDS	324590	327646	.	-	0	ID=cds-MBD7959139.1;transl_except=(pos:complement(327056..327058),aa:Sec)
        """
        matches = cls.PATTERN.findall(text)
        if len(matches) != 1:
            raise ValueError(f"Invalid transl_except: {text}=>{matches}")
        if matches[0][:3] == ("", "", ""):
            return matches[0][3:], -1
        if matches[0][3:] == ("", "", ""):
            return matches[0][:3], 1
        raise NotImplementedError(f"{matches[0][:3]} {matches[0][3:]}")

    def index(self, position: SimpleLocation, partial: bool):
        """
        >>> TranslExcept("(pos:3092589..3092591,aa:Sec)").index(SimpleLocation(3092001, 3095066, 1))
        (196, 'Sec')
        >>> TranslExcept("(pos:complement(327056..327058),aa:Sec)").index(SimpleLocation(324590, 327646, -1))
        (196, 'Sec')
        """
        assert self.strand == position.strand
        if self.strand == 1:
            return (self.start - 1 - position.start - partial) // 3, self.aa
        elif self.strand == -1:
            return (position.end - partial - self.end) // 3, self.aa
        raise NotImplementedError

    @classmethod
    def to_str(cls, transl_except: str, location: SimpleLocation, partial=False):
        return ";".join(
            "@".join(str(i) for i in cls(ia).index(location, partial))
            for ia in transl_except
        )

    @classmethod
    def un_str(cls, text: str | SeqRecord):
        if isinstance(text, SeqRecord):
            text = str(text.annotations.get("transl_except", ""))
        if not text:
            return
        for i in text.split(";"):
            i, a = i.split("@")
            yield int(i), a

    @classmethod
    def modify(cls, aa, text: str, how="X"):
        """
        U = "Sec";  selenocysteine
        O = "Pyl";  pyrrolysine
        """
        for j, ea in cls.un_str(text):
            if ea in IUPACData.protein_letters_3to1_extended:
                ea = IUPACData.protein_letters_3to1_extended[ea]
            if len(ea) == 1:
                aa = aa[:j] + ea + aa[j + 1 :]
        return aa


# endregion TranslExcept


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
        self.db: gffutils.FeatureDB = None
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
        len_rec = len(rec)
        fet: SeqFeature
        for fet in rec.features:
            if fet.type != fet_type or fet.location is None:
                continue
            # region extract seq
            if fet.location.end >= len_rec:
                # circular genome without overlap
                seq: SeqRecord = fet.extract(rec + rec)
            else:
                seq = fet.extract(rec)
            seq.features.clear()  # Drop features from `SeqFeature.extract`
            seq.features.append(fet)
            # endregion extract seq
            if fet_type == "CDS":
                if len(seq) < min_gene_length:
                    continue
                if translate:
                    seq = _translate(seq, fet, auto_fix=auto_fix)
                else:
                    seq.annotations["transl_table"] = fet.qualifiers.get(
                        "transl_table", ["Standard"]
                    )[0]
                    seq.annotations["frame"] = check_frame(fet.qualifiers)
                    if "transl_except" in fet.qualifiers:
                        seq.annotations["transl_except"] = check_transl_except(fet)
                seq.annotations["partial"] = fet.qualifiers.get("partial", ["00"])[0]
            seq.id = call_gene_id(rec.id, fet)
            seq.description = " # ".join(
                (
                    str(i).strip()
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
            )
            yield seq


def check_frame(annot: dict):
    frame_ = annot.get("frame", [0])
    frame_str = frame_[0] if isinstance(frame_, list) else frame_
    return int(frame_str) % 3 if frame_str and frame_str != "." else 0


def check_transl_except(fet: SeqFeature):
    if "transl_except" in fet.qualifiers:
        frame = check_frame(fet.qualifiers)
        assert isinstance(fet.location, SimpleLocation)
        return TranslExcept.to_str(
            fet.qualifiers["transl_except"], fet.location, partial=frame > 0
        )
    return ""


def translate(rec: SeqRecord, fet: SeqFeature | None = None, auto_fix=True):
    if fet is None:
        annotations = rec.annotations
        frame = check_frame(annotations)
    else:
        annotations = {
            k: fet.qualifiers[k][0]
            for k in ("transl_table", "partial")
            if k in fet.qualifiers
        }
        frame = annotations["frame"] = check_frame(fet.qualifiers)
    partial = annotations.get("partial", frame) in {0, "00"}
    if fet is not None and (transl_except := check_transl_except(fet)):
        annotations["transl_except"] = transl_except
    seq = rec[frame:].translate(table=str(annotations.get("transl_table", "Standard")))
    seq.annotations.update(annotations)
    if auto_fix:
        if partial:
            # here, if partial is "true" in refseq database, will not fix
            # elif frame is not 0, will not fix
            seq.seq = "M" + seq.seq[1:]
        if "transl_except" in annotations:
            seq.seq = TranslExcept.modify(seq.seq, str(annotations["transl_except"]))
    return seq


_translate = translate


def recover_qualifiers(description: str):
    _, qualifiers_s = description.rsplit(" # ", 1)
    qualifiers = {
        k: v.split(",") for k, v in (i.split("=", 1) for i in qualifiers_s.split(";"))
    }
    return qualifiers
