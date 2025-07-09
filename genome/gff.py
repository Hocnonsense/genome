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
    """
    Opens a file or returns a text stream for reading, supporting both plain and gzipped files.
    
    Parameters:
        data: A file path (string or Path) or an existing text stream.
    
    Returns:
        A readable text stream (TextIO) for the provided file or stream.
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
    Writes an iterable of SeqRecord objects to a GFF file.
    
    Parameters:
    	recs (Iterable[SeqRecord]): The sequence records to write.
    	out_handle (PathLike | TextIO): Output file path or writable file-like object. Supports gzip compression if the filename ends with '.gz'.
    	include_fasta (bool): If True, includes FASTA sequences in the output.
    
    Returns:
    	The result of GFFOutput.write, typically None.
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
        """
        Iterates over lines in a GFF file, yielding parsed features and detecting the start of the FASTA section.
        
        Yields:
            Parsed feature objects from valid GFF lines until the FASTA section is encountered. Updates the `fasta_start_pointer` attribute to mark the start of the FASTA section when detected.
        """
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
        """
        Parses and yields FASTA sequences from the point in the file where the FASTA section begins.
        
        Yields:
            SeqRecord: Each sequence record parsed from the FASTA section.
        """
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
        """
        Registers a gene ID inference function under its name in the class-level registry.
        
        Returns:
            The original function, unmodified.
        """
        cls.funcs[fn.__name__] = fn
        return fn

    @classmethod
    def get(cls, fn) -> FUNC_TYPE:
        """
        Retrieve a registered gene ID inference function by name or return the function if already provided.
        
        Parameters:
        	fn: The name of a registered inference function (str) or a callable function.
        
        Returns:
        	A callable gene ID inference function.
        """
        if isinstance(fn, str):
            return cls.funcs[fn]
        return fn


PRODIGAL_ID_PATERN = re.compile(r"^\d+_(\d+)$")


@InferGeneId.rec
def infer_gene_id(rec_id: str | None, fet: SeqFeature):
    """
    Infers a gene ID from a SeqFeature using the "Name" or "gene_group" qualifiers, or falls back to the feature's ID.
    
    Parameters:
    	rec_id (str | None): The record ID, used as a prefix if the "gene_group" qualifier is present.
    	fet (SeqFeature): The feature from which to infer the gene ID.
    
    Returns:
    	str: The inferred gene ID.
    """
    if fet.qualifiers and (name := fet.qualifiers.get("Name")) and name[0]:
        return name[0]
    if (group := fet.qualifiers.get("gene_group")) and group[0]:
        return f"{rec_id}_{group[0]}"
    return fet.id


@InferGeneId.rec
def infer_prodigal_gene_id(rec_id: str | None, fet: SeqFeature):
    """
    Generate a gene ID for a Prodigal-annotated feature by appending the numeric suffix from the feature's ID to the record ID.
    
    Parameters:
    	rec_id (str | None): The sequence record ID to prefix.
    	fet (SeqFeature): The feature whose ID contains the numeric suffix.
    
    Returns:
    	str: The constructed gene ID in the format "{rec_id}_{number}".
    """
    return f"{rec_id}_" + str(int(fet.id.rsplit("_", 1)[1]))


@InferGeneId.rec
def infer_refseq_gene_id(rec_id: str | None, fet: SeqFeature):
    """
    Infers a gene ID from a RefSeq-style feature by extracting the suffix after the last hyphen in the feature's ID.
    
    Parameters:
    	rec_id (str | None): The record ID (unused).
    	fet (SeqFeature): The feature from which to extract the gene ID.
    
    Returns:
    	str: The extracted gene ID suffix.
    """
    return fet.id.rsplit("-", 1)[1]


@InferGeneId.rec
def feat_id(rec_id: str | None, fet: SeqFeature):
    # assert fet.id.startswith(rec_id)
    """
    Returns the feature ID from a SeqFeature object.
    
    Parameters:
    	rec_id (str | None): The record ID, not used in this function.
    	fet (SeqFeature): The feature from which to extract the ID.
    
    Returns:
    	str: The ID of the provided feature.
    """
    return fet.id


# endregion InferGeneId


_biopython_strand = {"+": 1, "-": -1, ".": 0}


def to_seqfeature(feature: gffutils.feature.Feature):
    """
    Convert a gffutils Feature object to a Biopython SeqFeature.
    
    Transfers GFF fields (`source`, `score`, `seqid`, `frame`) and all attributes into the SeqFeature's qualifiers. Coordinates are converted from 1-based (GFF) to 0-based (Biopython), and strand is mapped to Biopython's convention.
    
    Parameters:
        feature: A gffutils Feature object to convert.
    
    Returns:
        SeqFeature: The corresponding Biopython SeqFeature with mapped qualifiers and location.
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
    """
    Create a `Parse` object for reading and processing GFF files, optionally with an associated FASTA file.
    
    Parameters:
        gff: Path to the GFF file.
        fa: Optional path to a FASTA file containing sequence data.
    
    Returns:
        Parse: An instance configured to parse the provided GFF (and FASTA, if given) files.
    """
    if fa is None:
        return Parse(gff)
    return Parse(fa, gff)


class Parse:
    @overload
    def __init__(self, gff: PathLike, /): """
Initialize the iterator for a GFF file, supporting both plain text and gzipped formats.

Parameters:
	gff (PathLike): Path to the GFF file to be iterated.
"""
...
    @overload
    def __init__(self, fa: PathLike, gff: PathLike, /): ...
    @overload
    def __init__(self, gff_or_fa: PathLike, create_now: bool): ...

    def __init__(self, gff_or_fa: PathLike, create_now: bool | PathLike = True):
        """
        Initialize a Parse object for parsing GFF and optionally FASTA files.
        
        Parameters:
            gff_or_fa: Path to a GFF or FASTA file.
            create_now: If True (default), immediately create a gffutils database from the provided file. If a path, use it as the GFF file for database creation. If False, defer database creation.
        """
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
        """
        Creates or loads a gffutils database from a GFF file iterator or the associated FASTA iterator.
        
        If a GFF iterator is provided, it is used to build the database; otherwise, the associated FASTA iterator is used. Supports in-memory or file-based databases and custom merge strategies. Issues a warning if no FASTA sequences are found when expected. Handles empty input files gracefully.
        """
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
        Yields SeqRecord objects parsed from the GFF file, each annotated with features from the gffutils database if available.
        
        If the GFF file does not contain embedded FASTA sequences, yields empty SeqRecords with IDs from the database.
        Otherwise, yields SeqRecords parsed from the FASTA section, with features attached from the database when present.
        
        Parameters:
            limit_info (str | None): Optional parameter for future extension; currently unused.
        
        Returns:
            Generator[SeqRecord, None, None]: An iterator over SeqRecord objects with features.
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
        """
        Extracts features of a specified type from parsed sequence records.
        
        Parameters:
            fet_type (str): The feature type to extract (default is "CDS").
            call_gene_id (InferGeneId.FUNC_TYPE | str): Function or name for inferring gene IDs.
            translate (bool): Whether to translate coding sequences to amino acids.
            min_aa_length (int): Minimum amino acid length for extracted features.
            auto_fix (bool): Whether to automatically fix translation issues.
        
        Returns:
            Generator[SeqRecord, None, None]: Yields extracted feature sequences as SeqRecord objects.
        """
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
    """
    Convert an iterable of SeqRecord objects into a dictionary keyed by unique sequence IDs.
    
    If multiple SeqRecords share the same ID but have different sequences, appends a numeric suffix to duplicate IDs to ensure uniqueness. Only one entry is kept per unique sequence per ID.
    
    Parameters:
        seqs: An iterable of SeqRecord objects to be indexed.
    
    Returns:
        A dictionary mapping unique sequence IDs to SeqRecord objects.
    """
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
    Extracts feature sequences of a specified type from an iterable of SeqRecords, with optional translation and filtering.
    
    Features of the given type (default "CDS") are extracted from each SeqRecord. For CDS features, sequences shorter than the minimum amino acid length (default 33) are skipped. If translation is enabled, extracted CDS sequences are translated to protein. Gene IDs are inferred using the provided function or function name.
    
    Parameters:
        fet_type (str): The feature type to extract (e.g., "CDS").
        call_gene_id (InferGeneId.FUNC_TYPE | str): Function or registered name for inferring gene IDs.
        translate (bool): If True and fet_type is "CDS", translates nucleotide sequences to protein.
        min_aa_length (int): Minimum amino acid length for CDS features to be included.
        auto_fix (bool): If True, attempts to auto-correct translation issues.
    
    Yields:
        SeqRecord: Extracted feature sequences, optionally translated and annotated.
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
    """
    Extracts a subsequence from a SeqRecord corresponding to a given feature, handling circular genomes and updating metadata.
    
    The extracted sequence is assigned a new ID inferred from the feature and record, its description is updated with coordinates and qualifiers, and its features list is replaced with the extracted feature. For CDS features, relevant annotations are updated.
    
    Parameters:
        rec (SeqRecord): The source sequence record.
        fet (SeqFeature): The feature to extract from the record.
        call_gene_id (Callable): Function to infer the gene ID for the extracted sequence.
    
    Returns:
        SeqRecord: The extracted subsequence as a new SeqRecord with updated metadata.
    """
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
    """
    Extracts qualifiers from a description string formatted with a trailing ' # ' segment.
    
    Parameters:
        description (str): A string containing a ' # ' followed by semicolon-separated key=value pairs.
    
    Returns:
        dict: A dictionary mapping qualifier keys to lists of values.
    """
    _, qualifiers_s = description.rsplit(" # ", 1)
    qualifiers = {
        k: v.split(",") for k, v in (i.split("=", 1) for i in qualifiers_s.split(";"))
    }
    return qualifiers
