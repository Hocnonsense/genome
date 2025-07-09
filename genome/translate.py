# -*- coding: utf-8 -*-
"""
* @Date: 2025-07-04 08:22:24
* @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
* @LastEditTime: 2025-07-09 19:37:50
* @FilePath: /genome/genome/translate.py
* @Description:
"""
# """

import re

from Bio.Data import IUPACData
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord


# region TranslExcept
class TranslExcept:
    PATTERN = re.compile(
        r"\(pos:([0-9]+)\.\.([0-9]+),aa:([A-Za-z]+)\)|"
        r"\(pos:complement\(([0-9]+)\.\.([0-9]+)\),aa:([A-Za-z]+)\)"
    )

    def __new__(cls, text: "str|TranslExcept"):
        """
        Return an existing TranslExcept instance if provided, or create a new instance for a transl_except string.
        
        Parameters:
            text (str or TranslExcept): A transl_except annotation string or an existing TranslExcept instance.
        
        Returns:
            TranslExcept: The provided instance or a new TranslExcept object.
        """
        if isinstance(text, cls):
            return text
        return super().__new__(cls)

    def __init__(self, text: str):
        """
        Initialize a TranslExcept instance by parsing a transl_except annotation string.
        
        Parameters:
            text (str): The transl_except annotation string to parse, specifying the codon position, strand, and amino acid.
        
        Raises:
            ValueError: If the annotation string cannot be parsed.
            NotImplementedError: If the annotation format is not supported.
        """
        self.text = text
        _region, self.strand = self.parse_transl_except(text)
        self.start = int(_region[0])
        self.end = int(_region[1])
        self.aa: str = _region[2]

    def __repr__(self):
        """
        Return a string representation of the translation exception, indicating position, strand orientation, and amino acid.
        
        Returns:
            str: String in standard or complement format based on the strand.
        """
        if self.strand == 1:
            return f"(pos:{self.start}..{self.end},aa:{self.aa})"
        elif self.strand == -1:
            return f"(pos:complement({self.start}..{self.end}),aa:{self.aa})"
        raise NotImplementedError(f"Unknown strand {self.strand}")

    @classmethod
    def parse_transl_except(cls, text: str):
        """
        Parses a transl_except annotation string to extract codon position and amino acid, determining strand orientation.
        
        Parameters:
            text (str): The transl_except annotation string to parse.
        
        Returns:
            tuple: A tuple containing the extracted position and amino acid, along with the strand direction (1 for forward, -1 for complement).
        
        Raises:
            ValueError: If the annotation does not match the expected pattern.
            NotImplementedError: If the pattern is unexpected or ambiguous.
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
        Return the zero-based codon index and amino acid for this translation exception relative to a given CDS location.
        
        Parameters:
            position (SimpleLocation): The CDS location to which the exception is relative.
            partial (bool): Whether the CDS is partial at the 5' end (affects codon indexing).
        
        Returns:
            tuple: (codon_index, amino_acid), where codon_index is the zero-based index of the codon with the exception and amino_acid is its code.
        
        Raises:
            AssertionError: If the strand of the exception and the CDS location do not match.
            NotImplementedError: If the strand is not +1 or -1.
        """
        assert self.strand == position.strand
        if self.strand == 1:
            return (self.start - 1 - position.start - partial) // 3, self.aa  # type: ignore[reportOperatorIssue]
        elif self.strand == -1:
            return (position.end - partial - self.end) // 3, self.aa  # type: ignore[reportOperatorIssue]
        raise NotImplementedError

    @classmethod
    def to_str(cls, transl_except: str, location: SimpleLocation, partial=False):
        """
        Convert a transl_except annotation string into a semicolon-separated string of codon indices and amino acids relative to a given location.
        
        Parameters:
            transl_except (str): The transl_except annotation string.
            location (SimpleLocation): The genomic location to which codon indices are relative.
            partial (bool): Whether the CDS is partial at the 5' end.
        
        Returns:
            str: A semicolon-separated string where each entry is in the format 'codon_index@amino_acid'.
        """
        return ";".join(
            "@".join(str(i) for i in cls(ia).index(location, partial))
            for ia in transl_except
        )

    @classmethod
    def un_str(cls, text: str | SeqRecord):
        """
        Yields codon index and amino acid pairs from a transl_except annotation string or SeqRecord.
        
        If a SeqRecord is provided, extracts the 'transl_except' annotation and parses it. Each entry is expected in the format 'index@amino_acid', separated by semicolons.
        
        Parameters:
            text (str or SeqRecord): A transl_except annotation string or a SeqRecord containing such an annotation.
        
        Yields:
            tuple[int, str]: Pairs of codon index and amino acid code.
        """
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
        Replaces amino acids in a protein sequence at positions specified by translational exceptions.
        
        Parameters:
            aa (str): The original protein sequence.
            text (str): A transl_except annotation string indicating codon indices and amino acids.
            how (str, optional): Placeholder for future modification modes (currently unused).
        
        Returns:
            str: The modified protein sequence with specified amino acids replaced according to transl_except annotations.
        """
        for j, ea in cls.un_str(text):
            if ea in IUPACData.protein_letters_3to1_extended:
                ea = IUPACData.protein_letters_3to1_extended[ea]
            if len(ea) == 1:
                aa = aa[:j] + ea + aa[j + 1 :]
        return aa


# endregion TranslExcept


def update_cds_annotations(fet: SeqFeature):
    """
    Extracts and returns key CDS annotation values from a SeqFeature as a dictionary.
    
    The returned dictionary includes the translation table, reading frame, translational exceptions (if present), and partial status, normalized for downstream translation processing.
    
    Returns:
        annot (dict): Dictionary with keys 'transl_table', 'frame', 'transl_except' (if present), and 'partial'.
    """
    annot: dict[str, str | int] = {}
    annot["transl_table"] = fet.qualifiers.get("transl_table", ["Standard"])[0]
    annot["frame"] = check_frame(fet.qualifiers)
    if "transl_except" in fet.qualifiers:
        annot["transl_except"] = check_transl_except(fet)
    annot["partial"] = fet.qualifiers.get("partial", ["00"])[0]
    return annot


def check_frame(annot: dict):
    """
    Extracts and normalizes the reading frame from annotation data.
    
    Returns:
        int: The reading frame as an integer modulo 3, or 0 if not specified or invalid.
    """
    frame_ = annot.get("frame", [0])
    frame_str = frame_[0] if isinstance(frame_, list) else frame_
    return int(frame_str) % 3 if frame_str and frame_str != "." else 0


def check_transl_except(fet: SeqFeature):
    """
    Extracts and encodes translational exceptions from a CDS feature as a string relative to its location.
    
    Returns:
        A semicolon-separated string representing codon indices and amino acids for each translational exception, or an empty string if none are present.
    """
    if "transl_except" in fet.qualifiers:
        frame = check_frame(fet.qualifiers)
        assert isinstance(fet.location, SimpleLocation)
        return TranslExcept.to_str(
            fet.qualifiers["transl_except"], fet.location, partial=frame > 0
        )
    return ""


def translate(rec: SeqRecord, fet: SeqFeature | None = None, auto_fix=True):
    """
    Translate a nucleotide sequence record into a protein sequence, applying CDS annotations and translation exceptions.
    
    If a CDS feature is provided, uses its qualifiers (translation table, frame, partial status, translational exceptions) to guide translation. Otherwise, uses the record's annotations. Applies frame offset, translation table, and modifies the resulting protein sequence for partial CDS and translational exceptions if `auto_fix` is enabled.
    
    Parameters:
        rec (SeqRecord): The nucleotide sequence record to translate.
        fet (SeqFeature, optional): The CDS feature containing translation qualifiers. If not provided, uses record annotations.
        auto_fix (bool, optional): Whether to apply automatic corrections for partial CDS and translational exceptions. Defaults to True.
    
    Returns:
        SeqRecord: The translated protein sequence record with updated annotations.
    """
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
        if transl_except := check_transl_except(fet):
            annotations["transl_except"] = transl_except
    seq = rec[frame:].translate(
        table=str(annotations.get("transl_table", "Standard")),
        id=True,
        name=True,
        description=True,
        dbxrefs=True,
    )
    seq.annotations.update(annotations)
    if auto_fix:
        if annotations.get("partial", frame) in {0, "00"}:
            # here, if partial is "true" in refseq database, will not fix
            # elif frame is not 0, will not fix
            seq.seq = "M" + seq.seq[1:]
        if "transl_except" in annotations:
            seq.seq = TranslExcept.modify(seq.seq, str(annotations["transl_except"]))
    return seq
