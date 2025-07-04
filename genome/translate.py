# -*- coding: utf-8 -*-
"""
* @Date: 2025-07-04 08:22:24
* @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
* @LastEditTime: 2025-07-04 08:31:34
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
            return (self.start - 1 - position.start - partial) // 3, self.aa  # type: ignore[reportOperatorIssue]
        elif self.strand == -1:
            return (position.end - partial - self.end) // 3, self.aa  # type: ignore[reportOperatorIssue]
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


def update_cds_annotations(fet: SeqFeature):
    annot: dict[str, str | int] = {}
    annot["transl_table"] = fet.qualifiers.get("transl_table", ["Standard"])[0]
    annot["frame"] = check_frame(fet.qualifiers)
    if "transl_except" in fet.qualifiers:
        annot["transl_except"] = check_transl_except(fet)
    annot["partial"] = fet.qualifiers.get("partial", ["00"])[0]
    return annot


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
