# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-12 19:53:55
* @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
* @LastEditTime: 2025-07-04 09:15:15
* @FilePath: /genome/tests/genome/test_gff.py
 * @Description:
__file__ = "tests/genome/test_prokka.py"
"""

from genome import gff
from Bio import SeqIO

from tests import Path, temp_output, test_files, test_temp


@temp_output
def test_gff_extract_protein_fa(test_temp: Path):
    """
    Test extracting protein sequences from a GFF file and writing them to a FASTA file with a minimum amino acid length filter.
    
    This test verifies that protein sequences parsed from the GFF file are correctly filtered by length, sorted by sequence ID, and written in two-line FASTA format to the specified output path.
    """
    gffi = test_files / "binny_contigs_4bins-top10-prodigal.gvmeta.gff"

    test_out = test_temp / "binny_contigs_4bins-top10-prodigal.gvmeta-ge33.faa"
    SeqIO.write(
        sorted(gff.parse(gffi).extract(min_aa_length=33), key=lambda x: str(x.id)),
        test_out,
        "fasta-2line",
    )


@temp_output
def test_gff_reset_reference(test_temp: Path):
    """
    Tests extraction of protein sequences from a GFF file using a reference genome, ensuring sequences of at least 33 amino acids are written to a FASTA file in sorted order.
    """
    gffi = test_files / "binny_contigs_4bins-top10-prodigal.gvmeta.gff"
    genome = test_files / "binny_contigs_4bins.fa"

    test_out = test_temp / "binny_contigs_4bins-top10-prodigal.gvmeta-ge33.faa"
    SeqIO.write(
        sorted(
            gff.parse(gffi, genome).extract(min_aa_length=33),
            key=lambda x: str(x.id),
        ),
        test_out,
        "fasta-2line",
    )


def test_to_dict():
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    seqs = [SeqRecord(Seq(i), "1") for i in ("AAA", "AAAT", "AAAT")]
    assert {k: v.seq for k, v in gff.to_dict(seqs).items()} == {
        "1": Seq("AAA"),
        "1-1": Seq("AAAT"),
    }
