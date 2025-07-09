# -*- coding: utf-8 -*-
"""
 * @Author: chapmanb https://github.com/chapmanb
 * @Date: 2023-02-15 commit 8a36af0
 * @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2025-01-06 22:07:49
 * @FilePath: /genome/genome/GFFOutput.py
 * @Description:
    Output Biopython SeqRecords and SeqFeatures to GFF3 format.

    The target format is GFF3, the current GFF standard:
        http://www.sequenceontology.org/gff3.shtml

    https://github.com/chapmanb/bcbb/blob/master/gff/BCBio/GFF/GFFParser.py
"""
# """
from urllib.parse import quote

from Bio import SeqIO


class _IdHandler:
    """Generate IDs for GFF3 Parent/Child relationships where they don't exist."""

    def __init__(self):
        self._prefix = "biopygen"
        self._counter = 1
        self._seen_ids = []

    def _generate_id(self, quals):
        """Generate a unique ID not present in our existing IDs."""
        gen_id = self._get_standard_id(quals)
        if gen_id is None:
            while 1:
                gen_id = "%s%s" % (self._prefix, self._counter)
                if gen_id not in self._seen_ids:
                    break
                self._counter += 1
        return gen_id

    def _get_standard_id(self, quals):
        """Retrieve standardized IDs from other sources like NCBI GenBank.

        This tries to find IDs from known key/values when stored differently
        than GFF3 specifications.
        """
        possible_keys = ["transcript_id", "protein_id"]
        for test_key in possible_keys:
            if test_key in quals:
                cur_id = quals[test_key]
                if isinstance(cur_id, tuple) or isinstance(cur_id, list):
                    return cur_id[0]
                else:
                    return cur_id
        return None

    def update_quals(self, quals, has_children):
        """Update a set of qualifiers, adding an ID if necessary."""
        cur_id = quals.get("ID", None)
        # if we have an ID, record it
        if cur_id:
            if not isinstance(cur_id, list) and not isinstance(cur_id, tuple):
                cur_id = [cur_id]
            for add_id in cur_id:
                self._seen_ids.append(add_id)
        # if we need one and don't have it, create a new one
        elif has_children:
            new_id = self._generate_id(quals)
            self._seen_ids.append(new_id)
            quals["ID"] = [new_id]
        return quals


class GFF3Writer:
    """Write GFF3 files starting with standard Biopython objects."""

    def __init__(self):
        pass

    def write(self, recs, out_handle, include_fasta=False):
        """
        Writes one or more Biopython SeqRecords to the output handle in GFF3 format, optionally including sequence data in FASTA format.
        
        Parameters:
        	recs: A single SeqRecord or an iterable of SeqRecords to be written.
        	out_handle: A writable file-like object where the GFF3 output will be written.
        	include_fasta (bool): If True, appends the corresponding sequence(s) in FASTA format after the GFF3 features.
        """
        id_handler = _IdHandler()
        self._write_header(out_handle)
        fasta_recs = []
        try:
            recs = iter(recs)
        except TypeError:
            recs = [recs]
        for rec in recs:
            seq_len = len(rec.seq) if rec.seq else 0
            self._write_rec(rec, out_handle)
            self._write_annotations(rec.annotations, rec.id, seq_len, out_handle)
            for sf in rec.features:
                sf = self._clean_feature(sf)
                id_handler = self._write_feature(sf, rec.id, out_handle, id_handler)
            if include_fasta and seq_len > 0:
                fasta_recs.append(rec)
        if len(fasta_recs) > 0:
            self._write_fasta(fasta_recs, out_handle)

    def _clean_feature(self, feature):
        quals = {}
        for key, val in feature.qualifiers.items():
            if not isinstance(val, (list, tuple)):
                val = [val]
            val = [str(x) for x in val]
            quals[key] = val
        feature.qualifiers = quals
        # Support for Biopython 1.68 and above, which removed sub_features
        if not hasattr(feature, "sub_features"):
            feature.sub_features = []
        clean_sub = [self._clean_feature(f) for f in feature.sub_features]
        feature.sub_features = clean_sub
        return feature

    def _write_rec(self, rec, out_handle):
        # if we have a SeqRecord, write out optional directive
        """
        Writes the GFF3 sequence-region directive for a SeqRecord if it contains a non-empty sequence.
        
        The directive specifies the region covered by the sequence, using 1-based coordinates.
        """
        if rec.seq and len(rec.seq) > 0:
            out_handle.write("##sequence-region %s 1 %s\n" % (rec.id, len(rec.seq)))

    def _get_phase(self, feature):
        """
        Determine the phase value for a feature in GFF3 format.
        
        For CDS features, calculates the phase from the 'codon_start' qualifier if present; otherwise, uses the 'phase' qualifier or defaults to '.'.
        Returns the phase as a string.
        """
        if "phase" in feature.qualifiers:
            phase = feature.qualifiers["phase"][0]
        elif feature.type == "CDS":
            phase = int(feature.qualifiers.get("codon_start", [1])[0]) - 1
        else:
            phase = "."
        return str(phase)

    def _write_feature(self, feature, rec_id, out_handle, id_handler, parent_id=None):
        """
        Writes a single feature and its subfeatures to the output in GFF3 format, including location, attributes, and parent-child relationships.
        
        Parameters:
            feature: The SeqFeature to write.
            rec_id: The sequence record identifier for the feature.
            out_handle: Output stream to write the GFF3 line.
            id_handler: Handler for generating and tracking unique feature IDs.
            parent_id: Optional parent feature ID for hierarchical relationships.
        
        Returns:
            Updated ID handler reflecting any new IDs assigned during the write process.
        """
        if feature.strand == 1:
            strand = "+"
        elif feature.strand == -1:
            strand = "-"
        else:
            strand = "."
        # remove any standard features from the qualifiers
        quals = feature.qualifiers.copy()
        for std_qual in ["source", "score", "phase"]:
            if std_qual in quals and len(quals[std_qual]) == 1:
                del quals[std_qual]
        # add a link to a parent identifier if it exists
        if parent_id:
            if not "Parent" in quals:
                quals["Parent"] = []
            quals["Parent"].append(parent_id)
        quals = id_handler.update_quals(quals, len(feature.sub_features) > 0)
        if feature.type:
            ftype = feature.type
        else:
            ftype = "sequence_feature"
        parts = [
            str(rec_id),
            feature.qualifiers.get("source", ["feature"])[0],
            ftype,
            str(feature.location.start + 1),  # 1-based indexing
            str(feature.location.end),
            feature.qualifiers.get("score", ["."])[0],
            strand,
            self._get_phase(feature),
            self._format_keyvals(quals),
        ]
        out_handle.write("\t".join(parts) + "\n")
        for sub_feature in feature.sub_features:
            id_handler = self._write_feature(
                sub_feature, rec_id, out_handle, id_handler, quals["ID"][0]
            )
        return id_handler

    def _format_keyvals(self, keyvals):
        format_kvs = []
        for key in sorted(keyvals.keys()):
            values = keyvals[key]
            key = key.strip()
            format_vals = []
            if not isinstance(values, list) or isinstance(values, tuple):
                values = [values]
            for val in values:
                val = quote(str(val).strip(), safe=":/ ")
                if (key and val) and val not in format_vals:
                    format_vals.append(val)
            format_kvs.append("%s=%s" % (key, ",".join(format_vals)))
        return ";".join(format_kvs)

    def _write_annotations(self, anns, rec_id, size, out_handle):
        """Add annotations which refer to an entire sequence."""
        format_anns = self._format_keyvals(anns)
        if format_anns:
            parts = [
                rec_id,
                "annotation",
                "remark",
                "1",
                str(size if size > 1 else 1),
                ".",
                ".",
                ".",
                format_anns,
            ]
            out_handle.write("\t".join(parts) + "\n")

    def _write_header(self, out_handle):
        """Write out standard header directives."""
        out_handle.write("##gff-version 3\n")

    def _write_fasta(self, recs, out_handle):
        """Write sequence records using the ##FASTA directive."""
        out_handle.write("##FASTA\n")
        SeqIO.write(recs, out_handle, "fasta")


def write(recs, out_handle, include_fasta=False):
    """High level interface to write GFF3 files from SeqRecords and SeqFeatures.

    If include_fasta is True, the GFF3 file will include sequence information
    using the ##FASTA directive.
    """
    writer = GFF3Writer()
    return writer.write(recs, out_handle, include_fasta)
