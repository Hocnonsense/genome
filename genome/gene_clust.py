# -*- coding: utf-8 -*-
"""
* @Date: 2022-10-15 21:29:41
* @LastEditors: hwrn hwrn.aou@sjtu.edu.cn
* @LastEditTime: 2025-07-01 16:40:24
* @FilePath: /genome/genome/gene_clust.py
* @Description:
"""


import os
import shutil
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import TYPE_CHECKING, Collection, Iterable, Literal, NamedTuple

import pandas as pd
from Bio import SeqIO, SeqRecord

from .pyrule import rules_dir, smk

PathLike = str | Path

GENE_CLUST_SMK = rules_dir / "gene_clust.smk"


class _CluBase:
    @classmethod
    def in_faa(cls, prefix: PathLike = "{prefix}"):
        """
        Return the path to the input FASTA (.faa) file corresponding to the given prefix.
        
        Parameters:
            prefix: The file prefix or base name for the FASTA file.
        
        Returns:
            Path to the input FASTA file with a `.faa` extension.
        """
        return Path(f"{prefix}.faa")

    @classmethod
    def from_in_faa(cls, faa: PathLike):
        """
        Create a new instance from the given input FASTA (.faa) file by inferring the prefix.
        
        Parameters:
            faa: Path to the input FASTA file. Must end with the '.faa' extension.
        
        Returns:
            An instance initialized with file paths derived from the input FASTA file's prefix.
        """
        assert f"{faa}".endswith(".faa")
        return cls.from_prefix(str(faa)[:-4])

    @classmethod
    def from_prefix(cls, prefix: PathLike):
        """
        Create an instance from a given file prefix by formatting all output file paths accordingly.
        
        Parameters:
        	prefix: The base path or prefix used to construct output file paths.
        
        Returns:
        	An instance with all file paths set based on the provided prefix.
        """
        return cls(*(Path(f"{i}".format(prefix=prefix)) for i in cls()))

    @classmethod
    def from_aout(cls, aout: PathLike):
        """
        Create an instance by inferring the file prefix from an output file path.
        
        Attempts to match the provided output file path against known suffix patterns to extract the prefix. Raises a KeyError if the suffix cannot be determined.
        """
        aout_ = str(aout)
        for suffix in (i.name.format(prefix="") for i in cls()):
            if aout_.endswith(suffix):
                return cls.from_prefix(aout_[: -len(suffix)])
        raise KeyError("Cannot determine suffix, please check file name")

    def _modify(self, *modify: PathLike | None):
        """
        Return a new instance with specified file paths replaced by new values.
        
        Parameters:
            modify: Optional new path values for each field; if a value is None, the original path is retained.
        
        Returns:
            A new instance of the same class with updated file paths.
        """
        return self.__class__(
            *(
                Path(modify) if modify else default
                for default, modify in zip(self, modify)
            )
        )

    TSV_COLS: tuple[str, ...] = ("All",)

    def _load_rep2all(self, *files):
        """
        Load and merge cluster membership TSV files into a single DataFrame.
        
        Each input file is read as a two-column TSV with column names determined by the class's `TSV_COLS` attribute. Files are merged sequentially on their shared columns to produce a unified mapping of representative to member sequences.
        
        Parameters:
            files: One or more file paths to TSV files containing cluster membership data.
        
        Returns:
            pd.DataFrame: Merged DataFrame mapping representative sequences to all members across clustering levels.
        """
        rep2all: pd.DataFrame
        for i, f in enumerate(files):
            names = self.TSV_COLS[i : i + 2][::-1]
            df = pd.read_csv(f, sep="\t", header=None, names=names)
            rep2all = df if i == 0 else rep2all.merge(df)
        return rep2all

    def load_rep2all(self, keep: Literal[True] = True) -> pd.DataFrame: """
Load and return the merged cluster membership DataFrame for this clustering result.

Parameters:
	keep (Literal[True] or tuple of str): If True, returns all columns; if a tuple of column names, returns only those columns.

Returns:
	pd.DataFrame: DataFrame containing cluster membership information, filtered by the specified columns.
"""
...

    @property
    def rep2all(self):
        """
        Load and return the full cluster membership DataFrame for all available columns.
        """
        return self.load_rep2all(keep=True)

    if TYPE_CHECKING:

        def __iter__(self):
            """
            Yields a single empty Path object when the instance is iterated over.
            """
            yield Path()

    @classmethod
    def exec_rep2all(
        cls,
        files: Iterable[PathLike] | None = None,
        faas: Iterable[SeqRecord.SeqRecord] | None = None,
        keep_prefix: PathLike | Literal[False] | None = None,
        threads=4,
        profile: PathLike | None = None,
    ):
        """
        Executes the gene clustering Snakemake workflow on provided sequence files or records and returns the resulting cluster membership data as a DataFrame.
        
        Parameters:
            files (Iterable[PathLike] | None): Input FASTA file paths containing sequences to cluster.
            faas (Iterable[SeqRecord.SeqRecord] | None): Sequence records to cluster, written to a temporary FASTA file if provided.
            keep_prefix (PathLike | Literal[False] | None): If set, moves the temporary input FASTA to this prefix directory and uses it as the output location; if False or None, uses a temporary file.
            threads (int): Number of threads to use for Snakemake execution.
            profile (PathLike | None): Optional Snakemake profile directory for workflow execution.
        
        Returns:
            pandas.DataFrame: DataFrame containing cluster membership information as defined by the class's `rep2all` property.
        
        Raises:
            RuntimeError: If the Snakemake workflow fails to execute successfully.
            NotImplementedError: If execution path is incomplete or not implemented.
        """
        assert files or faas
        with NamedTemporaryFile("w", suffix=".faa", delete=True) as tmpf:
            with open(tmpf.name) as fi:
                if faas:
                    SeqIO.write(faas, tmpf.name, "fasta-2line")
                    keep_prefix = keep_prefix or False
                if files:
                    for i, file in enumerate(files):
                        with open(file) as fi:
                            while True:  # read 16 Kib words one time
                                if not (block := fi.read(65536)):
                                    break
                                tmpf.write(block)
                            tmpf.flush()
                    if keep_prefix is None and i == 0:
                        # infer filename automatically if not given in some cases
                        prefix = Path(file).with_suffix("")
                        if str(file) == cls.in_faa(prefix):
                            keep_prefix = False
                            tmpf.name = str(file)
            if keep_prefix:
                Path(keep_prefix).mkdir(parents=True, exist_ok=True)
                shutil.move(tmpf.name, cls.in_faa(keep_prefix))
                tpmf_out = cls.from_prefix(keep_prefix)
            else:
                tpmf_out = cls.from_prefix(tmpf.name[:-4])
            tpmf_out_str = " ".join([str(i) for i in tpmf_out])
            if profile:
                profile_str = f"--profile {profile}"
            smk_params = (
                f"-s {GENE_CLUST_SMK} "
                f"{tpmf_out_str} "
                f"--use-conda "
                f"{profile_str} "
                f"-c{threads} -p "
            )

            try:
                os.system(f"ls {tmpf.name}")
                print("params:", "snakemake", smk_params)
                smk(smk_params)
            except SystemExit as se:
                if se.code:
                    print(se.code, se.with_traceback(None))
                    raise RuntimeError("snakemake seems not run successfully.")
                else:
                    return tpmf_out.rep2all

        raise NotImplementedError("")


class _UniRefClu(NamedTuple):
    u100: Path = Path("{prefix}-uniref100.tsv")
    u90: Path = Path("{prefix}-uniref90.tsv")
    u50: Path = Path("{prefix}-uniref50.tsv")


class UniRefClu(_UniRefClu, _CluBase):
    @classmethod
    def from_prefix_modify(
        cls,
        prefix: PathLike,
        u100: PathLike | None = None,
        u90: PathLike | None = None,
        u50: PathLike | None = None,
    ):
        """
        Create a new instance from a prefix, optionally overriding the default output file paths.
        
        Parameters:
            prefix: The base path used to construct default output file paths.
            u100: Optional custom path for the UniRef100 output file.
            u90: Optional custom path for the UniRef90 output file.
            u50: Optional custom path for the UniRef50 output file.
        
        Returns:
            An instance with file paths set according to the prefix and any provided overrides.
        """
        return cls.from_prefix(prefix)._modify(u100, u90, u50)

    TSV_COLS = "All", "U100", "U90", "U50"

    def load_rep2all(
        self,
        keep: Iterable[Literal["All", "U100", "U90", "U50"]] | Literal[True] = (
            "All",
            "U50",
        ),
    ):
        """
        Load and return UniRef cluster membership data as a DataFrame, selecting specified columns.
        
        Parameters:
            keep: If True, returns all columns; otherwise, returns only the columns listed.
        
        Returns:
            pandas.DataFrame containing cluster membership information for the selected columns.
        """
        keep_ = list(self.TSV_COLS if keep is True else keep)
        return self._load_rep2all(self.u100, self.u90, self.u50)[keep_]


def extract(
    subset: Collection[str],
    files: Iterable[PathLike] | None = None,
    faas: Iterable[SeqRecord.SeqRecord] | None = None,
):
    """
    Yield sequence records whose IDs are present in the given subset.
    
    Either sequence files or an iterable of sequence records must be provided. For each sequence, yields only those whose ID is found in `subset`.
    
    Parameters:
        subset (Collection[str]): Set of sequence IDs to extract.
        files (Iterable[PathLike], optional): Sequence file paths to search for matching records.
        faas (Iterable[SeqRecord.SeqRecord], optional): Iterable of sequence records to filter.
    
    Yields:
        SeqRecord.SeqRecord: Sequence records with IDs in `subset`.
    """
    assert files or faas
    if faas:
        for faa in faas:
            if faa.id in subset:
                yield faa
    if files:
        for file in files:
            for faa in SeqIO.parse(file, "fasta"):
                if faa.id in subset:
                    yield faa


class _MmseqOut(NamedTuple):
    all_100: Path = Path("{prefix}-clu_100.tsv")
    all_clu: Path = Path("{prefix}-clu.tsv")
    all_clu_faa: Path = Path("{prefix}-clu_rep.faa")


class MmseqOut(_MmseqOut, _CluBase):
    @classmethod
    def from_prefix_modify(
        cls,
        prefix: PathLike,
        all_100: PathLike | None = None,
        all_clu: PathLike | None = None,
        all_clu_faa: PathLike | None = None,
    ):
        """
        Create a new instance using a prefix, with optional overrides for output file paths.
        
        Parameters:
            prefix: The base path used to construct default output file paths.
            all_100: Optional path to override the default 100% identity cluster TSV file.
            all_clu: Optional path to override the default cluster TSV file.
            all_clu_faa: Optional path to override the default representative sequences FASTA file.
        
        Returns:
            An instance with file paths set according to the prefix and any provided overrides.
        """
        return cls.from_prefix(prefix)._modify(all_100, all_clu, all_clu_faa)

    TSV_COLS = "All", "Rep100", "Rep"

    def load_rep2all(
        self,
        keep: Iterable[Literal["All", "Rep100", "Rep"]] | Literal[True] = (
            "All",
            "Rep",
        ),
    ):
        """
        Load and return cluster membership data for MMseqs2 clustering as a DataFrame.
        
        Parameters:
            keep: Iterable of column names to retain in the output DataFrame, or True to keep all columns. Valid options are "All", "Rep100", and "Rep".
        
        Returns:
            pandas.DataFrame: Cluster membership data with the specified columns.
        """
        keep_ = list(self.TSV_COLS if keep is True else keep)
        return self._load_rep2all(self.all_100, self.all_clu)[keep_]


def mmseq_clust(
    files: Iterable[PathLike] | None = None,
    faas: Iterable[SeqRecord.SeqRecord] | None = None,
    out_prefix: PathLike | MmseqOut = "gene",
    threads=4,
) -> MmseqOut:
    # infer gff_out automatically if not given in some cases
    """
    Run MMseqs2 clustering on input sequence files or records and return output file paths.
    
    Parameters:
        files: Iterable of input FASTA file paths, or None if using sequence records.
        faas: Iterable of SeqRecord objects, or None if using input files.
        out_prefix: Output prefix path or MmseqOut instance specifying output file locations.
        threads: Number of threads to use for clustering.
    
    Returns:
        MmseqOut: An instance containing paths to the MMseqs2 clustering output files.
    """
    if not isinstance(out_prefix, tuple):
        _out_prefix = MmseqOut.from_prefix(out_prefix)
    else:
        _out_prefix = MmseqOut(*out_prefix)
    with NamedTemporaryFile("w", suffix=".faa", delete=True) as tmpf:
        MmseqOut.exec_rep2all(files, faas, tmpf.name[:-4], threads)
        tpmf_out = MmseqOut.from_prefix(tmpf.name[:-4])
        os.system(f"ls {tmpf.name}")
        for ffrom, fto in zip(tpmf_out, _out_prefix):
            Path(fto).parent.mkdir(parents=True, exist_ok=True)
            shutil.move(ffrom, fto)
        return _out_prefix


class _MmFamily(NamedTuple):
    mf100: Path = Path("{prefix}-mf100.tsv")
    mfamily: Path = Path("{prefix}-mfamily.tsv")


class MmFamily(_MmFamily, _CluBase):
    @classmethod
    def from_prefix_modify(
        cls,
        prefix: PathLike,
        mf100: PathLike | None = None,
        mfamily: PathLike | None = None,
    ):
        """
        Create an instance using a prefix, with optional overrides for the `mf100` and `mfamily` file paths.
        
        Parameters:
        	mf100: Optional path to override the default 100% family cluster TSV file.
        	mfamily: Optional path to override the default family cluster TSV file.
        
        Returns:
        	Instance with file paths set according to the prefix and any provided overrides.
        """
        return cls.from_prefix(prefix)._modify(mf100, mfamily)

    TSV_COLS = "All", "F100", "Family"

    def load_rep2all(
        self,
        keep: Iterable[Literal["All", "F100", "Family"]] | Literal[True] = (
            "All",
            "Family",
        ),
    ):
        """
        Load and return the merged family clustering membership DataFrame, optionally filtering columns.
        
        Parameters:
            keep: If True, returns all columns; otherwise, returns only the specified columns ("All", "F100", "Family").
        
        Returns:
            pandas.DataFrame: DataFrame containing cluster membership information for the selected columns.
        """
        keep_ = list(self.TSV_COLS if keep is True else keep)
        return self._load_rep2all(self.mf100, self.mfamily)[keep_]


class _MmSpecies(NamedTuple):
    """https://www.nature.com/articles/s41559-024-02357-0

    All gene sequences were clustered using MMseqs2 (ref. 42) with minimum overlap of 50%,
        minimum identity threshold of 80% and clustering mode 0.
        The rest of the parameters were left as default.
    """

    mf100: Path = Path("{prefix}-mf100.tsv")
    mspecies: Path = Path("{prefix}-mspecies.tsv")


class MmSpecies(_MmSpecies, _CluBase):
    @classmethod
    def from_prefix_modify(
        cls,
        prefix: PathLike,
        mf100: PathLike | None = None,
        mspecies: PathLike | None = None,
    ):
        """
        Create an instance from a prefix, optionally overriding the `mf100` and `mspecies` file paths.
        
        Parameters:
            prefix: The base path used to construct default output file paths.
            mf100: Optional custom path for the 100% family cluster TSV file.
            mspecies: Optional custom path for the species cluster TSV file.
        
        Returns:
            An instance with file paths set according to the prefix and any provided overrides.
        """
        return cls.from_prefix(prefix)._modify(mf100, mspecies)

    TSV_COLS = "All", "F100", "Species"

    def load_rep2all(
        self,
        keep: Iterable[Literal["All", "F100", "Species"]] | Literal[True] = (
            "All",
            "Species",
        ),
    ):
        """
        Load and return the merged DataFrame of species-level cluster memberships.
        
        Parameters:
            keep: If True, include all columns; otherwise, include only the specified columns ("All", "F100", "Species").
        
        Returns:
            pandas.DataFrame: DataFrame containing cluster membership information for the selected columns.
        """
        keep_ = list(self.TSV_COLS if keep is True else keep)
        return self._load_rep2all(self.mf100, self.mspecies)[keep_]
