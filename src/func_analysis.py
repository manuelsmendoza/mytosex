import numpy as np
import os
import pandas as pd
import subprocess as sp


def build_index(sequence, index, threads):
    """ Build bowtie2 index
    build_index(/path/sequence.fasta, /path/index_basename, threads)

    Parameters
    ----------
    sequence : str
        Path to the sequence
    index : str
        Path to the index
    threads : int
        Number of threads to use
    """
    cmd = "bowtie2-build " \
          + "--threads " + str(threads) + " " \
          + sequence + " " \
          + index
    out = sp.run(
        cmd,
        shell=True,
        capture_output=True
    )


def reads_align(index, alias, outdir, threads, freads=None, rreads=None, sreads=None):
    """ Align reads to a reference using bowtie2

    Parameters
    ----------
    index : str
        Path to the sequence index
    alias : str
        A human-fiendly name to call the reference
    outdir : str
        Path to the directory to write the alignments
    threads : int
        Number of threads to use
    freads : str
        Path to the forward reads if the layout is paired-end
    rreads : str
        Path to the reverse reads if the layout is paired-end
    sreads : str
        Path to the reads if the layout is single-end
    """
    if freads is not None and rreads is not None and sreads is None:
        cmd = "bowtie2 " \
              + "--threads " + str(threads) + " " \
              + "--end-to-end " \
              + "--very-sensitive " \
              + "-x " + index + " " \
              + "-1 " + freads + " " \
              + "-2 " + rreads + " " \
              + "-S " + os.path.join(outdir, alias + ".sam")
    elif freads is None and rreads is None and sreads is not None:
        cmd = "bowtie2 " \
              + "--threads " + str(threads) + " " \
              + "--end-to-end " \
              + "--very-sensitive " \
              + "-x " + index + " " \
              + "-U " + sreads + " " \
              + "-S " + os.path.join(outdir, alias + ".sam")
    out = sp.run(
        cmd,
        shell=True,
        capture_output=True
    )


def filter_alignment(alignment, threads, layout, require=None, exclude=None):
    """ Filter reads alignment and mark duplications

    Parameters
    ----------
    alignment : str
        Path to the raw alignment
    threads : int
        Number of threads to use for filtering the alignments
    layout : str {single, paired}
        Sequencing layout
    require : int
        Flags required for the alignments
    exclude : int
        Flags to exclude the alignments
    """
    if not os.path.exists(alignment):
        raise FileNotFoundError("Alignment file not found")
    elif os.path.exists(alignment) and not os.access(alignment, os.R_OK):
        raise PermissionError("Permission denied to read the alignment")

    sample_prefix = os.path.splitext(alignment)[0]
    if layout == "paired":
        fixmate = "samtools fixmate " \
                  + "-r " \
                  + "-m " \
                  + "-O bam " \
                  + "--threads " + str(threads) + " " \
                  + sample_prefix + ".collate.bam " \
                  + sample_prefix + ".fixmate.bam"
        if require is None:
            require = 3
        if exclude is None:
            exclude = 1804
    elif layout == "single" and require is None:
        fixmate = "samtools fixmate " \
                  + "-r " \
                  + "-m " \
                  + "-p " \
                  + "-O bam " \
                  + "--threads " + str(threads) + " " \
                  + sample_prefix + ".collate.bam " \
                  + sample_prefix + ".fixmate.bam"
        if require is None:
            require = 1
        if exclude is None:
            exclude = 1796
    compress_alignment = "samtools view " \
                         + "--output-fmt BAM " \
                         + "-@ " + str(threads) + " " \
                         + "-o  " + sample_prefix + ".bam " \
                         + alignment
    filter_mapped = "samtools view " \
                    + "--output-fmt BAM " \
                    + "-q 30 " \
                    + "-F " + str(exclude) + " " \
                    + "-f " + str(require) + " " \
                    + "-@ " + str(threads) + " " \
                    + "-o  " + sample_prefix + ".filtered.bam " \
                    + alignment
    group_names = "samtools collate " \
                  + "--output-fmt BAM " \
                  + "-@ " + str(threads) + " " \
                  + "-o " + sample_prefix + ".collate.bam " \
                  + sample_prefix + ".filtered.bam"
    sort_position = "samtools sort " \
                    + "-O bam " \
                    + "-@ " + str(threads) + " " \
                    + "-o " + sample_prefix + ".sort.bam " \
                    + sample_prefix + ".fixmate.bam"
    deduplicate = "samtools markdup " \
                  + "-r " \
                  + "-S " \
                  + "--threads " + str(threads) + " " \
                  + sample_prefix + ".sort.bam " \
                  + sample_prefix + ".markdup.bam"
    for cmd in [compress_alignment, filter_mapped, group_names, fixmate, sort_position, deduplicate]:
        out = sp.run(
            cmd,
            shell=True,
            capture_output=True
        )
    os.remove(alignment)


def extract_stats(alignment, features, threads):
    """ Extract alignment statistics

    Parameters
    ----------
    alignment : str
        Path to the reads alignment after filtering
    features : str
        Path to genome features coordinates (bed file)
    threads: int
        Number of threads to use
    """
    prefix = os.path.splitext(os.path.splitext(alignment)[0])[0]
    indx = "samtools index " \
           + "-b " \
           + "-@ " + str(threads) + " " \
           + alignment
    cov = "samtools coverage " \
          + "-q 30 " \
          + "-o " + prefix + ".cov.tsv" + " " \
          + alignment
    bcov = "samtools bedcov " \
           + "-Q 30 " \
           + "-d 0 " \
           + features + " " \
           + alignment + " > " \
           + prefix + ".bedcov.tsv"
    deep = "samtools depth " \
           + "-aa " \
           + "-o " + prefix + ".depth.tsv" + " " \
           + alignment
    for cmd in [indx, cov, bcov, deep]:
        out = sp.run(
            cmd,
            shell=True,
            capture_output=True
        )


def gini_index(variable):
    """ Calculate the Gini's index

    Parameters
    ----------
    variable : numpy.ndarray
        A numeric variable indexed in non-decreasing order

    Returns
    -------
    gi : float
        Gini's index value
    """
    if np.all(variable == 0):
        return np.NaN
    else:
        gindex = 2 * np.multiply(np.array(range(1, variable.size + 1)), np.sort(variable)).sum()
        gindex = gindex / (variable.size * variable.sum())
        gindex = gindex - (variable.size + 1) / variable.size

        return gindex


def alignment_stats(coverage, feat_coverage, depth, salias, ralias):
    """ Calculate alignment statistics

    Parameters
    ----------
    coverage : str
        Path to the genome coverage statistics
    feat_coverage : str
        Path to genomic features coverage
    depth : str
        Path to the sequencing depth base-per-pase
    salias : str
        The alias of the sample
    ralias: str
        The alias of the reference

    Returns
    -------
    align_stats : pandas.core.frame.DataFrame
        Alignment statistics
    """

    allcov = pd.read_csv(coverage, sep="\t")
    col_names = ["seqname", "start", "end", "name", "score", "strand", "counts", "length"]
    fetcov = pd.read_csv(feat_coverage, sep="\t", names=col_names)
    fetcov["depth"] = 100 * fetcov.loc[:, "counts"] / fetcov.loc[:, "length"]
    alldep = pd.read_csv(depth, sep="\t", names=["seqname", "position", "depth"])
    align_stats = {
        "sample": [salias],
        "mtfcov": [list(allcov.loc[allcov["#rname"] == ralias + "_mtf", "coverage"])[0]],
        "mtmcov": [list(allcov.loc[allcov["#rname"] == ralias + "_mtm", "coverage"])[0]],
        "mtfmd": [list(allcov.loc[allcov["#rname"] == ralias + "_mtf", "meandepth"])[0]],
        "mtmmd": [list(allcov.loc[allcov["#rname"] == ralias + "_mtm", "meandepth"])[0]],
        "mtfgi": [gini_index(np.array(list(alldep.loc[alldep["seqname"] == ralias + "_mtf", "depth"])))],
        "mtmgi": [gini_index(np.array(list(alldep.loc[alldep["seqname"] == ralias + "_mtm", "depth"])))],
    }
    align_stats = pd.DataFrame.from_dict(align_stats)

    return align_stats


def sex_round(x):
    """ Round the sex prediction index

    Parameters
    ----------
    x : int

    Returns
    ----------
    Round x
    """
    if x is np.nan:
        return np.nan
    else:
        return int(x.round())
