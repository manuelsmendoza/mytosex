import numpy as np
import os
import pandas as pd
import subprocess as sp
import sys


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
        elif exclude is None:
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
        elif exclude is None:
            exclude = 1796

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
                  + sample_prefix + ".fixmate.bam " \
                  + sample_prefix + ".bam"
    for cmd in [filter_mapped, group_names, fixmate, sort_position, deduplicate]:
        print(cmd)
        out = sp.run(
            cmd,
            shell=True,
            capture_output=True
        )
        print(out)
