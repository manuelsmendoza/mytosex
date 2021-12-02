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
