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

    :param index: Reference index
    :param alias: Sample alias
    :param outdir: Directory to store the alignment
    :param threads: Number of threads to use
    :param freads: Forward reads, is paired-end
    :param rreads: Reverse reads, if paired-end
    :param sreads: Reads, if single-end
    """
    if freads is not None and rreads is not None and sreads is None:
        cmd = "bowtie2 " \
              + "--threads " + str(threads) + " " \
              + "--end-to-end " \
              + "--very-sensitive " \
              + "--no-mixed " \
              + "--no-discordant " \
              + "--no-unal " \
              + "--al-conc-gz " + os.path.join(outdir, alias) + " " \
              + "-x " + index + " " \
              + "-1 " + freads + " " \
              + "-2 " + rreads + " " \
              + "-S " + os.path.join(outdir, alias + ".sam")
    elif freads is None and rreads is None and sreads is not None:
        cmd = "bowtie2 " \
              + "--threads " + str(threads) + " " \
              + "--end-to-end " \
              + "--very-sensitive " \
              + "--no-unal " \
              + "--al-gz " + os.path.join(outdir, alias) + " " \
              + "-x " + index + " " \
              + "-U " + sreads + " " \
              + "-S " + os.path.join(outdir, alias + ".sam")
    else:
        print(tnow() + " FAIL: Wrong layout")
        raise ValueError("Forward and reverse reads are not compatible with single-end sequence samples")
    out = sp.run(
        cmd,
        shell=True,
        capture_output=True
    )
