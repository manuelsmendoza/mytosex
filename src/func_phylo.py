import os
import subprocess as sp
from src.func_util import check_file


def compress_file(file, threads=1):
    """ Compress a file

    Parameters
    ----------
    file : str
        Path to the file to compress
    threads : int
        Number of threads to use
    """
    if check_file(file):
        cmd = "pigz " \
              + "--force " \
              + "--best " \
              + "--processes " + str(threads) + " " \
              + file
    else:
        raise FileExistsError("File not found")
    
    out = sp.run(
        cmd,
        shell=True,
        capture_output=True
    )


def extract_reads(alignment, output_dir, layout, threads=1):
    """ Extract the reads mapped propperly to the genome

    Parameters
    ----------
    alignment : str
        Path to the alignment after filtering process
    output_dir : str
        Path to the directory to write the output into
    layout : str {single, paired}
        Sequencing layout
    threads : int
        Number of threads to use
    """
    sample_name = os.path.basename(os.path.splitext(os.path.splitext(alignment)[0])[0])
    if layout == "single":
        cmd = "samtools fasta " \
              + "-@ " + str(threads) + " " \
              + "-o " + os.path.join(output_dir,  sample_name + ".fasta") + " " \
              + alignment
        reads = [os.path.join(output_dir,  sample_name + ".fasta")]
    elif layout == "paired":
        cmd = "samtools fasta " \
              + "-@ " + str(threads) + " " \
              + "-1 " + os.path.join(output_dir,  sample_name + "_1.fasta") + " " \
              + "-2 " + os.path.join(output_dir,  sample_name + "_2.fasta") + " " \
              + alignment
        reads = [os.path.join(output_dir,  sample_name + "_1.fasta"), 
                 os.path.join(output_dir,  sample_name + "_2.fasta")]

    out = sp.run(
        cmd,
        shell=True,
        capture_output=True
    )
    for reads_file in reads:
        compress_file(reads_file, threads=threads)
    