import os
import pandas as pd
import subprocess as sp
from Bio import SeqIO
from src.func_util import check_file


def extract_cds(sequence, annotation, outdir, alias):
    """ Extract the CDS sequences from a genome and its annotation

    Parameters
    ----------
    sequence : str
        Path to the genome sequences
    annotation : str
        Path to the genome annotation
    outdir : str
        Path to directory to write the CDS sequences
    alias : str
        Human friendly sequence name
    """
    cmd = "gffread " \
          + "-g " + sequence + " " \
          + "-x " + os.path.join(outdir, alias + "_cds.fasta") + " " \
          + annotation
    out = sp.run(
        cmd,
        shell=True,
        capture_output=True
    )


def create_db(sequences, outdir):
    """ Create a local BLAST database

    Parameters
    ----------
    sequences : str
        Path to the sequences to build the database
    outdir : str
        Path to the directory to wirte te output
    alias : str
        Sequences common name
    """
    cmd = "makeblastdb " \
          + "-dbtype nucl " \
          + "-in " + sequences + " " \
          + "-out " + os.path.join(outdir, os.path.splitext(os.path.basename(sequences))[0])
    out = sp.run(
        cmd,
        shell=True,
        capture_output=True
    )


def extract_reads(alignment, output_dir, layout, alias, threads=1):
    """ Extract the reads mapped propperly to the genome

    Parameters
    ----------
    alignment : str
        Path to the alignment after filtering process
    output_dir : str
        Path to the directory to write the output into
    layout : str {single, paired}
        Sequencing layout
    alias : str
        Sample common name
    threads : int
        Number of threads to use
    """
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

    if layout == "single":
        cmd = "samtools fasta " \
              + "-@ " + str(threads) + " " \
              + "-o " + os.path.join(output_dir,  alias + ".fasta") + " " \
              + alignment
        reads = [os.path.join(output_dir,  alias + ".fasta")]
    elif layout == "paired":
        cmd = "samtools fasta " \
              + "-@ " + str(threads) + " " \
              + "-1 " + os.path.join(output_dir,  alias + "_1.fasta") + " " \
              + "-2 " + os.path.join(output_dir,  alias + "_2.fasta") + " " \
              + alignment
        reads = [os.path.join(output_dir,  alias + "_1.fasta"),
                 os.path.join(output_dir,  alias + "_2.fasta")]

    out = sp.run(
        cmd,
        shell=True,
        capture_output=True
    )
    for reads_file in reads:
        compress_file(reads_file, threads=threads)


def transcripts_assembly(alignment, outdir, threads, maxmem, layout, alias, freads=None, rreads=None, sreads=None):
    """ Assembly and identify the mitogenes sequences

    Parameters
    ----------
    freads : str
        Path to the forward reads (if paired-end)
    rreads : str
        Path to the reverse reads (if paired-end)
    sreads : str
        Path to the reads (is single-end)
    alignment : str
        Path to the alignment after filtering process
    refseq : str
        Sequence of the mitogenome of reference
    refann : str
        Annotation of the mitogenome of reference
    outdir : str
        Path to the output directory
    alias : str
        Human-friendly name to handle the sample
    threads : int
        Number of threads to use
    maxmem : int
        Maximum amount of memory to use
    layout : str {single, paired}
        Sequencing layout
    """
    assembly_dir = os.path.join(outdir, alias + "_trinity")
    if not os.path.exists(assembly_dir):
        os.mkdir(assembly_dir)

    assembly_cmd = "Trinity " \
                   + "--bypass_java_version_check " \
                   + "--seqType fa " \
                   + "--max_memory " + str(maxmem) + "G " \
                   + "--CPU " + str(threads) + " " \
                   + "--genome_guided_bam " + alignment + " " \
                   + "--genome_guided_max_intron 300 " \
                   + "--output " + assembly_dir
    if layout == "single":
        assembly_cmd += " --single " + sreads
    elif layout == "paired":
        assembly_cmd += " --left " + freads + " --right " + rreads

    identify_orf_cmd = "TransDecoder.LongOrfs " \
                       + "-G Mitochondrial-Invertebrates " \
                       + "-t " + os.path.join(assembly_dir, "Trinity-GG.fasta") + " " \
                       + "-O " + assembly_dir
    predict_orf_cmd = "TransDecoder.Predict " \
                      + "-G Mitochondrial-Invertebrates " \
                      + "-t " + os.path.join(assembly_dir, "Trinity-GG.fasta") + " " \
                      + "-O " + assembly_dir

    os.chdir(assembly_dir)
    for cmd in [assembly_cmd, identify_orf_cmd, predict_orf_cmd]:
        out = sp.run(
            cmd,
            shell=True,
            capture_output=True
        )


def annotate_cds(codseq, database, alias, outdir, threads):
    """ Identify the gene encoded by a sequence by homology

    Parameters
    ----------
    codseq : str
        Path to the file containing all the protein-coding sequences
    database : str
        Nucleotide database to use
    alias : str
        Human-friendly name to identify the sample
    outdir : str
        Path to the directo to write the output
    threads : int
        Number of threads to use
    """
    cmd = "blastn " \
          + "-query " + codseq + " " \
          + "-db " + database + " " \
          + "-out " + os.path.join(outdir, alias + ".outfmt6") + " " \
          + "-evalue 1e-16 " \
          + "-outfmt 6 " \
          + "-max_hsps 1 " \
          + "-max_target_seqs 1 " \
          + "-num_threads " + str(threads)

    out = sp.run(
        cmd,
        shell=True,
        capture_output=True
    )


def build_annotation(codseq, codann, ref_alias, sample_alias):
    """ Rebuild the coding sequences annotation to convert into gff format

    Parameters
    ----------
    codseq : str
        Path to the codging sequences
    codann : str
        Path to the annotation results (blastn output)
    ref_alias : str
        Human-friendly name to call the reference
    sample_alias : str
        Sample human-friendly name
    """
    colnames = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send",
                "evalue", "bitscore"]
    all_annotation = pd.read_csv(codann, sep="\t", names=colnames)
    seq_annotation = []
    for seq in SeqIO.parse(codseq, "fasta"):
        cds_info = seq.description.split(" ")
        if cds_info[0] in list(all_annotation.loc[:, "qseqid"]):
            seq_att = list(all_annotation.loc[all_annotation["qseqid"] == cds_info[0], "sseqid"])[0]
            seq_annotation.append(
                {
                    "seqname": cds_info[-1].split(":")[0],
                    "source": "BLAST",
                    "feature": "CDS",
                    "start": cds_info[-1].split(":")[-1].split("-")[0],
                    "end": cds_info[-1].split(":")[-1].split("-")[1].split("(")[0],
                    "score": ".",
                    "strand": cds_info[-1].split(":")[-1].split("(")[-1].replace(")", ""),
                    "frame": int(0),
                    "attribute": "ID=" + seq_att.replace(ref_alias, sample_alias)
                }
            )

    return pd.DataFrame(seq_annotation)


def multiple_alignment(seq_in, alg_out, threads=1, iterations=1000):
    """ Align multiple sequences

    Parameters
    ----------
    seq_in : str
        File containing the sequences to align
    alg_out : str
        Path to output file with the alignments
    threads : int
        Number of threads to use
    iterations : int
        Maximum number of iterative refinement
    """
    cmd = "mafft-ginsi " \
          + "--maxiterate " + str(iterations) + " " \
          + "--thread " + str(threads) + " " \
          + seq_in + " > " + alg_out

    out = sp.run(
        cmd,
        shell=True,
        capture_output=True
    )
