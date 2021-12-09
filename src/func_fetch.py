import os.path

import pandas as pd
import subprocess as sp
from Bio import Entrez
from Bio import SeqIO
from Bio import SeqRecord


def fetch_reads(accession, outdir, alias, layout, threads=1):
    """ Download reads from the NCBI SRA database
    fetch_reads(accession, /path/reads_dir)

    Parameters
    ----------
    accession : str
        The sample accession number in the NCBI SRA database
    outdir : str
        Directory to store the reads
    alias : str
        A user-friendly name to handle the samples
    layout : str {single, paired}
        Sequencing layout
    threads : int
        Number of threads to download the sample in parallel
    """
    def compress_file(file, threads):
        """ Compress a file

        Parameters
        ----------
        file : str
            Path to the file to compress
        threads : int
            Number of threads to use for the compression
        """
        cmd = "pigz " \
              + "--best " \
              + "--force " \
              + "--processes " + str(threads) + " " \
              + file
        out = sp.run(
            cmd,
            shell=True,
            capture_output=True
        )

    cmd = "fasterq-dump " \
          + "--split-files " \
          + "--outdir " + outdir + " " \
          + "--temp " + outdir + " " \
          + "--threads " + str(threads) + " " \
          + accession
    out = sp.run(
        cmd,
        shell=True,
        capture_output=True
    )

    if layout == "paired":
        for strand in ["1", "2"]:
            os.rename(
                os.path.join(outdir, accession + "_" + strand + ".fastq"),
                os.path.join(outdir, alias + "_" + strand + ".fastq")
            )
            compress_file(os.path.join(outdir, alias + "_" + strand + ".fastq"), threads)
    else:
        os.rename(
            os.path.join(outdir, accession + ".fastq"),
            os.path.join(outdir, alias + ".fastq")
        )
        compress_file(os.path.join(outdir, alias + ".fastq"), threads)


def fetch_sequence(accession):
    """ Download a sequence from NCBI Nucleotide
    fetch_sequence(accession)

    Parameters
    ----------
    accession : str
        The sequence accession number

    Returns
    -------
    sequence : Bio.SeqRecord
        The sequence and its features
    """
    Entrez.email = "mytisexdev@gmail.com"
    sequence_gb = Entrez.efetch(
        db="Nucleotide",
        id=accession,
        rettype="gb",
        retmode="text"
    )
    sequence_record = SeqIO.read(sequence_gb, 'genbank')

    return sequence_record


def annotate(record, alias):
    """ Extract the annotation from a genbank record

    Parameters
    ----------
    record : Bio.SeqRecord
        A genbank record containing the sequence and its annotation
    alias : str
        A common name to use to identify the sequence

    Returns
    -------
    dict :
        A dictionary containing the sequence, its annotation in gff and the features coordinates
    """
    def find_frame(qualifiers):
        """ Find the starting frame
        find_frame(qualifiers)

        Parameters
        ----------
        qualifiers : dict
            Information about the features from a sequence record

        Returns
        -------
        int
            The starting frame
        """
        qualifiers = list(qualifiers.keys())
        for item in range(0, len(qualifiers) - 1):
            if qualifiers[item] == "codon_start":
                return item

    # Modify the description line removing redundant fields
    record = record
    record.id = alias
    record.name = ""
    record.desc = ""

    # Extract the sequence annotation and transform it into a DataFrame
    annotation = []
    for feature in record.features:
        feature_annotation = {
            "seqname": record.id,
            "sourcer": "GenBank",
            "feature": feature.type,
            "start": feature.location.nofuzzy_start,
            "end": feature.location.nofuzzy_end,
            "score": ".",
            "strand": feature.location.strand
        }
        if feature.type in ["rRNA", "tRNA"]:
            feature_annotation.update(
                {
                    "frame": ".",
                    "attribute": "ID=" + dict(feature.qualifiers)["product"][0]
                }
            )
        elif feature.type == "CDS":
            feature_annotation.update(
                {
                    "frame": find_frame(feature.qualifiers),
                    "attribute": "ID=" + dict(feature.qualifiers)["gene"][0]
                }
            )
        annotation.append(feature_annotation)
    annotation = pd.DataFrame(annotation)

    # Rename the genes to have all of them with the same nomenclature
    annotation.replace(
        {"strand": [1, -1]},
        {"strand": ["+", "-"]},
        regex=True,
        inplace=True
    )
    annotation.replace(
        {"attribute": [r"^ID=COI$", r"^ID=COII$", r"^ID=COIII$"]},
        {"attribute": [r"ID=COX1", r"ID=COX2", r"ID=COX3"]},
        regex=True,
        inplace=True
    )
    annotation.replace(
        {"attribute": [r"^ID=cob"]},
        {"attribute": [r"ID=CYTB"]},
        regex=True,
        inplace=True
    )
    annotation.replace(
        {"attribute": [r"^ID=cox1", r"^ID=cox2", r"^ID=cox3"]},
        {"attribute": [r"ID=COX1", r"ID=COX2", r"ID=COX3"]},
        regex=True,
        inplace=True
    )
    annotation.replace(
        {"attribute": [r"ID=NAD", r"ID=nad"]},
        {"attribute": [r"ID=ND", r"ID=ND"]},
        regex=True,
        inplace=True
    )
    annotation.replace(
        {"attribute": [r"ID=ND4l"]},
        {"attribute": [r"ID=ND4L"]},
        regex=True,
        inplace=True
    )
    annotation.replace(
        {"attribute": [r"ID=ATPase", r"ID=atp"]},
        {"attribute": [r"ID=ATP", r"ID=ATP"]},
        regex=True,
        inplace=True
    )
    annotation.replace(
        {"attribute": [r" ribosomal RNA$"]},
        {"attribute": [r""]},
        regex=True,
        inplace=True
    )
    annotation.replace(
        {"attribute": [r"16S$", r"12S$", r"small subunit", r"large subunit"]},
        {"attribute": [r"RNR2", r"RNR1", r"RNR1", r"RNR2"]},
        regex=True,
        inplace=True
    )
    annotation.replace(
        {"attribute": [r"tRNA-Tyr$", r"tRNA-Lys$", r"tRNA-Met$", r"tRNA-Leu$", r"tRNA-Val$", r"tRNA-Ser$"]},
        {"attribute": [r"TRNY", r"TRNK", r"TRNM", r"TRNL", r"TRNV", r"TRNS"]},
        regex=True,
        inplace=True
    )
    annotation.replace(
        {"attribute": [r"tRNA-Arg$", r"tRNA-Trp$", r"tRNA-Ala$", r"tRNA-His$", r"tRNA-Pro$", r"tRNA-Thr$"]},
        {"attribute": [r"TRNR", r"TRNW", r"TRNA", r"TRNH", r"TRNP", r"TRNT"]},
        regex=True,
        inplace=True
    )
    annotation.replace(
        {"attribute": [r"tRNA-Phe$", r"tRNA-Gly$", r"tRNA-Asn$", r"tRNA-Glu$", r"tRNA-Cys$", r"tRNA-Ile$"]},
        {"attribute": [r"TRNF", r"TRNG", r"TRNN", r"TRNE", r"TRNC", r"TRNI"]},
        regex=True,
        inplace=True
    )
    annotation.replace(
        {"attribute": [r"tRNA-Gln$", r"tRNA-Asp$"]},
        {"attribute": [r"TRNQ", r"TRND"]},
        regex=True,
        inplace=True
    )
    annotation.replace(
        {"attribute": [r"^ID="]},
        {"attribute": [r"ID=" + alias + "_"]},
        regex=True,
        inplace=True
    )
    filter_condition = [feature in ["CDS", "tRNA", "rRNA"] for feature in list(annotation["feature"])]
    annotation = annotation[filter_condition]

    coordinates_fields = ["seqname", "start", "end", "attribute", "score", "strand"]
    coordinates = annotation.loc[:, coordinates_fields]
    coordinates.replace(
        {"attribute": "^ID="},
        {"attribute": ""},
        regex=True,
        inplace=True
    )

    sequence = SeqRecord.SeqRecord(
        record.seq,
        id=alias,
        name="",
        description=""
    )

    result = {
        "sequence": sequence,
        "annotation": annotation,
        "coordinates": coordinates
    }
    return result


def export_record(sequence_info, outdir, alias):
    """ Export the sequence, its annotation and the features coordinates

    Parameters
    ----------
    sequence_info : dict
        Sequence information fetched from the NCBI
    outdir : str
        Directory to write the information
    alias : str
        A human-friendly name to use to handle the samples information
    """
    # Export the sequence
    with open(os.path.join(outdir, alias + ".fasta"), "w") as out_sequence:
        SeqIO.write(sequence_info["sequence"], out_sequence, "fasta")
    out_sequence.close()

    # Export the annotation and the features coordinates
    sequence_info["annotation"].to_csv(os.path.join(outdir, alias + ".gff"), sep="\t", header=False, index=False)
    sequence_info["coordinates"].to_csv(os.path.join(outdir, alias + ".bed"), sep="\t", header=False, index=False)
