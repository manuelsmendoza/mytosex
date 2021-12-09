import os
import sys
from Bio import SeqIO
from src.func_phylo import *
from src.func_setup import load_settings
from src.func_util import tnow

# Load settings
settings = load_settings(os.getenv("MYTOSEX_SETTINGS"))
tmp_dir = os.path.join(settings["output_dir"], "tmp")
data_dir = os.path.join(settings["output_dir"], "data")

if not os.path.exists(os.path.join(data_dir, "msa")):
    os.mkdir(os.path.join(data_dir, "msa"))

# Create a local database to annotate the genes
# print(tnow() + " INFO: Creating a local database to annotate the genes", file=sys.stdout)
# extract_cds(
#     sequence=os.path.join(tmp_dir, settings["reference"]["alias"] + ".fasta"),
#     annotation=os.path.join(tmp_dir, settings["reference"]["alias"] + ".gff"),
#     outdir=tmp_dir,
#     alias=settings["reference"]["alias"]
# )
# create_db(
#     sequences=os.path.join(tmp_dir, settings["reference"]["alias"] + "_cds.fasta"),
#     outdir=tmp_dir
# )

# Assemble the transcriptomes and extract the CDS
for sample in list(settings["samples"].keys()):
    sname = settings["samples"][sample]["alias"]
    # print(tnow() + " INFO: Extracting mitochondrial reads of " + sname, file=sys.stdout)
    # extract_reads(
    #     alignment=os.path.join(tmp_dir, settings["samples"][sample]["alias"] + ".fixmate.bam"),
    #     alias=settings["samples"][sample]["alias"],
    #     output_dir=tmp_dir,
    #     layout=settings["samples"][sample]["layout"],
    #     threads=settings["numb_threads"]
    # )

    # print(tnow() + " INFO: Assembling the mitogenome of " + sname, file=sys.stdout)
    # if settings["samples"][sample]["layout"] == "single":
    #     transcripts_assembly(
    #         layout="single",
    #         sreads=os.path.join(tmp_dir, settings["samples"][sample]["alias"] + ".fasta.gz"),
    #         alignment=os.path.join(tmp_dir, settings["samples"][sample]["alias"] + ".markdup.bam"),
    #         outdir=tmp_dir,
    #         threads=settings["numb_threads"],
    #         maxmem=settings["max_memory"],
    #         alias=settings["samples"][sample]["alias"]
    #     )
    # elif settings["samples"][sample]["layout"] == "paired":
    #     transcripts_assembly(
    #         layout="paired",
    #         freads=os.path.join(tmp_dir, settings["samples"][sample]["alias"] + "_1.fasta.gz"),
    #         rreads=os.path.join(tmp_dir, settings["samples"][sample]["alias"] + "_2.fasta.gz"),
    #         alignment=os.path.join(tmp_dir, settings["samples"][sample]["alias"] + ".markdup.bam"),
    #         outdir=tmp_dir,
    #         threads=settings["numb_threads"],
    #         maxmem=settings["max_memory"],
    #         alias=settings["samples"][sample]["alias"]
    #     )
    
    # print(tnow() + " INFO: Annotating the genes of " + settings["samples"][sample]["alias"], file=sys.stdout)
    # wd = os.path.join(tmp_dir, settings["samples"][sample]["alias"] + "_trinity")
    # annotate_cds(
    #     codseq=os.path.join(wd, "Trinity-GG.fasta.transdecoder.cds"),
    #     database=os.path.join(tmp_dir, settings["reference"]["alias"] + "_cds"),
    #     outdir=tmp_dir,
    #     alias=settings["samples"][sample]["alias"],
    #     threads=settings["numb_threads"]
    # )
    # 
    # ann = build_annotation(
    #     codseq=os.path.join(wd, "Trinity-GG.fasta.transdecoder.cds"),
    #     codann=os.path.join(tmp_dir, settings["samples"][sample]["alias"] + ".outfmt6"),
    #     ref_alias=settings["reference"]["alias"],
    #     sample_alias=settings["samples"][sample]["alias"]
    # )
    # ann.to_csv(
    #     os.path.join(tmp_dir, settings["samples"][sample]["alias"] + ".gff"),
    #     sep="\t",
    #     index=False,
    #     header=False
    # )
    
    # print(tnow() + " INFO: Extracting coding sequences from " + settings["samples"][sample]["alias"], file=sys.stdout)
    # extract_cds(
    #     sequence=os.path.join(tmp_dir, settings["samples"][sample]["alias"] + "_trinity", "Trinity-GG.fasta"),
    #     annotation=os.path.join(tmp_dir, settings["samples"][sample]["alias"] + ".gff"),
    #     outdir=tmp_dir,
    #     alias=settings["samples"][sample]["alias"]
    # )


# Merge the mitogenomes of other species in a single file
# print(tnow() + " INFO: Merging the mitogenomes of other species in a single file", file=sys.stdout)
# for specie in list(settings["other_spp"].keys()):
#     for mt in ["mt", "mtf", "mtm"]:
#         if mt in list(settings["other_spp"][specie].keys()):
#             prefix = settings["other_spp"][specie]["alias"] + "_" + mt
#             extract_cds(
#                 sequence=os.path.join(tmp_dir, prefix + ".fasta"),
#                 annotation=os.path.join(tmp_dir, prefix + ".gff"),
#                 outdir=tmp_dir,
#                 alias=prefix
#             )
#
# for mt in ["mtf", "mtm"]:
#     prefix = settings["reference"]["alias"] + "_" + mt
#     extract_cds(
#         sequence=os.path.join(tmp_dir, prefix + ".fasta"),
#         annotation=os.path.join(tmp_dir, prefix + ".gff"),
#         outdir=tmp_dir,
#         alias=prefix
#     )

# Store each gene into a single file
cds_files = list(filter(lambda f: "_cds.fasta" in f, os.listdir(tmp_dir)))
cds_files.remove(settings["reference"]["alias"] + "_cds.fasta")
cds_files = [os.path.join(tmp_dir, x) for x in cds_files]
for gene in ["ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND5", "ND6", "ND4L"]:
    print(tnow() + " INFO: Compiling the sequences all the sequences realated to " + gene, file=sys.stdout)
    with open(os.path.join(tmp_dir, gene + ".fasta"), "w") as gene_file:
        for file in cds_files:
            for rec in SeqIO.parse(file, "fasta"):
                if rec.id.endswith(gene):
                    SeqIO.write(rec, gene_file, "fasta")
    gene_file.close()

    print(tnow() + " INFO: Aligning the multiple sequences of " + gene, file=sys.stdout)
    multiple_alignment(
        seq_in=os.path.join(tmp_dir, gene + ".fasta"),
        alg_out=os.path.join(data_dir, "msa", gene + ".fasta"),
        threads=settings["numb_threads"]
    )
