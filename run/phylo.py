import os
import sys
from src.func_phylo import *
from src.func_setup import load_settings
from src.func_util import tnow

# Load settings
settings = load_settings(os.getenv("MYTOSEX_SETTINGS"))
tmp_dir = os.path.join(settings["output_dir"], "tmp")

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
    
    print(tnow() + " INFO: Annotating the genes of " + settings["samples"][sample]["alias"])
    wd = os.path.join(tmp_dir, settings["samples"][sample]["alias"] + "_trinity")
    annotate_cds(
        codseq=os.path.join(wd, "Trinity-GG.fasta.transdecoder.cds"),
        database=os.path.join(tmp_dir, settings["reference"]["alias"] + "_cds"),
        outdir=tmp_dir,
        alias=settings["samples"][sample]["alias"],
        threads=settings["numb_threads"]
    )