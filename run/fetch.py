import os
import pandas as pd
import sys
from Bio import SeqIO
from src.func_fetch import *
from src.func_setup import load_settings
from src.func_util import tnow

# Load settings
settings = load_settings(os.getenv("MYTOSEX_SETTINGS"))
tmp_dir = os.path.join(settings["output_dir"], "tmp")

# Download the sequences and export them together with the annotation and features coordinates
if settings["from_ncbi"]["seqs"]:
    for mt_type in ["mtf", "mtm"]:
        if settings["reference"][mt_type + "_ncbi"]:
            print(tnow() + " INFO: Downloading the " + mt_type + " of reference", file=sys.stdout)
            alias = settings["reference"]["alias"] + "_" + mt_type
            record = annotate(fetch_sequence(settings["reference"][mt_type]), alias)
            export_record(record, tmp_dir, alias)

    for mt_type in ["mt", "mtf", "mtm"]:
        for specie in list(settings["other_spp"].keys()):
            if mt_type in list(settings["other_spp"][specie].keys()):
                if settings["other_spp"][specie][mt_type + "_ncbi"]:
                    specie_name = settings["other_spp"][specie]["alias"]
                    print(tnow() + " INFO: Downloading the " + mt_type + " of " + specie_name, file=sys.stdout)
                    alias = specie_name + "_" + mt_type
                    record = annotate(fetch_sequence(settings["other_spp"][specie][mt_type]), alias)
                    export_record(record, tmp_dir, alias)


# Download the samples reads
if settings["from_ncbi"]["reads"]:
    for sample in list(settings["samples"].keys()):
        if settings["samples"][sample]["ncbi"]:
            print(tnow() + " INFO: Downloading the reads of " + settings["samples"][sample]["alias"], file=sys.stdout)
            fetch_reads(
                settings["samples"][sample]["accession"],
                tmp_dir,
                settings["samples"][sample]["alias"],
                settings["samples"][sample]["layout"],
                settings["numb_threads"]
            )

open(os.path.join(settings["output_dir"], ".download.ok"), "w").close()
