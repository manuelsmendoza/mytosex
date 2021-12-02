import os
import sys
from src.func_analysis import *
from src.func_setup import load_settings
from src.func_util import tnow, pass_file

# Load settings
settings = load_settings(os.getenv("MYTOSEX_SETTINGS"))
tmp_dir = os.path.join(settings["output_dir"], "tmp")
data_dir = os.path.join(settings["output_dir"], "data")
figs_dir = os.path.join(settings["output_dir"], "figures")

# Building the references of
print(tnow() + " INFO: Merging both mitogenomes of reference", file=sys.stdout)
for ext in [".fasta", ".gff", ".bed"]:
    open(os.path.join(tmp_dir, settings["reference"]["alias"] + ext), "w").close()
    for mt_type in ["mtf", "mtm"]:
        if settings["reference"][mt_type + "_ncbi"]:
            pass_file(
                os.path.join(tmp_dir, settings["reference"]["alias"] + "_" + mt_type + ext),
                os.path.join(tmp_dir, settings["reference"]["alias"] + ext)
            )
        else:
            pass_file(
                os.path.join(settings["reference"][mt_type]),
                os.path.join(tmp_dir, settings["reference"]["alias"] + ext)
            )

# Build the sequence index
print(tnow() + " INFO: Building the reference index", file=sys.stdout)
build_index(
    os.path.join(tmp_dir, settings["reference"]["alias"] + ".fasta"),
    os.path.join(tmp_dir, settings["reference"]["alias"]),
    settings["numb_threads"]
)
