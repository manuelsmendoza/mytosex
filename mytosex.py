""" Description:
MyToSex: In silico sex determination of bivalve with mitochondrial double uniparental inheritance based on the
mitochondrial genomes content.
"""

import argparse as arg
import os
import shutil as sh
import sys
from src.func_setup import *
from src.func_util import *

__author__ = "Manuel Mendoza"
__version__ = "0.1"
__license__ = "MIT"


def main():
    # command-line arguments
    parser = arg.ArgumentParser(description="Sex inference from mitotypes content using RNA-Seq data")
    parser.add_argument("settings", metavar="settings.yaml", type=str, nargs=1, help="settings file")
    parser.add_argument("-v", "--version", action="version", version=__version__, help="show the current version")
    args = parser.parse_args()
    os.environ["MYTOSEX_SETTINGS"] = os.path.abspath(args.settings[0])
    os.environ["CUDA_VISIBLE_DEVICES"] = ""
    os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"

    # Common values
    output_dir = os.path.join(load_settings(os.path.abspath(args.settings[0]))["output_dir"], "mytosex_result")
    settings_checked = os.path.join(output_dir, "settings.json")

    print(tnow() + " INFO: Starting mitochondrial analysis", file=sys.stdout)

    # Check for the dependencies and settings values
    if not check_dir(output_dir):
        import run.setup
    elif check_dir(output_dir) and not check_file(os.path.join(output_dir, ".settings.ok")):
        import run.setup
    else:
        os.environ["MYTOSEX_SETTINGS"] = settings_checked
        print(tnow() + " WARN: Settings checked previously", file=sys.stdout)

    # Download the reads and sequences (if needed)
    settings = load_settings(settings_checked)
    if any(list(settings["from_ncbi"].values())) and not check_file(os.path.join(output_dir, ".download.ok")):
        import run.fetch
    elif any(list(settings["from_ncbi"].values())) and check_file(os.path.join(output_dir, ".download.ok")):
        print(tnow() + " WARN: Samples reads and mitogenomes downloaded previously", file=sys.stdout)

    # Perform the analysis
    if not check_file(os.path.join(output_dir, ".analysis.ok")):
        import run.analysis
    else:
        print(tnow() + " WARN: Sex determination done", file=sys.stdout)

    # Predict the sex of the samples
    if not check_file(os.path.join(output_dir, ".prediction.ok")):
        import run.prediction
    else:
        print(tnow() + " WARN: Sex inference done", file=sys.stdout)

    # Run the phylogenetic analysis
    if "other_spp" in list(settings.keys()) and not check_file(os.path.join(output_dir, ".phylo.ok")):
        import run.phylo
    else:
        print(tnow() + " WARN: Skipping phylogenetic analysis", file=sys.stdout)

    sh.rmtree(os.path.join(settings["output_dir"], "tmp"))
    print(tnow() + " INFO: Analysis finished", file=sys.stdout)
    sys.exit()


if __name__ == "__main__":
    main()
