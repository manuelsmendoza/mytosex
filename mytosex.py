""" Description:
MyToSex: In silico sex determination of bivalve with mitochondrial double uniparental inheritance based on the
mitochondrial genomes content.
"""

import argparse as arg
import os
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

    # Check for the dependencies and settings values
    output_dir = os.path.join(load_settings(os.path.abspath(args.settings[0]))["output_dir"], "mytosex_result")
    if not check_dir(output_dir) and not check_file(os.path.join(output_dir, ".settings.ok")):
        import run.setup
    elif check_dir(output_dir) and not check_file(os.path.join(output_dir, ".settings.ok")):
        import run.setup
    else:
        print(tnow() + " INFO: Settings checked previously", file=sys.stdout)


if __name__ == "__main__":
    main()
