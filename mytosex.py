""" Description:
MyToSex: In silico sex determination of bivalve with mitochondrial double uniparental inheritance based on the
mitochondrial genomes content.
"""

import argparse as arg
import os
import yaml

__author__ = "Manuel Mendoza"
__version__ = "0.1"
__license__ = "MIT"

# Configuration of command-line arguments
parser = arg.ArgumentParser(
    description="Sex inference from mitotypes content using RNA-Seq data"
)
parser.add_argument(
    "settings",
    metavar="settings.yaml",
    type=str,
    nargs=1,
    help="settings file"
)
parser.add_argument(
    "-v", "--version",
    action="version",
    version=__version__,
    help="show the current version"
)
args = parser.parse_args()

# Load settings values
with open(args.settings[0], "r") as settings_file:
    try:
        settings = yaml.load(
            stream=settings_file,
            Loader=yaml.FullLoader
        )
    except ImportError:
        print("Settings file could not be loaded")

print(settings)
