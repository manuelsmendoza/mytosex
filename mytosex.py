""" Description:
MyToSex: In silico sex determination of bivalve with mitochondrial double uniparental inheritance based on the
mitochondrial genomes content.
"""


import argparse as arg
import os

__author__ = "Manuel Mendoza"
__version__ = "0.1"
__license__ = "MIT"


def main():
    # command-line arguments
    parser = arg.ArgumentParser(description="Sex inference from mitotypes content using RNA-Seq data")
    parser.add_argument("settings", metavar="settings.yaml", type=str, nargs=1, help="settings file")
    parser.add_argument("-v", "--version", action="version", version=__version__, help="show the current version")

    args = parser.parse_args()
    print(args)
    #os.environ["MYTOSEX_SETTINGS"] = args[0]

    # Check for the dependencies and settings values
    import chunks.setup
    print(os.getenv("MYTOSEX_SETTINGS"))


if __name__ == "__main__":
    main()
