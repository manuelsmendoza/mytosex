""" Description:
MyToSex: In silico sex determination of bivalve with mitochondrial double uniparental inheritance based on the
mitochondrial genomes content.
"""

import argparse as arg

__author__ = "Manuel Mendoza"
__version__ = "0.1"
__license__ = "MIT"


def main():
    """
    Insert here more information
    """
    # command-line arguments
    parser = arg.ArgumentParser(description="Sex inference from mitotypes content using RNA-Seq data")
    parser.add_argument("settings", metavar="settings.yaml", type=str, nargs=1, help="settings file")
    parser.add_argument("-v", "--version", action="version", version=__version__, help="show the current version")

    args = parser.parse_args()

    # Check if the dependencies are installed
    from analysis.check_dep import check_depmod


if __name__ == "__main__":
    main()
