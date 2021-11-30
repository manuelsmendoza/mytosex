import sys
import os
import yaml
from Bio import Entrez
from xml.etree import ElementTree


def check_dep(dep_name, dep_type):
    """ Check if a dependency is installed
    check_dep(dep_name, dep_type)

    Parameters
    ----------
    dep_name: str
        The name of the dependency to check
    dep_type: str {"module", "tool"}
        The type of dependency, tool or module

    Returns
    -------
    bool
        The result of searching for the dependency in the system
    """
    if dep_type == "module":
        path = sys.path
    elif dep_type == "tool":
        path = os.getenv("PATH").split(":")
    else:
        raise ValueError("Dependencies only can be modules or tools")

    detection = []
    for directory in path:
        detection.append(os.path.exists(os.path.join(directory, dep_name)))
    if any(detection):
        return True
    else:
        return False


def load_settings(settings_file):
    """ Load the settings values from a yaml file
    laod_settings(settings.yaml)

    Parameters
    ----------
    settings_file: str
        Path to the settings value

    Returns
    -------
    dict
        The settings to run the analysis
    """
    def check_file(path):
        """ Check if a regular file exists and the USER has permission to read it
        check_file(path)

        Parameters
        ----------
        path: str
            The path to the file to check

        Returns
        -------
        bool
            True if the file exists and it is readable by the user
        """
        if not os.path.exists(path):
            raise FileExistsError("File not found")
        elif os.path.exists(path) and not os.access(path, os.R_OK):
            raise PermissionError("Permission denied to read the file")
        else:
            return True

    if check_file(settings_file):
        with open(settings_file) as settings_file:
            settings = yaml.load(stream=settings_file, Loader=yaml.FullLoader)

    return settings


def check_settings(settings_file):
    """ Verify if the settings file exists and is readable by the user
    check_settings(settings.yaml)

    Parameters
    ----------
    settings_file: str
        Path to the settings file in YAML format

    Returns
    -------
    dict
        A settings checked file. That contains more information after run the multiples checking steps:
            - Computational resources to use if there were not specified previously
            - Output directory: add mytosex_result to the original output directory
            - Samples from NCBI: bool
            - Sequencing layout: {paired_end, single_end}
    """
    def check_file(path):
        """ Check if a regular file exists and the USER has permission to read it
        check_file(path)

        Parameters
        ----------
        path: str
            The path to the file to check

        Returns
        -------
        bool
            True if the file exists and it is readable by the user
        """
        if not os.path.exists(path):
            raise FileExistsError("File not found")
        elif os.path.exists(path) and not os.access(path, os.R_OK):
            raise PermissionError("Permission denied to read the file")
        else:
            return True

    def check_dir(path):
        """ Check if a directory exists and the USER has permission to write inside
        check_dir(path)

        Parameters
        ----------
        path: str
            The path to the directory to check

        Returns
        -------
        bool
            True if the directory exists and is writeable
        """
        if not os.path.exists(path):
            raise FileExistsError("Directory not found")
        elif not os.path.isdir(path):
            raise NotADirectoryError("Directory does not exists")
        elif not os.access(path, os.W_OK):
            raise PermissionError("Permission denied to write inside the directory")
        else:
            return True

    def check_reference(sequence):
        """ Check if the sequence and its annotation is available and readable
        check_reference(sequence.fasta)

        Parameters
        ----------
        sequence: str
            Path to the mitogenome sequence

        Returns
        -------
        bool
            True if the sequence and its annotation is avaible
        """
        if check_file(sequence) and any([check_file(os.path.splitext(sequence)[0] + ext) for ext in [".gff", ".gtf"]]):
            return True
        else:
            raise FileNotFoundError("Reference annotation missed")

    def check_accession(accession, ncbi_database):
        """ Verify if an accession number is correct
        check_accession(acc, db)

        Parameters
        ----------
        accession: str
            The sequence or sample accession ID
        ncbi_database: str {"Nucleotide", "SRA"}
            Database where the accession belongs to

        Returns
        -------
        bool
            True if the accession number is correct
        """

        Entrez.email = "mytisexdev@gmail.com"
        ncbi_search = Entrez.esearch(db=ncbi_database, term=accession + "[ACCN]", retmax=100, idtype="acc")
        if int(Entrez.read(ncbi_search)["Count"]) > 0:
            return True
        else:
            raise ValueError("Wrong accession number")

    def check_layout(accession):
        """ Check the sequencing layout without downloading the NCBI SRA sample
        check_layout(acc)

        Parameters
        ----------
        accession: str
            Sample accession number

        Returns
        -------
        str
            The sequencing layout, paired-end ot single-end
        """
        ncbi_result = Entrez.efetch(db="sra", id=accession, rettype="full", retmode="xml")
        result_root = ElementTree.parse(ncbi_result).getroot()
        result_tags = [element.tag for element in result_root.iter()]
        seq_layout = result_tags[result_tags.index("LIBRARY_LAYOUT") + 1]
        return seq_layout.lower()

    # Check if the file exists and can be read and load them
    if check_file(settings_file):
        with open(settings_file) as settings_file:
            settings = yaml.load(stream=settings_file, Loader=yaml.FullLoader)

    # Check the output directory
    if "output_dir" in list(settings.keys()) and check_dir(settings["output_dir"]):
        settings["output_dir"] = os.path.join(settings["output_dir"], "mytosex_result")
    else:
        raise ValueError("Output directory missed")

    # Checking the resources to use
    if "numb_threads" not in list(settings.keys()) or settings["numb_threads"] > os.cpu_count():
        settings["numb_threads"] = os.cpu_count()

    tot_mem = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES') / 1024 ** 3
    if "max_memory" not in list(settings.keys()) or settings["max_memory"] > tot_mem:
        settings["max_memory"] = tot_mem

    # Checking the mitogenomes of reference
    if "reference" not in list(settings.keys()):
        raise ValueError("Reference information not provided")
    else:
        if "alias" not in list(settings["reference"].keys()):
            raise ValueError("Reference alias missed")
        for mt_type in ["mtf", "mtm"]:
            if mt_type in list(settings["reference"].keys()):
                if settings["reference"][mt_type].count("/") > 0 and check_reference(settings["reference"][mt_type]):
                    settings["reference"].update({mt_type + "_ncbi": False})
                elif check_accession(settings["reference"][mt_type], "Nucleotide"):
                    settings["reference"].update({mt_type + "_ncbi": True})
            else:
                raise ValueError("Some mitotype of reference missed")

    # Checking samples
    if "samples" not in list(settings.keys()):
        raise ValueError("Samples information nor provided")
    else:
        for sample in list(settings["samples"].keys()):
            sample_keys = list(settings["samples"][sample].keys())
            if "alias" not in sample_keys:
                raise ValueError("Samples alias missed")
            elif "accession" in sample_keys:
                accession = settings["samples"][sample]["accession"]
                if check_accession(accession, "SRA"):
                    settings["samples"][sample].update({"ncbi": True})
                    settings["samples"][sample].update({"layout": check_layout(accession)})
            elif "single" in sample_keys:
                reads = settings["samples"][sample]["single"]
                if check_file(reads):
                    settings["samples"][sample].update({"ncbi": False})
                    settings["samples"][sample].update({"layout": "single"})
            elif "forward" in sample_keys and "reverse" in sample_keys:
                reads = [settings["samples"][sample]["forward"], settings["samples"][sample]["reverse"]]
                if check_file(reads[0]) and check_file(reads[1]):
                    settings["samples"][sample].update({"ncbi": False})
                    settings["samples"][sample].update({"layout": "paired"})

    # Check for other species information
    if "other_spp" in list(settings.keys()):
        for specie in list(settings["other_spp"].keys()):
            if "alias" not in list(settings["other_spp"][specie].keys()):
                raise ValueError("Additional specie alias missed")
            else:
                for mt_type in ["mt", "mtf", "mtm"]:
                    if mt_type in list(settings["other_spp"][specie].keys()):
                        reference = settings["other_spp"][specie][mt_type]
                        if reference.count("/") > 0 and check_reference(reference):
                            settings["other_spp"][specie].update({mt_type + "_ncbi": False})
                        elif check_accession(reference, "Nucleotide"):
                            settings["other_spp"][specie].update({mt_type + "_ncbi": True})

    # Add new item to download the sequences or reads
    ncbi_seqs = []
    for mt_type in ["mt", "mtf", "mtm"]:
        if mt_type in list(settings["reference"].keys()):
            ncbi_seqs.append(settings["reference"][mt_type + "_ncbi"])
        for specie in list(settings["other_spp"].keys()):
            if mt_type in list(settings["other_spp"][specie].keys()):
                ncbi_seqs.append(settings["other_spp"][specie][mt_type + "_ncbi"])

    ncbi_reads = []
    for sample in list(settings["samples"].keys()):
        ncbi_reads.append(settings["samples"][sample]["ncbi"])

    settings.update({"from_ncbi": {"seqs": any(ncbi_seqs), "reads": any(ncbi_reads)}})

    # Return the new settings values
    return settings
