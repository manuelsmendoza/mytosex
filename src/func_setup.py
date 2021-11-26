import sys
import os
import yaml


def check_dep(dep_name, dep_type):
    """ Check if a dependency is installed
    check_dep(dep_name, dep_type)

    Parameters
    ----------
    dep_name: str
        The name of the dependency to check
    dep_type: {"module", "tool"}
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

    res = []
    for p in path:
        res.append(os.path.exists(os.path.join(p, dep_name)))
    if any(res):
        return True
    else:
        return False


def check_settings(set_file):
    """ Check the settings provided to perform the analysis
    check_settings(settings.yaml)

    Parameters
    ----------
    set_file: str
        Path to the settings file in YAML format

    Returns
    -------
    dict
        A dictionary containing all the settings information
    """
    if os.path.exists(set_file) and os.path.isfile(set_file) and os.access(set_file, os.R_OK):
        with open(set_file) as set_values:
            settings = yaml.load(set_values, Loader=yaml.FullLoader)
    else:
        raise FileNotFoundError("Check setting file")
