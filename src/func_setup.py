import sys
import os
from src.func_util import tnow


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


def check_settings(set_file=os.getenv("MYTOSEX_SETTINGS")):
    """ Check the settings provided to perform the analysis

    """
