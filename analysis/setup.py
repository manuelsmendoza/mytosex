"""
Checking the parameters provided in the settings file
"""

import sys
import os
from misc.util import tnow


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


# Checking for the dependencies
print(tnow() + " INFO: Checking for the dependencies")

mod_name = ["yaml"]
tool_name = ["bowtie2"]
dep_res = {}
for dep in mod_name:
    dep_res.update({dep: check_dep(dep_name=dep, dep_type="module")})
for dep in tool_name:
    dep_res.update({dep: check_dep(dep_name=dep, dep_type="tool")})

if not all(list(dep_res.values())):
    for dep in list(dep_res.keys()):
        if not dep_res[dep]:
            print(tnow() + " FAIL: " + dep + " is not installed")
    raise ModuleNotFoundError
