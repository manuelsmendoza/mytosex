"""
Checking the parameters provided in the settings file
"""

import site
import sys
import os
from misc.util import tnow


def check_depmod(module):
    """
    check_depmod(module)

    Check if a python module is installed

    Parameters
    ----------
    module: str
        The name of python module to check

    Returns
    -------
    bool
        Result of searching if the module is installed
    """

    path = site.getsitepackages()
    res = []
    for p in path:
        res.append(os.path.exists(os.path.join(p, module)))

    if any(res):
        return True
    else:
        return False


resmod = {}
modules_required = [
    "yaml",
    "numpy"
]
for mod in modules_required:
    resmod.update({mod: check_depmod(mod)})
print(resmod)
