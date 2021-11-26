"""
Checking the parameters provided in the settings file
"""

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

    res = module in sys.modules.keys()
    return res


def check_settings(settings_file):
    """
    check_settings(settings.yaml)

    Check all the values provided in the settings file.

    Parameters
    ----------
    settings_file: str
        Absolut path to the settings file

    See Also
    --------
    More information about yaml format https://yaml.org
    """
    with open(settings_file, "r") as infile:
        try:
            settings = yaml.load(
                stream=infile,
                Loader=yaml.FullLoader
            )
        except ImportError:
            print("Settings file could not be loaded")


resmod = {}
modules_required = [
    "yaml",
    "numpy"
]
for mod in modules_required:
    resmod.update({mod: check_depmod(mod)})
print(resmod)
