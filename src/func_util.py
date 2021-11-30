"""
Some useful functions:
    - tnow(): Show the current time
"""

import datetime as dt
import os


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
    if os.path.exists(path) and os.access(path, os.R_OK):
        return True
    else:
        return False


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
    if os.path.exists(path) and os.path.isdir(path) and os.access(path, os.W_OK):
        return True
    else:
        return False


def tnow(time_format="%d %b %Y %H:%M:%S %Z"):
    """
    tnow(format)

    Capture the current time and export it as human-readable string.

    Parameters
    ----------
    time_format: str, optional
        Format to print the time

    Returns
    -------
    str
        Human readable time

    See Also
    --------
    datetime.datetime.strftime:
        Convert object to a string according to a given format using C standard codes: https://strftime.org

    Examples
    --------
    tnow()
    """
    ctime = "[" + dt.datetime.now().strftime(time_format) + "]"
    ctime = ctime.replace(" ]", "]")
    return ctime
