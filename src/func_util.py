"""
Some useful functions:
    - tnow: Show the current time
"""

import datetime as dt


def tnow(format="%d %b %Y %H:%M:%S %Z"):
    """
    tnow(format)

    Capture the current time and export it as human-readable string.

    Parameters
    ----------
    format: str, optional
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
    ctime = "[" + dt.datetime.now().strftime(format) + "]"
    ctime = ctime.replace(" ]", "]")
    return ctime
