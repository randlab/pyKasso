"""
Common functions that can be wrote and grouped here.
"""

import os

def get_settings(example=None):
    """
    Provides the datafiles settings.

    Parameters
    ----------
    example : str, optional
        Examples available : 'betteraz'.

    Examples
    --------
    >>> pk.get_settings(example='betteraz')
    """

    # copying defaults file from source package
    import shutil
    import glob

    path = os.path.dirname(os.path.abspath(__file__)) + '/' + 'default_files' + '/'

    if example is not None:

        # create inputs directory
        dst = 'inputs'

        if not os.path.exists(dst):
            try:
                os.makedirs(dst)
                print("Directory named '", dst, "' has been created.")
            except FileExistsError:
                raise

        path = path + example + '/'

        srcs = glob.glob(path + '*')
        dst  = dst + '/'

        for src in srcs:
            shutil.copy(src, dst)

    else:
        src = path + 'settings.yaml'
        dst  = 'settings.yaml'

        shutil.copy(src, dst)

    return None

def opendatafile(file_location):
    """
    Reads a file.

    Parameters
    ----------

    Returns
    -------

    Examples
    --------
    """
    try:
        file = open(file_location,"r")
        text = file.readlines()
        file.close()
    except:
        print("Error : impossible to open datafile.")
        raise
    return text

def loadpoints(text):
    """
    Loads points in a text. Generally combined with 'opendatafile()'.

    Parameters
    ----------

    Returns
    -------

    Examples
    --------
    """
    points = []
    try:
        for data_line in text:
            x, y = data_line.strip().split()
            points.append((float(x), float(y)))
    except ValueError:
        print("Dimensions of data frame does not match with number of variables. Missing values or values in excess.")
        raise
    except:
        print("Unexpected error:", sys.exc_info()[0])
        raise
    return points
