import os

def get_settings(example=False):
    """
    Provide the datafiles settings.

    Parameters
    ----------
    example : bool, optionnal
        If True, pyKasso will provide you the Betteraz's files example

    Examples
    --------
        >>> pk.get_settings(example=True)
    """

    # copying defaults file from source package
    import shutil
    import glob

    path = os.path.dirname(os.path.abspath(__file__)) + '/' + 'default_files' + '/'

    if example == True:

        # create inputs directory
        dst = 'inputs'

        if not os.path.exists(dst):
            try:
                os.makedirs(dst)
                print("Directory named '", dst, "' has been created.")
            except FileExistsError:
                raise

        path = path + 'betteraz' + '/'

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
    Private function to simply open a file.
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
    Private function to load points in a text. Generally combined with 'opendatafile()'.
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
