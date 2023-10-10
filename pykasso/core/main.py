"""
Module storing functions giving access to the public content of pyKasso.
"""

### Local dependencies
from pykasso.core.application import Application
from pykasso._utils.datareader import DataReader


def pykasso() -> Application:
    """
    Create and return an ``Application``.
    
    TODO
    
    Returns
    -------
    Application
        TODO
        
    See Also
    --------
    Application, Project, Grid
    
    Examples
    --------
    >>> import pykasso as pk
    >>> app = pk.pykasso()
    """
    out = Application()
    return out


def create_datareader() -> DataReader:
    """
    Create and return a ``DataReader``.
    
    TODO

    Returns
    -------
    DataReader
        TODO
        
    See Also
    --------
    DataReader
        
    Examples
    --------
    >>> import pykasso as pk
    >>> data_reader = pk.create_datareader()
    """
    out = DataReader()
    return out
