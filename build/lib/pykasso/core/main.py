"""
Module containing functions for accessing the public content of pyKasso.
"""

### Local dependencies
from .application import Application
from .grid import Grid
from .._utils.datareader import DataReader


def pykasso() -> Application:
    """
    Create and return an ``Application``.
    
    Returns
    -------
    Application
        
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


def create_datareader(grid: Grid = None) -> DataReader:
    """
    Create and return a ``DataReader``.

    Returns
    -------
    DataReader
    
    Parameters
    ----------
    Grid
        
    See Also
    --------
    DataReader
        
    Examples
    --------
    >>> import pykasso as pk
    >>> data_reader = pk.create_datareader()
    """
    out = DataReader(grid=grid)
    return out
