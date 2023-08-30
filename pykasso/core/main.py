"""
Module loading the default available contents from pyKasso.
"""

### Local dependencies
from pykasso.core.application import Application

    
def pykasso() -> Application:
    """
    Create and return a pyKasso application.
    
    Examples
    --------
    >>> import pykasso as pk
    >>> app = pk.pykasso()
    """
    out = Application()
    return out
