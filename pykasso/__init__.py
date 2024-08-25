"""
pyKasso
=======

pyKasso is a python3 open-source package intended to simulate easily and
quickly karst networks using a geological model, hydrogeological, and
structural data.

Documentation
-------------
Complete API documentation with a variety of examples are available here :
TODO

License
-------
Released under the GPL-3.0 license.
Copyright (C) 2023 University of Neuchâtel - CHYN.
 - François Miville <francois@miville.org>
 - Chloé Fandel     <cfandel@email.arizona.edu>
 - Philippe Renard  <philippe.renard@unine.ch>

Available subpackages
---------------------
core
   Kartic conduit network generator
analysis
   Karstic conduit network analysis tool
visualization
   Karst conduit network visualization tool
   
Utilities
---------
__version__
    pyKasso version string
"""

__all__ = []

# Import pyKasso's core
from . import core
from .core import *
__all__.extend(core.__all__)

# Import DataReader class
# from pykasso._utils.datareader import DataReader
# __all__.extend(['DataReader'])

# Import pyKasso version string
from ._version import __version__
__all__.extend(['__version__'])
