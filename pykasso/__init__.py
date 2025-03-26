"""
pyKasso
=======

pyKasso is a python3 open-source package intended to simulate easily and
quickly karst networks using a geological model, hydrogeological, and
structural data.

License
-------
Released under the GPL-3.0 license.
Copyright (C) 2025 University of Neuchâtel - CHYN.
 - François Miville <francois.miville@ikmail.com>
 - Philippe Renard  <philippe.renard@unine.ch>
 - Chloé Fandel     <cfandel@email.arizona.edu>

Available subpackages
---------------------
core
   Karstic conduit network generator
analysis
   Karstic conduit network analysis tool
visualization
   Karstic conduit network visualization tool
   
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

# Import pyKasso version string
from ._version import __version__
__all__.extend(['__version__'])
