"""
Contains the core of pyKasso: application, project, grid, etc.

This subpackage contains the core of pyKasso. It constructs an
application class able to manage a project and to communicate between
the different other subpackages.

Please note that this module is private.  All functions and objects
are available in the main ``pykasso`` namespace - use that instead.

"""

__all__ = []

# Import pykasso function
from pykasso.core.main import pykasso
__all__.extend(['pykasso'])
