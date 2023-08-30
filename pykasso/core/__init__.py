"""
Contains the core of pyKasso.

Please note that this module is private.  All functions and objects
are available in the main ``pykasso`` namespace - use that instead.

"""

__all__ = []

# Import pykasso function
from pykasso.core.main import pykasso
__all__.extend(['pykasso'])

# Import pykasso namespaces
# from pykasso.core._namespaces import (GRID_PARAMETERS,
#                                       GEOLOGICAL_FEATURES,
#                                       SURFACE_FEATURES,
#                                       DOMAIN_FEATURES,
#                                       ISOTROPIC_FEATURES,
#                                       ANISOTROPIC_FEATURES,
#                                       AUTHORIZED_FEATURES)
# __all__.extend(['GRID_PARAMETERS', 'GEOLOGICAL_FEATURES', 'SURFACE_FEATURES',
#                 'DOMAIN_FEATURES', 'ISOTROPIC_FEATURES',
#                 'ANISOTROPIC_FEATURES', 'AUTHORIZED_FEATURES'])
