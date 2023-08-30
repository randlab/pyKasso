"""
Module defining some constants in pyKasso.
"""

# GRID_PARAMETERS = [
#     'x0',
#     'y0',
#     'z0',
#     'nx',
#     'ny',
#     'nz',
#     'dx',
#     'dy',
#     'dz'
# ]

GEOLOGICAL_FEATURES = [
    'domain',
    'geology',
    'faults',
    'fractures',
]

SURFACE_FEATURES = [
    'topography',
    'bedrock',
    'water_level',
]

DOMAIN_FEATURES = [
    'delimitation',
    'topography',
    'bedrock',
    'water_level',
]

ISOTROPIC_FEATURES = [
    'cost',
    'time',
    'karst',
]

ANISOTROPIC_FEATURES = [
    'cost',
    'alpha',
    'beta',
    'time',
    'karst',
    'gradient',
]

features = [
    GEOLOGICAL_FEATURES,
    DOMAIN_FEATURES,
    ANISOTROPIC_FEATURES
]

AUTHORIZED_FEATURES = [f for list_ in features for f in list_]
