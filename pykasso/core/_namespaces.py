"""
Module defining some constants in pyKasso.
"""

### Internal dependencies
from PIL import Image


MISC_DIR_PATH = '/../_misc/'
DEFAULT_PARAMETERS_FILENAME = 'parameters.yaml'
DEFAULT_PROJECT_FILENAME = 'project.yaml'
DEFAULT_LOG_FILENAME = 'project.log'

GRID_PARAMETERS = [
    'x0',
    'y0',
    'z0',
    'nx',
    'ny',
    'nz',
    'dx',
    'dy',
    'dz'
]

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

VALID_EXTENSIONS_DATAFRAME = [
    'gslib',
    'vox',
    # 'grd',  # TODO
]

VALID_EXTENSIONS_IMAGE = [key.strip('.') for key in Image.EXTENSION.keys()]

VALID_EXTENSIONS_DATA = [
    'gslib',
    'vox',
    'csv',
    'txt',
    'npy',
    'tif',  # TODO - rasterio
    'tiff',  # TODO - rasterio
    'asc',  # TODO - rasterio
]
VALID_EXTENSIONS_DATA.extend(VALID_EXTENSIONS_IMAGE)

# Define default fast-marching costs
DEFAULT_FMM_COSTS = {
    # 'out': 0.999,
    'out': 10,  # TODO
    'geology': 0.4,
    # 'aquifer': 0.4,
    # 'aquiclude': 0.8,
    'beddings': 0.35,  # TODO
    'faults': 0.2,
    'fractures': 0.2,
    'karst': 0.1,
    'conduits': 0.1,
    'ratio': 0.5,
}

DEFAULT_FEATURE_PARAMETERS = {
    'geology': {
        'nodata': 1,
        'name': 'unit {}',
        'model': True,
    },
    'faults': {
        'nodata': 0,
        'name': 'fault {}',
        'model': True,
    },
    'fractures': {
        'nodata': 0,
        'name': 'family {}',
        'model': True,
    },
}