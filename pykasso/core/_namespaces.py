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
