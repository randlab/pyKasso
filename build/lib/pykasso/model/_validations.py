"""
Input validation functions.
"""

import PIL
import sys
import logging
import rasterio
import numpy as np

from .._utils import datareader
from .._utils import validation as val

this = sys.modules[__name__]


#################
### FUNCTIONS ###
#################

def read_file(path: str, attribute: str) -> np.ndarray:
    extension = path.split('.')[-1]
    try:
        ### GSLIB
        if extension == 'gslib':
            data = np.genfromtxt(path, skip_header=3, dtype=np.int8)

        ### Numpy_pickle
        elif extension == 'npy':
            data = np.load(path)

        ### Images
        elif extension in ['jpg', 'png']:
            data = np.asarray(PIL.Image.open(path).convert('L')).T

        ### CSV
        elif extension == 'csv':
            data = np.genfromtxt(path, delimiter=',').T
            
        ### TIF, TIFF
        elif extension in ['tif', 'tiff']:
            data = rasterio.open(path).read(1).T

        ### Others
        else:
            data = np.genfromtxt(path)

    except Exception as err:
        msg = ("Impossible to read the file designated by the '{}' attribute."
               " Location : {}".format(attribute, path))
        this.logger.error(msg)
        raise err
    else:
        return data


def is_list_length_valid(data: list, value: int, attribute: str) -> bool:
    if len(data) < value:
        msg = ("'{}' data length is too short ({} elements minimum)."
               .format(attribute, value))
        this.logger.critical(msg)
        raise ValueError(msg)
    else:
        return True

     
def is_coordinate_type_valid(coordinate: tuple, types: tuple,
                             attribute: str) -> bool:
    if not isinstance(coordinate, types):
        msg = ("The values of the '{}' attribute contains at least one invalid"
               " vertex. Coordinates must be of type : {}."
               .format(attribute, types))
        this.logger.critical(msg)
        raise TypeError(msg)
    else:
        return True
                   

def is_surface_dimensions_valid(attribute: str, array: np.ndarray,
                                grid) -> bool:
    nx, ny, nz = grid.shape
    if not (array.shape == (nx, ny)):
        msg = ("The '{}' array shape does not match with grid surface."
               " Array shape: {}, Grid surface shape: {}"
               .format(attribute, array.shape, (nx, ny)))
        this.logger.critical(msg)
        raise ValueError(msg)
    else:
        return True
    

def is_costs_dictionnary_valid(costs_dictionnary: dict, ids_data: list):
    """ """
    for i in ids_data:
        if i not in costs_dictionnary:
            msg = ("The data id ({}) is not within 'costs' dictionnary keys"
                   " ({})".format(i, list(costs_dictionnary.keys())))
            this.logger.error(msg)
            raise KeyError(msg)
    return True
    
    
###############################################################################

########################
### OUTLETS - INLETS ###
########################

def validate_settings_points(settings: dict,
                             attribute: str,
                             ) -> dict:
    """
    Validate the parameters of ``outlets`` and ``inlets`` keys.
    
    TODO
    
    Tested parameters:
     - ``seed``: must be of type int;
     - ``number``: must be of type int, value must be greater than zero;
     - ``data``: TODO
    """
    # Set logger
    logger = logging.getLogger("{}.validation".format(attribute))
    
    ### 'seed' ###
    
    # Test if 'seed' is of type int or None
    try:
        if settings['seed'] is not None:
            val.is_variable_type_valid(variable_name='seed',
                                       variable_value=settings['seed'],
                                       valid_types=(int))
    except TypeError as error:
        logger.error(error)
        raise
    
    ### 'subdomain' ###
    
    # Test if 'subdomain' is of type str
    try:
        val.is_variable_type_valid(variable_name='subdomain',
                                   variable_value=settings['subdomain'],
                                   valid_types=(str))
    except TypeError as error:
        logger.error(error)
        raise
    
    ### 'shuffle' ###
    
    # Test if 'shuffle' is of type bool
    try:
        val.is_variable_type_valid(variable_name='shuffle',
                                   variable_value=settings['shuffle'],
                                   valid_types=(bool))
    except TypeError as error:
        logger.error(error)
        raise
    
    ### 'number' ###
    
    # Test if 'number' is of type int
    try:
        val.is_variable_type_valid(variable_name='number',
                                   variable_value=settings['number'],
                                   valid_types=(int))
    except TypeError as error:
        logger.error(error)
        raise

    # Test if 'number' is greater than zero
    try:
        val.is_parameter_comparison_valid(parameter_name='number',
                                          parameter_value=settings['number'],
                                          logical_test='>',
                                          compared_to=0)
    except ValueError as error:
        logger.error(error)
        raise
    
    ### 'importance' ###
    
    # Test if 'importance' is of type list
    try:
        val.is_variable_type_valid(variable_name='importance',
                                   variable_value=settings['importance'],
                                   valid_types=(list))
    except TypeError as error:
        logger.error(error)
        raise
    
    # Test if length of the list is adequate with declared number of points
    try:
        if len(settings['importance']) > settings['number']:
            # TODO - write error msg
            msg = ""
            raise Exception(msg)
    # TODO - custom exception ?
    except Exception as error:
        logger.error(error)
        raise
    
    ### 'data'
    
    # Test if data is empty
    if isinstance(settings['data'], str) and (settings['data'] == ''):
        settings['data'] = []
    
    # If 'data' type is str, try to read the file
    if isinstance(settings['data'], str):
        filepath = settings['data']

        # Test if filepath is valid
        val.is_filepath_valid(filepath)

        # Try to read filepath
        dr = datareader.DataReader()
        settings['data'] = dr.get_data_from_file(filepath)
    
    # If 'data' type is np.ndarray, transform it to list
    if isinstance(settings['data'], np.ndarray):
        
        # If the list of points contains only one element
        if len(settings['data'].shape) == 1:
            settings['data'] = np.array([settings['data']])
        
        settings['data'] = settings['data'].tolist()
        
    return settings
