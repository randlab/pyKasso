"""
Input validations functions
"""

import PIL
import os
import sys
import logging
import rasterio
import numpy as np
from shapely.geometry import Point

from pykasso._utils import datareader
from pykasso.core._namespaces import DEFAULT_FMM_COSTS

from pykasso._utils import validation as val




# to remove ?
this = sys.modules[__name__]  # useful ?


#################
### FUNCTIONS ###
#################

def read_file(path: str, attribute: str) -> np.ndarray:
    """
    TODO
    """
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

# def validate_settings(feature: str, feature_settings: dict, grid=None) -> dict:
#     """
#     TODO
#     """
    
#     elif feature in ['geology', 'faults', 'fractures']:
#         settings = validate_geologic_feature_settings(feature_settings)
#         settings = validate_fractures_settings(feature_settings)
#     elif feature == 'domain':
#         settings = validate_domain_settings(feature_settings)
#     elif feature in ['inlets', 'outlets']:
#     settings = validate_points_feature_settings(feature_settings)

#     return settings

###########
### SKS ###
###########

def validate_settings_sks(settings: dict) -> dict:
    """TODO
    
    - 'seed'
    - 'algorithm'
    - 'costs'
    - 'mode'
    - 'factors'
    """
    logger = logging.getLogger("sks.validation")
    
    ### 'seed'
    # Must be an int
    try:
        val.is_variable_type_valid(variable_name='seed',
                                   variable_value=settings['seed'],
                                   valid_types=(int))
    except TypeError as error:
        logger.error(error)
        raise
    
    ### 'algorithm'
    # Can be 'Isotropic3' or 'Riemann3'
    # SKS_VALID_ALGORITHM = ['Isotropic3', 'Riemann3']
    # try:
        # val.is_element_in_list()
    # if settings['algorithm'] not in SKS_ALGORITHM:
        
    
    # # Update the default travel cost dictionary
    # for parameter, cost in settings['costs'].items():
    #     DEFAULT_FMM_COSTS[parameter] = cost
        
    # # Complete the settings dictionary
    # for parameter, cost in DEFAULT_FMM_COSTS.items():
    #     if parameter not in settings['costs']:
    #         settings['costs'][parameter] = DEFAULT_FMM_COSTS[parameter]
    
    # 'mode'
    if settings['mode'] not in ['A', 'B', 'C', 'D']:
        pass
    
    # 'factors'
    
    return settings


#########################
### GEOLOGIC FEATURES ### ( Geology - Faults)
#########################

# def validate_geologic_feature_settings(settings: dict,
#                                        attribute: str,
#                                        grid) -> dict:
 
#     # If 'data' is empty
#     if isinstance(settings['data'], (str)) and (settings['data'] == ''):
#         settings['costs'] = {1: 0.4}
#         return settings
    
#     # Checks if type is str or numpy.ndarray
#     is_attribute_type_valid('data', settings['data'], (str, np.ndarray))

#     # Type is str
#     if isinstance(settings['data'], (str)):
#         path = settings['data']

#         # Checks if the datafile exist
#         is_path_valid(path, 'data')

#         # Checks if extension is valid
#         valid_extensions = ['gslib', 'npy', 'png', 'jpg', 'txt', 'csv', 'tif', 'tiff']
#         is_extension_valid(path, valid_extensions)
            
#         # Tries to open file
#         data = read_file(path, 'data')

#     # Type is np.ndarray
#     if isinstance(settings['data'], (np.ndarray)):
#         data = settings['data']

#     # Checks if data dimensions are valid
#     if isinstance(data, (np.ndarray)):
#         is_volume_dimensions_valid(data, attribute, grid, settings['axis'])
    
#     # Checks if provided 'costs' are valid
#     ids_data = np.unique(data)
#     if 'costs' in settings:
#         pass
#         # is_costs_dictionnary_valid(settings['costs'], ids_data)
#     else:
#         settings['costs'] = {}
#         if attribute == 'geology':
#             for i in range(len(ids_data)):
#                 settings['costs'][i + 1] = DEFAULT_FMM_COSTS['aquifer']
#         # elif attribute == 'fractures':  # TODO
#         #     for i in ids:
#         #         if i == 0:
#         #             continue
#         #         settings['costs'][i] = this.default_fmm_costs[self.label] + ((i - 1) * (1/100))
#         else:
#             settings['costs'][1] = DEFAULT_FMM_COSTS[attribute]
#         msg = "'costs' dictionary has been set automatically."
#         this.logger.warning(msg)
        
#     return settings


# def is_volume_dimensions_valid(array: np.ndarray, attribute: str, grid,
#                                axis: str) -> bool:
#     nx, ny, nz = grid.shape
#     axis_surface_size = {
#         'x': (ny * nz),
#         'y': (nx * nz),
#         'z': (nx * ny),
#     }
#     axis_surface_shapes = {
#         'x': (ny, nz),
#         'y': (nx, nz),
#         'z': (nx, ny),
#     }
    
#     data_size = len(array)
#     ### GSLIB ndarray case
#     if len(array.shape) == 1:
#         if data_size == nx * ny * nz:
#             return True
#         else:
#             msg = "The '{}' data dimensions do not match with the volume "
#             if data_size == axis_surface_size[axis]:
#                 msg += ("but match with the {}-side surface of the grid. Data "
#                         "will be replicated on {}-axis."
#                         .format(attribute, axis))
#                 this.logger.debug(msg)
#                 return True
#             else:
#                 msg += "neither with the {}-side surface of the grid (data : {}, grid_volume : {}, grid_{}_surface : {}).".format(attribute, data_size, nx * ny * nz, axis, axis_surface_size[axis])
#                 this.logger.critical(msg)
#                 raise ValueError(msg)

#     ### 2D ndarray case
#     elif len(array.shape) == 2:
#         if array.shape == axis_surface_shapes[axis]:
#             msg = ("The '{}' data dimensions match with the {}-side surface of"
#                    " the grid. Dat will be replicated on {}-axis."
#                    .format(attribute, axis, axis))
#             this.logger.debug(msg)
#             return True
#         else:
#             msg = ("The '{}' data dimensions do not match with the {}-side "
#                    "surface of the grid (data : {}, grid_volume : {}, "
#                    "grid_{}_surface : {})."
#                    .format(attribute, axis, array.shape, (nx, ny, nz),
#                            axis, axis_surface_shapes[axis]))
#             this.logger.critical(msg)
#             raise ValueError(msg)

#     ### 3D ndarray case
#     elif len(array.shape) == 3:
#         if array.shape == (nx, ny, nz):
#             return True
#         else:
#             msg = ("The '{}' data dimensions do not match neither with the"
#                    " volume nor the surface of the grid (data : {}, grid : {})"
#                    ".".format(attribute, array.shape, grid.shape))
#             this.logger.critical(msg)
#             raise ValueError(msg)
#     ### Other cases
#     else:
#         msg = ("The '{}' data dimensions do not match neither with the volume"
#                " nor the surface of the grid (data : {}, grid : {})."
#                .format(attribute, array.shape, grid.shape))
#         this.logger.critical(msg)
#         raise ValueError(msg)















































##############
### DOMAIN ###
##############

# def validate_domain_settings(settings: dict, grid) -> dict:
#     """TODO
    
#     - 'delimitation'
#     - 'topography'
#     - 'bedrock'
#     - 'water_table'
#     """
    
#     # Checks data validity for each subattribute
#     for attribute in ['delimitation', 'topography', 'bedrock', 'water_table']:
        
#         # Checks if data exists
#         test_a = isinstance(settings[attribute], (str))
#         if test_a and (settings[attribute] == ''):
#             continue
        
#         # Validates 'delimitation' attribute settings
#         if attribute == 'delimitation':
#             settings = validate_delimitation_settings(settings, grid)
#         else:
#             settings = validate_surface_settings(settings, attribute, grid)
    
#     return settings
 
 
# def validate_delimitation_settings(settings: dict, grid) -> dict:
#     ### Checks if type is valid
#     feature = 'delimitation'
#     valid_types = (str, list, np.ndarray)
#     is_attribute_type_valid(feature, settings[feature], valid_types)
    
#     # Type is str
#     if isinstance(settings[feature], (str)):
        
#         path = settings[feature]

#         # Checks if the datafile exist
#         is_path_valid(path, feature)

#         # Tries to open file
#         data = read_file(path, feature)
        
#     # Type is np.ndarray
#     elif isinstance(settings[feature], (np.ndarray)):
#         data = settings[feature].tolist()
        
#     # Type is list
#     if isinstance(settings[feature], (list)):
#         data = settings[feature]
        
#     ### Checks validity of data
#     # Checks if there is at least 3 points declared
#     is_list_length_valid(data, 3, feature)
    
#     # Checks if each vertex contains 2 coordinates
#     for vertex in data:
#         if len(vertex) != 2:
#             msg = ("The values of the 'delimitation' attribute contains at"
#                    "least one invalid vertex. Format must be like:"
#                    "[[x0, y0], ..., [xn, yn]].")
#             this.logger.critical(msg)
#             raise ValueError(msg)
    
#     # Checks type of coordinates
#     for vertex in data:
#         for coordinate in vertex:
#             is_coordinate_type_valid(coordinate, (int, float), feature)
            
#     # Checks if vertex are well in grid limits
#     validated_vertices = (
#         [k for k, (x, y) in enumerate(data)
#          if (grid.polygon.contains(Point(x, y)) is True)]
#     )
    
#     # Checks if there is enough vertices in order to create a delimitation
#     if len(validated_vertices) < 3:
#         msg = ("Not enough vertices inside the grid limits to create a"
#                "delimitation (3 minimum).")
#         this.logger.critical(msg)
#         raise ValueError(msg)

#     # Prints a warning when some points has not been validated
#     if len(validated_vertices) < len(data):
#         msg = ("{}/{} vertices validated inside the grid limits."
#                .format(len(validated_vertices), len(data)))
#         this.logger.warning(msg)
        
#     return settings


# def validate_surface_settings(settings: dict, attribute: str, grid) -> dict:
    
#     # Checks if type is valid
#     is_attribute_type_valid(attribute, settings[attribute], (str, np.ndarray))
    
#     # Type is str
#     if isinstance(settings[attribute], (str)):
#         path = settings[attribute]

#         # Checks if the datafile exist
#         is_path_valid(path, attribute)

#         # Tries to open file
#         data = read_file(path, attribute)
        
#     else:
#         data = settings[attribute]
        
#     # Checks if data surface is valid
#     is_surface_dimensions_valid(attribute, data, grid)
    
#     return settings
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
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
    
    # Test if 'seed' is of type int
    try:
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
            
    # # Control array size
    # if settings['data'].shape[1] not in [2, 3]:
    #     # TODO
    #     pass

    # # ### Checks validity of data

    # # # Checks if each points contains 2 or 3 coordinates
    # # for point in points:
    # #     if len(point) not in [2, 3]:
    # #         msg = "The values of the 'data' attribute contains at least one invalid point. Format must be like : [[x0, y0], ..., [xn, yn]] or [[x0, y0, z0], ..., [xn, yn, zn]]."
    # #         this.logger.critical(msg)
    # #         raise ValueError(msg)

    # # # Checks type of coordinates
    # # for point in points:
    # #     for coordinate in point:
    # #         if np.isnan(coordinate) or (not isinstance(coordinate, (int, float))):
    # #             msg = "The values of the 'data' attribute contains at least one invalid point. Coordinates must be of type int or float."
    # #             this.logger.critical(msg)
    # #             raise TypeError(msg)

    # # Checks if 'shuffle' is of type bool
    # is_attribute_type_valid('shuffle', settings['shuffle'], (bool))

    # # Checks if 'importance' is of type list
    # is_attribute_type_valid('importance', settings['importance'], (list))

    # # Checks if there is enough points to associate
    # if settings['number'] < len(settings['importance']):
    #     # TODO
    #     msg = "TODO - importance !"
    #     this.logger.critical(msg)
    #     raise ValueError(msg)
    
    # if attribute == 'inlets':
    #     # Checks if 'per_outlet' is of type list
    #     is_attribute_type_valid('per_outlet', settings['per_outlet'], (list))

    #     # Checks if there is enough points to associate
    #     if settings['number'] < len(settings['importance'])*len(settings['per_outlet']):
    #         # TODO - warning ou véritable erreur ???
    #         msg = "TODO - importance - per outlet !"
    #         this.logger.critical(msg)
    #         raise ValueError(msg)




    return settings


#################
### FRACTURES ###
#################

### TODO


# def validate_fractures_settings(settings: dict, grid) -> dict:
    
#     settings = validate_geologic_feature_settings(settings, 'fractures', grid)
    
#     if 'settings' in settings:
        
#         if (settings['settings'] == '') or (settings['settings'] == {}):
#             del settings['settings']
#         else:
#             costs = {}
#             for (i, frac_family) in enumerate(settings['settings']):
#                 for elem in ['alpha', 'density', 'orientation', 'dip', 'length']:
#                     if elem not in settings['settings'][frac_family]:
#                         msg = "The '{}' attribute in frac family '{}' is missing (required)".format(elem, frac_family)
#                         this.logger.critical(msg)
#                         raise KeyError(msg)
#                     if elem == 'density':
#                         settings['settings'][frac_family]['density'] = float(settings['settings'][frac_family]['density'])
#                 if 'cost' not in settings['settings'][frac_family]:
#                     msg = "The 'cost' attribute is missing (optional) and has been set with default value(s)."
#                     this.logger.warning(msg)
#                     settings['settings'][frac_family]['cost'] = DEFAULT_FMM_COSTS['fractures']
#                 costs[i+1] = settings['settings'][frac_family]['cost']
#             settings['costs'] = costs
        
#     return settings
