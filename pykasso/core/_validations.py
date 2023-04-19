"""
Input validations functions
"""

import PIL
import os
import sys
import logging
import numpy as np
from shapely.geometry import Point

############
### TODO ###
############
# validations verbosity
# revoir la VERBOSITY des messages
# message different quand la feature est vide ??

sks = sys.modules['pykasso.core.sks']
this = sys.modules[__name__]

this.ATTRIBUTES = {
    'sks': {
        'seed': ['optional', 0],
        'algorithm': ['optional', 'Isotropic3'],
        'costs': ['optional', {}],
        'mode': ['optional', 'A'],
    },
    'grid': {
        'x0': ['required', ''],
        'y0': ['required', ''],
        'z0': ['required', ''],
        'nx': ['required', ''],
        'ny': ['required', ''],
        'nz': ['required', ''],
        'dx': ['required', ''],
        'dy': ['required', ''],
        'dz': ['required', ''],
    },
    'domain': {
        'delimitation': ['optional', ''],
        'topography': ['optional', ''],
        'bedrock': ['optional', ''],
        'water_level': ['optional', ''],
    },
    'geology': {
        'data': ['optional', ''],
        'axis': ['optional', 'z'],
    },
    'beddings': {
        'data': ['optional', ''],
    },
    'faults': {
        'data': ['optional', ''],
        'axis': ['optional', 'z'],
    },
    'outlets': {
        'number': ['required', ''],
        'data': ['optional', ''],
        'shuffle': ['optional', False],
        'importance': ['required', []],
        'mode': ['optional', 'surface_down'],
        'geology': ['optional', None],
        'seed': ['optional', 0],
    },
    'inlets': {
        'number': ['required', ''],
        'data': ['optional', ''],
        'shuffle': ['optional', False],
        'importance': ['required', []],
        'per_outlet': ['required', []],
        'mode': ['optional', 'surface_up'],
        'geology': ['optional', None],
        'seed': ['optional', 0],
    },
    'tracers': {},
    'fractures': {
        'data': ['optional', ''],
        'axis': ['optional', 'z'],
        'seed': ['optional', 0],
    },
}

#################
### FUNCTIONS ###
#################

def validate_attribute_presence(settings, attribute, kind, default_value={}, is_subattribute=False):
    """
    TODO
    """
    if not is_subattribute:
        this.logger = logging.getLogger("{}.validation".format(attribute))

    if attribute not in settings:
        msg = "The '{}' attribute is missing ({}).".format(attribute, kind)
        if kind == 'required':
            this.logger.critical(msg)
            raise KeyError(msg)
        elif kind == 'optional':
            this.logger.warning(msg)
            settings[attribute] = default_value
    else:
        if settings[attribute] is None:
            msg = "The '{}' attribute is not defined ({}).".format(attribute, kind)
            if kind == 'required':
                this.logger.critical(msg)
                raise ValueError(msg)
            elif kind == 'optional':
                this.logger.warning(msg)
                settings[attribute] = default_value
        else:
            # TODO
            # Should we log ?
            pass
    return settings


def is_attribute_type_valid(attribute:str, value, types):
    """
    TODO
    """
    if not isinstance(value, types):
        msg = "The value of the '{}' attribute must be of type : {}".format(attribute, types)
        this.logger.critical(msg)
        raise TypeError(msg)
    else:
        return True


def is_path_valid(path:str, attribute:str) -> bool:
    """
    TODO
    """
    if not os.path.exists(path):
        msg = "The path from '{}' attribute is not valid. '{}' does not exist.".format(attribute, path)
        this.logger.error(msg)
        raise FileNotFoundError(msg)
    else:
        return True

def is_extension_valid(path:str, valid_extensions:list) -> bool:
    """
    TODO
    """
    extension = path.split('.')[-1]
    if extension not in valid_extensions:
        msg = "The extension '.{}' from '{}' location is not valid. Valid extensions : {}.".format(extension, path, valid_extensions)
        this.logger.error(msg)
        raise ValueError(msg)
    else:
        return True


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

        ### Others
        else:
            data = np.genfromtxt(path)

    except Exception as err:
        msg = "Impossible to read the file designated by the '{}' attribute. Location : {}".format(attribute, path)
        this.logger.error(msg)
        raise err
    else:
        return data


def is_attribute_value_valid(attribute: str, value: float, logical_test: str,
                             compared_to: float) -> bool:
    """
    TODO
    """
    logical_test_text = {
        '>': 'greater than',
        '>=': 'greater than or equal to',
        '<': 'less than',
        '<=': 'less than or equal to',
        '==': 'equal to',
        '!=': 'not equal to'
    }
    test = str(value) + logical_test + str(compared_to)
    if not eval(test):
        msg = "The value of the '{}' attribute must be {} {}.".format(attribute, logical_test_text[logical_test], compared_to)
        this.logger.critical(msg)
        raise ValueError(msg)
    else:
        return True


def is_list_length_valid(data: list, value: int, attribute: str) -> bool:
    if len(data) < value:
        msg = "'{}' data length is too short ({} elements minimum).".format(attribute, value)
        this.logger.critical(msg)
        raise ValueError(msg)
    else:
        return True

     
def is_coordinate_type_valid(coordinate: tuple, types: tuple,
                             attribute: str) -> bool:
    if not isinstance(coordinate, types):
        msg = "The values of the '{}' attribute contains at least one invalid vertex. Coordinates must be of type : {}.".format(attribute, types)
        this.logger.critical(msg)
        raise TypeError(msg)
    else:
        return True
                   

def is_surface_dimensions_valid(array: np.ndarray, grid) -> bool:
    nx, ny, nz = grid.shape
    if not (array.shape == (nx, ny)):
        msg = 'TODO - array shape : {} / nx, ny : {}'.format(array.shape, (nx, ny))
        this.logger.critical(msg)
        raise ValueError(msg)
    else:
        return True
    

def is_costs_dictionnary_valid(costs_dictionnary: dict, ids_data: list):
    """ """
    for i in ids_data:
        if i not in costs_dictionnary:
            msg = "The data id ({}) is not within 'costs' dictionnary keys ({})".format(i, list(costs_dictionnary.keys()))
            this.logger.error(msg)
            raise KeyError(msg)
    return True
    

####################################################################################

def validate_settings(feature: str, feature_settings: dict, grid=None) -> dict:
    """
    TODO
    """
    if feature == 'sks': 
        settings = validate_sks_settings(feature_settings)
    elif feature == 'grid':
        settings = validate_grid_settings(feature_settings)
    elif feature == 'domain':
        settings = validate_domain_settings(feature_settings, grid)
    elif feature in ['geology', 'beddings', 'faults']:
        settings = validate_geologic_feature_settings(feature_settings,
                                                      feature, grid)
    elif feature in ['inlets', 'outlets']:
        settings = validate_points_feature_settings(feature_settings, feature)
    elif feature == 'tracers':
        settings = validate_tracers_settings(feature_settings)
    elif feature == 'fractures':
        settings = validate_fractures_settings(feature_settings, grid)
        
    return settings

###########
### SKS ### (Seed - Algorithm - Costs)
###########

def validate_sks_settings(settings: dict) -> dict:
    # Checks attributes presence
    for attribute in this.ATTRIBUTES['sks']:
        kind, default_value = this.ATTRIBUTES['sks'][attribute]
        settings = validate_attribute_presence(settings, attribute, kind,
                                               default_value,
                                               is_subattribute=True)
        
    # Checks the presence of essential fmm costs
    costs = ['ratio', 'conduits', 'out']
    for cost in costs:
        if cost not in settings['costs']:
            settings['costs'][cost] = sks.default_fmm_costs[cost]
            
    return settings

############
### GRID ###
############

def validate_grid_settings(settings: dict) -> dict:
    # Checks attributes presence
    for attribute in this.ATTRIBUTES['grid']:
        kind, default_value = this.ATTRIBUTES['grid'][attribute]
        settings = validate_attribute_presence(settings, attribute, kind,
                                               default_value,
                                               is_subattribute=True)

    # Checks if the values of attributes are of type int or float
    for attribute in ['x0', 'y0', 'z0', 'dx', 'dy', 'dz']:
        is_attribute_type_valid(attribute, settings[attribute], (int, float))

    # Checks if the values of attributes are of type int
    for attribute in ['nx', 'ny', 'nz']:
        is_attribute_type_valid(attribute, settings[attribute], (int))

    # Checks if the values of attributes are well upper 0
    for attribute in ['nx', 'ny', 'nz']:
        is_attribute_value_valid(attribute, settings[attribute], '>', 0)

    return settings

##############
### DOMAIN ### (Delimitation - Topography - Bedrock Elevation - Water Level Elevation)
##############
def validate_domain_settings(settings:dict, grid) -> dict:
    # Checks attributes presence
    for attribute in this.ATTRIBUTES['domain']:
        kind, default_value = this.ATTRIBUTES['domain'][attribute]
        settings = validate_attribute_presence(settings, attribute, kind,
                                               default_value,
                                               is_subattribute=True)
    
    for attribute in ['delimitation', 'topography', 'bedrock', 'water_level']:
        
        # Checks if data exists
        if isinstance(settings[attribute], (str)) and (settings[attribute] == ''):
            continue
        
        # Validates 'delimitation' attribute settings
        if attribute == 'delimitation':
            settings = validate_delimitation_settings(settings, grid)
        else:
            settings = validate_surface_settings(settings, attribute, grid)
    
    return settings
 
def validate_delimitation_settings(settings: dict, grid) -> dict:
    ### Checks if type is valid
    is_attribute_type_valid('delimitation', settings['delimitation'], (str, list, np.ndarray))
    
    # Type is str
    if isinstance(settings['delimitation'], (str)):
        path = settings['delimitation']

        # Checks if the datafile exist
        is_path_valid(path, 'delimitation')

        # Tries to open file
        data = read_file(path, 'delimitation')
        
    # Type is np.ndarray
    if isinstance(settings['delimitation'], (np.ndarray)):
        data = settings['delimitation'].tolist()
        
    # Type is list
    if isinstance(settings['delimitation'], (list)):
        data = settings['delimitation']
        
    ### Checks validity of data
    # Checks if there is at least 3 points declared
    is_list_length_valid(data, 3, 'delimitation')
    
    # Checks if each vertex contains 2 coordinates
    for vertex in data:
        if len(vertex) != 2:
            msg = "The values of the 'delimitation' attribute contains at least one invalid vertex. Format must be like : [[x0, y0], ..., [xn, yn]]."
            this.logger.critical(msg)
            raise ValueError(msg)
    
    # Checks type of coordinates
    for vertex in data:
        for coordinate in vertex:
            is_coordinate_type_valid(coordinate, (int, float), 'delimitation')
            
    # Checks if vertex are well in grid limits
    validated_vertices = [k for k,(x,y) in enumerate(data) if (grid.polygon.contains(Point(x,y)) == True)]

    if len(validated_vertices) < 3:
        msg = "Not enough vertices inside the grid limits to create a delimitation (3 minimum)."
        this.logger.critical(msg)
        raise ValueError(msg)

    if len(validated_vertices) < len(data):
        msg = '{}/{} vertices validated inside the grid limits.'.format(len(validated_vertices), len(data))
        this.logger.warning(msg)
        
    return settings

def validate_surface_settings(settings: dict, attribute: str, grid) -> dict:
    
    # Checks if type is valid
    is_attribute_type_valid(attribute, settings[attribute], (str, np.ndarray))
    
    # Type is str
    if isinstance(settings[attribute], (str)):
        path = settings[attribute]

        # Checks if the datafile exist
        is_path_valid(path, attribute)

        # Tries to open file
        data = read_file(path, attribute)
        
    else:
        data = settings[attribute]
        
    # Checks if data surface is valid
    is_surface_dimensions_valid(data, grid)
    
    return settings

#########################
### GEOLOGIC FEATURES ### ( Geology - Beddings - Faults)
#########################

def validate_geologic_feature_settings(settings: dict, attribute: str,
                                       grid) -> dict:

    # Checks attributes presence
    for attribute_ in this.ATTRIBUTES[attribute]:
        kind, default_value = this.ATTRIBUTES[attribute][attribute_]
        settings = validate_attribute_presence(settings, attribute_, kind, default_value, is_subattribute=True)
 
    # If 'data' is empty
    if isinstance(settings['data'], (str)) and (settings['data'] == ''):
        settings['costs'] = {1 : 0.4}
        return settings
    
    # Checks if type is str or numpy.ndarray
    is_attribute_type_valid('data', settings['data'], (str, np.ndarray))

    # Type is str
    if isinstance(settings['data'], (str)):
        path = settings['data']

        # Checks if the datafile exist
        is_path_valid(path, 'data')

        # Checks if extension is valid
        valid_extensions = ['gslib', 'npy', 'png', 'jpg', 'txt', 'csv']
        is_extension_valid(path, valid_extensions)
            
        # Tries to open file
        data = read_file(path, 'data')

    # Type is np.ndarray
    if isinstance(settings['data'], (np.ndarray)):
        data = settings['data']

    # Checks if data dimensions are valid
    if attribute in ['beddings']:
        is_surface_dimensions_valid(data, grid)
    else:
        if isinstance(data, (np.ndarray)):
            is_volume_dimensions_valid(data, attribute, grid, settings['axis'])
    
    # Checks if provided 'costs' are valid
    ids_data = np.unique(data)
    if 'costs' in settings:
        is_costs_dictionnary_valid(settings['costs'], ids_data)
    else:
        settings['costs'] = {}
        if attribute == 'geology':
            for i in range(len(ids_data)):
                settings['costs'][i+1] = sks.default_fmm_costs['aquifer']
        # elif attribute == 'fractures': # TODO 
        #     for i in ids:
        #         if i == 0:
        #             continue
        #         settings['costs'][i] = this.default_fmm_costs[self.label] + ((i - 1) * (1/100))
        else:
            settings['costs'][1] = sks.default_fmm_costs[attribute]
        msg = "'costs' dictionary has been set automatically."
        this.logger.warning(msg)
        
    return settings


def is_volume_dimensions_valid(array:np.ndarray, attribute:str, grid, axis:str) -> bool:
    nx, ny, nz = grid.shape
    axis_surface_size = {
        'x' : (ny * nz),
        'y' : (nx * nz),
        'z' : (nx * ny),
    }
    axis_surface_shapes = {
        'x' : (ny, nz),
        'y' : (nx, nz),
        'z' : (nx, ny),
    }
    
    data_size = len(array)
    ### GSLIB ndarray case
    if len(array.shape) == 1:
        if data_size == nx * ny * nz:
            return True
        else:
            msg = "The '{}' data dimensions do not match with the volume "
            if data_size == axis_surface_size[axis]:
                msg += "but match with the {}-side surface of the grid. Data will be replicated on {}-axis.".format(attribute, axis, axis)
                this.logger.debug(msg)
                return True
            else:
                msg += "neither with the {}-side surface of the grid (data : {}, grid_volume : {}, grid_{}_surface : {}).".format(attribute, axis, data_size, nx * ny * nz, axis, axis_surface_size[axis])
                this.logger.critical(msg)
                raise ValueError(msg)

    ### 2D ndarray case
    elif len(array.shape) == 2:
        if array.shape == axis_surface_shapes[axis]:
            msg = "The '{}' data dimensions match with the {}-side surface of the grid. Dat will be replicated on {}-axis.".format(attribute, axis, axis)
            this.logger.debug(msg)
            return True
        else:
            msg = "The '{}' data dimensions do not match with the {}-side surface of the grid (data : {}, grid_volume : {}, grid_{}_surface : {}).".format(attribute, axis, array.shape, (nx,ny,nz), axis, axis_surface_shapes[axis])
            this.logger.critical(msg)
            raise ValueError(msg)

    ### 3D ndarray case
    elif len(array.shape) == 3:
        if array.shape == (nx, ny, nz):
            return True
        else:
            msg = "The '{}' data dimensions do not match neither with the volume nor the surface of the grid (data : {}, grid : {}).".format(attribute, array.shape, grid.shape)
            this.logger.critical(msg)
            raise ValueError(msg)
    ### Other cases
    else:
        msg = "The '{}' data dimensions do not match neither with the volume nor the surface of the grid (data : {}, grid : {}).".format(attribute, array.shape, grid.shape)
        this.logger.critical(msg)
        raise ValueError(msg)
    
########################
### OUTLETS - INLETS ###
########################
def validate_points_feature_settings(settings:dict, attribute:str) -> dict:
    """
    TODO
    """
    points = ''

    # Checks attributes presence
    for attribute_ in this.ATTRIBUTES[attribute]:
        kind, default_value = this.ATTRIBUTES[attribute][attribute_]
        settings = validate_attribute_presence(settings, attribute_, kind, default_value, is_subattribute=True)
        
    # Checks if 'number' is of type int
    is_attribute_type_valid('number', settings['number'], (int))

    # Checks if 'number' attribute value is valid
    is_attribute_value_valid('number', settings['number'], '>', 0)

    # Checks data when provided
    if not ((settings['data'] == []) or (settings['data'] == '')):

        # Checks if data type is valid:
        is_attribute_type_valid('data', settings['data'], (str, list))

        # Type is str:
        if isinstance(settings['data'], (str)):
            path = settings['data']

            # Checks if the datafile exist
            is_path_valid(path, 'data')

            # Tries to open file
            points = read_file(path, 'data')
            if len(points.shape) == 1: # If points contains only one element
                points = np.array([points])

        # Type is list:
        if isinstance(settings['data'], (list)):
            points = settings['data']

        ### Checks validity of data

        # Checks if each points contains 2 or 3 coordinates
        for point in points:
            if len(point) not in [2, 3]:
                msg = "The values of the 'data' attribute contains at least one invalid point. Format must be like : [[x0, y0], ..., [xn, yn]] or [[x0, y0, z0], ..., [xn, yn, zn]]."
                this.logger.critical(msg)
                raise ValueError(msg)

        # Checks type of coordinates
        for point in points:
            for coordinate in point:
                if np.isnan(coordinate) or (not isinstance(coordinate, (int, float))):
                    msg = "The values of the 'data' attribute contains at least one invalid point. Coordinates must be of type int or float."
                    this.logger.critical(msg)
                    raise TypeError(msg)

    # Checks if 'shuffle' is of type bool
    is_attribute_type_valid('shuffle', settings['shuffle'], (bool))

    # Checks if 'importance' is of type list
    is_attribute_type_valid('importance', settings['importance'], (list))

    # Checks if there is enough points to associate
    if settings['number'] < len(settings['importance']):
        # TODO
        msg = "TODO - importance !"
        this.logger.critical(msg)
        raise ValueError(msg)
    
    if attribute == 'inlets':
        # Checks if 'per_outlet' is of type list
        is_attribute_type_valid('per_outlet', settings['per_outlet'], (list))

        # Checks if there is enough points to associate
        if settings['number'] < len(settings['importance'])*len(settings['per_outlet']):
            # TODO - warning ou véritable erreur ???
            msg = "TODO - importance - per outlet !"
            this.logger.critical(msg)
            raise ValueError(msg)

    return settings


### TODO
def validate_tracers_settings(settings:dict) -> dict:
    return settings

#################
### FRACTURES ###
#################

### TODO
def validate_fractures_settings(settings:dict, grid) -> dict:
    
    settings = validate_geologic_feature_settings(settings, 'fractures', grid)
    
    if 'settings' in settings:
        
        if (settings['settings'] == '') or (settings['settings'] == {}):
            del settings['settings']
        else:
            costs = {}
            for (i, frac_family) in enumerate(settings['settings']):
                for elem in ['alpha', 'density', 'orientation', 'dip', 'length']:
                    if elem not in settings['settings'][frac_family]:
                        msg = "The '{}' attribute in frac family '{}' is missing (required)".format(elem, frac_family)
                        this.logger.critical(msg)
                        raise KeyError(msg)
                    if elem == 'density':
                        settings['settings'][frac_family]['density'] = float(settings['settings'][frac_family]['density'])
                if 'cost' not in settings['settings'][frac_family]:
                    msg = "The 'cost' attribute is missing (optional) and has been set with default value(s)."
                    this.logger.warning(msg)
                    settings['settings'][frac_family]['cost'] = sks.default_fmm_costs['fractures']
                costs[i+1] = settings['settings'][frac_family]['cost']
            settings['costs'] = costs
        
    return settings
