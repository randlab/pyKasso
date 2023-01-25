"""
Input validations functions
"""

import os
import sys
import logging
import numpy as np

############
### TODO ###
############
# VERBOSITY MGMT
# message different quand la feature est vide ?? 
# validation en fonction de 'axis'
# msg = "The '{}' attribute was missing. Topography will be set from geological data.".format(attribute)
# geologic_feature : Karst
# geologic_feature : Field
# validate_tracers_settings()
# validate_fractures_settings()


this = sys.modules[__name__]

this.ATTRIBUTES = {
    'sks' : {
        'seed' : ['optional', 0],
    },
    'grid' : {
        'x0' : ['required', ''],
        'y0' : ['required', ''],
        'z0' : ['required', ''],
        'nx' : ['required', ''],
        'ny' : ['required', ''],
        'nz' : ['required', ''],
        'dx' : ['required', ''],
        'dy' : ['required', ''],
        'dz' : ['required', ''],
    },
    'mask' : {
        'data' : ['optional', ''],
    },
    'topography' : {
        'data' : ['optional', ''],
    },
    'geology' : {
        'data' : ['optional', ''],
        'axis' : ['optional', 'z'],
        'cost' : ['optional', {1 : 0.4}],
    },
    'faults' : {
        'data' : ['optional', ''],
        'axis' : ['optional', 'z'],
        'cost' : ['optional', {1 : 0.2}],
    },
    'outlets' : {
        'number'     : ['required', ''],
        'data'       : ['optional', ''],
        'shuffle'    : ['optional', False],
        'importance' : ['required', []],
        'x'          : ['optional', 'lambda : grid.xmin + grid.nx * rng.random() * grid.dx'],
        'y'          : ['optional', 'lambda : grid.ymin + grid.ny * rng.random() * grid.dy'],
        'z'          : ['optional', 'lambda : grid.zmin + grid.dz/1000'],
        'geology'    : ['optional', None],
        'seed'       : ['optional', 0],
    },
    'inlets'  : {
        'number'     : ['required', ''],
        'data'       : ['optional', ''],
        'shuffle'    : ['optional', False],
        'importance' : ['required', []],
        'per_outlet' : ['required', []],
        'x'          : ['optional', 'lambda : grid.xmin + grid.nx * rng.random() * grid.dx'],
        'y'          : ['optional', 'lambda : grid.ymin + grid.ny * rng.random() * grid.dy'],
        'z'          : ['optional', 'lambda : grid.zmax - grid.dz/1000'],
        'geology'    : ['optional', None],
        'seed'       : ['optional', 0],
    },
    'tracers'   : {},
    'fractures' : {
        'data'     : ['optional', ''],
        'seed'     : ['optional', 0],
        'cost'     : ['optional', {1 : 0.2}],
        'settings' : ['optional', ''],
        
    },
    'fmm' : {
        'algorithm' : ['optional', 'Isotropic3'],
        'cost'      : ['optional', {'ratio' : 0.5}]
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


def is_attribute_type_valid(attribute, value, types):
    """
    TODO
    """
    if not isinstance(value, types):
        msg = "The value of the '{}' attribute must be of type : {}".format(attribute, types)
        this.logger.critical(msg)
        raise TypeError(msg)
    else:
        return True


def is_path_valid(path, attribute):
    """
    TODO
    """
    if not os.path.exists(path):
        msg = "The path from '{}' attribute is not valid. '{}' does not exist.".format(attribute, path)
        this.logger.error(msg)
        raise FileNotFoundError(msg)
    else:
        return True

def is_extension_valid(path, valid_extensions):
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


def read_file(path, attribute):
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

        # Images
        elif extension in ['jpg', 'png']:
            return None

        ### Others
        else:
            data = np.genfromtxt(path)

    except Exception as err:
        msg = "Impossible to read the file designated by the '{}' attribute. Location : {}".format(attribute, path)
        this.logger.error(msg)
        raise err
    else:
        return data


def is_attribute_value_valid(attribute, value, logical_test, compared_to):
    """
    TODO
    """
    logical_test_text = {
        '>'  : 'greater than',
        '>=' : 'greater than or equal to',
        '<'  : 'less than',
        '<=' : 'less than or equal to',
        '==' : 'equal to',
        '!=' : 'not equal to'
    }
    test = str(value) + logical_test + str(compared_to)
    if not eval(test):
        msg = "The value of the '{}' attribute must be {} {}.".format(attribute, logical_test_text[logical_test], compared_to)
        this.logger.critical(msg)
        raise ValueError(msg)
    else:
        return True


################################################################################################
### CORE ###
############

###############################################################################################
### SKS ###
###########

def validate_settings(feature, feature_settings, grid=None):
    """ 
    TODO
    """
    if feature == 'grid':
        settings = validate_grid_settings(feature_settings)
    elif feature == 'mask':
        settings = validate_mask_settings(feature_settings, grid)
    elif feature in ['topography', 'geology', 'faults']:
        settings = validate_geologic_feature_settings(feature, feature_settings, grid)
    elif feature == 'sks':
        settings = validate_seed_settings(feature_settings)
    elif feature in ['inlets', 'outlets']:
        settings = validate_points_feature_settings(feature, feature_settings)
    # elif feature == 'tracers':
    #     settings = validate_tracers_settings()
    elif feature == 'fractures':
        settings = validate_fractures_settings(feature_settings)
    elif feature == 'fmm':
        settings = validate_fmm_settings(feature_settings)
        
    return settings

############
### GRID ###
############
def validate_grid_settings(settings):
    """
    TODO
    """
    # Checks attributes presence
    for attribute in this.ATTRIBUTES['grid']:
        settings = validate_attribute_presence(settings, attribute, this.ATTRIBUTES['grid'][attribute][0], is_subattribute=True)

    # TODO
    # Checks global size of the model
    # grid_surface = sks_settings[attribute]['nx'] * sks_settings[attribute]['ny']
    # grid_volume  = grid_surface * sks_settings[attribute]['nz']
    # if grid_volume > core_settings.grid_size_limit:

    # Checks if the values of attributes are of type int or float
    attributes = ['x0', 'y0', 'z0', 'dx', 'dy', 'dz']
    for attribute in attributes:
        is_attribute_type_valid(attribute, settings[attribute], (int, float))

    # Checks if the values of attributes are of type int
    attributes = ['nx', 'ny', 'nz']
    for attribute in attributes:
        is_attribute_type_valid(attribute, settings[attribute], (int))

    # Checks if the values of attributes are well upper 1
    attributes = ['nx', 'ny', 'nz']
    for attribute in attributes:
        is_attribute_value_valid(attribute, settings[attribute], '>', 0)

    return settings


############
### MASK ###
############
def validate_mask_settings(settings, grid):
    """
    TODO
    """
    # Checks attributes presence
    for attribute in this.ATTRIBUTES['mask']:
        kind, default_value = this.ATTRIBUTES['mask'][attribute]
        settings = validate_attribute_presence(settings, attribute, kind, default_value, is_subattribute=True)

    # If 'data' is not empty
    if settings['data'] != '':

        # Checks if type is str or list
        if is_attribute_type_valid('data', settings['data'], (str, list)):

            # Type is str
            if isinstance(settings['data'], (str)):
                path = settings['data']

                # Checks if the datafile exist
                if is_path_valid(path, 'data'):

                    # Tries to open file
                    data = read_file(path, 'data')
            
            # Type is list
            if isinstance(settings['data'], (list)):
                data = settings['data']

        ### Checks validity of data
        # Checks if there is at least 3 points declared
        if len(data) < 3:
            msg = "Not enough vertices to create a mask delimitation (3 minimum)."
            this.logger.critical(msg)
            raise ValueError(msg)

        # Checks if each vertex contains 2 coordinates
        for vertex in data:
            if len(vertex) != 2:
                msg = "The values of the 'data' attribute contains at least one invalid vertex. Format must be like : [[x0, y0], ..., [xn, yn]]."
                this.logger.critical(msg)
                raise ValueError(msg)
            
        # Checks type of coordinates
        for vertex in data:
            for coordinate in vertex:
                if not isinstance(coordinate, (int, float)):
                    msg = "The values of the 'data' attribute contains at least one invalid vertex. Coordinates must be of type int or float."
                    this.logger.critical(msg)
                    raise TypeError(msg)

        # Checks if vertex are well in grid limits
        validated_vertices = [k for k,(x,y) in enumerate(data) if (grid.path.contains_point((x,y)) == True)]

        if len(validated_vertices) < 3:
            msg = "Not enough vertices inside the grid limits to create a mask delimitation (3 minimum)."
            this.logger.critical(msg)
            raise ValueError(msg)

        if len(validated_vertices) < len(data):
            msg = '{}/{} vertices validated inside the grid limits.'.format(len(validated_vertices), len(data))
            this.logger.warning(msg)

    return settings


#########################
### GEOLOGIC FEATURES ###
#########################
# TODO - KARST
# TODO - FIELD

def validate_geologic_feature_settings(geologic_feature, settings, grid):
    """
    TODO
    """
    # Checks attributes presence
    for attribute in this.ATTRIBUTES[geologic_feature]:
        kind, default_value = this.ATTRIBUTES[geologic_feature][attribute]
        settings = validate_attribute_presence(settings, attribute, kind, default_value, is_subattribute=True)

    # If 'data' is not empty
    if settings['data'] != '':

        # Checks if type is str or numpy.ndarray
        if is_attribute_type_valid('data', settings['data'], (str, np.ndarray)):

            # Type is str
            if isinstance(settings['data'], (str)):
                path = settings['data']

                # Checks if the datafile exist
                if is_path_valid(path, 'data'):

                    # Checks if extension is valid
                    valid_extensions = ['gslib', 'npy', 'png', 'jpg']
                    if is_extension_valid(path, valid_extensions):
                        
                        # Tries to open file
                        data = read_file(path, 'data')

            # Type is np.ndarray
            if isinstance(settings['data'], (np.ndarray)):
                data = settings['data']

        ### Checks if data dimensions are valid
        nx, ny, nz = grid.nx, grid.ny, grid.nz

        # Numpy array case
        if isinstance(data, (np.ndarray)):
            # GSLIB case
            if len(data.shape) == 1:
                data_size = len(data)
                if data_size != nx * ny * nz:
                    if data_size != nx * ny:
                        msg = "The '{}' data dimensions do not match neither with the volume nor the surface of the grid (data : {}, grid : {})".format(attribute, data_size, nx * ny * nz)
                        this.logger.critical(msg)
                        raise ValueError(msg)
                    else:
                        msg = "The '{}' data dimensions do not match with the volume but with the surface of the grid. Data will be replicated on z-axis.".format(geologic_feature)
                        this.logger.warning(msg)

            else:
                if data.shape != (nx, ny, nz):
                    if data.shape[0:2] != (nx, ny):
                        msg = "The '{}' data dimensions do not match neither with the volume nor the surface of the grid (data : {}, grid : {}).".format(geologic_feature, data.shape, (nx, ny, nz))
                        this.logger.critical(msg)
                        raise ValueError(msg)
                    else:
                        msg = "The '{}' data dimensions do not match with the volume but with the surface of the grid. Data will be replicated on z-axis.".format(geologic_feature)
                        this.logger.warning(msg)

    return settings

#############
### SEEDS ###
#############
# todo - rename ?
def validate_seed_settings(settings):
    """ 
    TODO
    """
    settings = validate_attribute_presence(settings, 'seed', 'optional', default_value=0, is_subattribute=True)
    return settings

########################
### OUTLETS - INLETS ###
########################
def validate_points_feature_settings(points_feature, settings):
    """
    TODO
    """
    points = ''

    # Checks attributes presence
    for attribute in this.ATTRIBUTES[points_feature]:
        kind, default_value = this.ATTRIBUTES[points_feature][attribute]
        settings = validate_attribute_presence(settings, attribute, kind, default_value, is_subattribute=True)
        
    # Checks if 'number' is of type int
    is_attribute_type_valid('number', settings['number'], (int))

    # Checks if 'number' attribute value is valid
    is_attribute_value_valid('number', settings['number'], '>', 0)

    # Checks data when provided
    if not ((settings['data'] == []) or (settings['data'] == '')):

        # Checks if data type is valid:
        if is_attribute_type_valid('data', settings['data'], (str, list)):

            # Type is str:
            if isinstance(settings['data'], (str)):
                path = settings['data']

                # Checks if the datafile exist
                if is_path_valid(path, 'data'):

                    # Tries to open file
                    points = read_file(path, 'data')

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
    
    if points_feature == 'inlets':
        # Checks if 'per_outlet' is of type list
        is_attribute_type_valid('per_outlet', settings['per_outlet'], (list))

        # Checks if there is enough points to associate
        if settings['number'] < len(settings['importance'])*len(settings['per_outlet']):
            # TODO - warning ou véritable erreur ???
            msg = "TODO - importance - per outlet !"
            this.logger.critical(msg)
            raise ValueError(msg)

    return settings


# TODO 
def validate_tracers_settings(settings):
    """
    TODO
    """
    return settings

#################
### FRACTURES ###
#################

# TODO
def validate_fractures_settings(settings):
    """
    TODO
    """

    # Checks attributes presence
    for attribute in this.ATTRIBUTES['fractures']:
        kind, default_value = this.ATTRIBUTES['fractures'][attribute]
        settings = validate_attribute_presence(settings, attribute, kind, default_value, is_subattribute=True)

    return settings


###########
### FMM ###
###########

def validate_fmm_settings(settings) :
    """
    TODO
    """

    # Checks attributes presence
    for attribute in this.ATTRIBUTES['fmm']:
        kind, default_value = this.ATTRIBUTES['fmm'][attribute]
        settings['fmm'] = validate_attribute_presence(settings['fmm'], attribute, kind, default_value, is_subattribute=True)

    # TODO - develop more test
    # if len(settings['inlets']['per_outlet']) != settings['outlets']['number']:
    #     msg = "_validate_fmm_settings... TODO"
    #     this.logger.critical(msg)
    #     raise ValueError(msg)

    return settings