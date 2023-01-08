"""
Input validations functions
"""

import os
import sys
import logging
import numpy as np

# import _exceptions
# import core_settings

############
### TODO ###
############
# VERBOSITY MGMT
# Count the warnings ?
# message different quand la feature est vide
# validation en fonction de 'axis'
# msg = "The '{}' attribute was missing. Topography will be set from geological data.".format(attribute)
# geologic_feature : Karst
# geologic_feature : Field
# validate_tracers_settings()
# validate_fractures_settings()
# validate_seed ?


this = sys.modules[__name__]
this.logger = logging.getLogger("validation.settings")

this.ERRORS   = []
this.WARNINGS = 0 

this.MANDATORY_ATTRIBUTES_SKS = {
    'grid' : ['x0', 'y0', 'z0', 'nx', 'ny', 'nz', 'dx', 'dy', 'dz'],
    'fmm'  : [],
}
this.OPTIONAL_ATTRIBUTES_SKS  = {
    'mask' : {
        'data' : ''
    },
    'topography' : {
        'data' : ''
    },
    'geology' : {
        'data' : '',
        'axis' : 'z',
        'cost' : '',
    },
    'faults' : {
        'data' : '',
        'axis' : 'z',
        'cost' : '',
    },
}
this.MANDATORY_ATTRIBUTES_SIM = {
    'seed'    : [],
    'outlets' : {
        'number'     : 1,
        'data'       : '',
        'shuffle'    : False,
        'importance' : '',
    },
    'inlets'  : {
        'number'     : 1,
        'data'       : '',
        'shuffle'    : False,
        'importance' : '',
        'per_outlet' : '',
    },
}
this.OPTIONAL_ATTRIBUTES_SIM  = {
    'tracers'   : {},
    'fractures' : {},
}

def validate_settings_structure(settings, kind):
    """
    TODO
    """
    if kind == 'SKS':
        this.mandatory_attributes = this.MANDATORY_ATTRIBUTES_SKS
        this.optional_attributes  = this.OPTIONAL_ATTRIBUTES_SKS
    if kind == 'SIM':
        this.mandatory_attributes = this.MANDATORY_ATTRIBUTES_SIM
        this.optional_attributes  = this.OPTIONAL_ATTRIBUTES_SIM
    
    for attribute in this.mandatory_attributes:
        settings = validate_attribute_presence(settings, attribute)

    for attribute in this.optional_attributes:
        settings = validate_attribute_presence(settings, attribute, 'optional')

    return settings

#################
### FUNCTIONS ###
#################

def validate_attribute_presence(settings, attribute, kind='mandatory'):
    """
    TODO
    """
    if attribute not in settings:
        msg = "The '{}' attribute is missing ({}).".format(attribute, kind)
        if kind == 'mandatory':
            this.logger.critical(msg)
            this.ERRORS.append(KeyError(msg))
        if kind == 'optional':
            this.logger.warning(msg)
            settings[attribute] = this.optional_attributes[attribute]
    
    if (attribute in settings) and (settings[attribute] is None):
        msg = "The value of the '{}' attribute is not defined.".format(attribute, kind)
        if kind == 'mandatory':
            this.logger.critical(msg)
            this.ERRORS.append(KeyError(msg))
        if kind == 'optional':
            this.logger.warning(msg)
            settings[attribute] = this.optional_attributes[attribute]

    return settings


def is_attribute_type_valid(attribute, value, types):
    """
    TODO
    """
    if not isinstance(value, types):
        msg = "The value of the '{}' attribute must be of type : {}".format(attribute, types)
        this.logger.critical(msg)
        this.ERRORS.append(TypeError(msg))
        return False
    else:
        return True


def is_path_valid(path, attribute):
    """
    TODO
    """
    if not os.path.exists(path):
        msg = "The path from '{}' attribute is not valid. '{}' does not exist.".format(attribute, path)
        this.logger.error(msg)
        this.ERRORS.append(FileNotFoundError(msg))
        return False
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
        this.ERRORS.append(ValueError(msg))
        return False
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
        this.ERRORS.append(err)
        return None

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
        this.ERRORS.append(ValueError(msg))
        return False
    else:
        return True


def is_state_valid(feature, log=False):
    """
    TODO
    """
    if len(this.ERRORS) == 0:
        if log:
            this.logger.info('The {} settings has been validated'.format(feature))
        return True
    else:
        this.logger.error("{} error(s) detected. The {} settings has not been validated".format(len(this.ERRORS), feature))
        raise this.ERRORS[0]


################################################################################################
### CORE ###
############

###############################################################################################
### SKS ###
###########

############
### GRID ###
############
def validate_grid_settings(settings):
    """
    TODO
    """
    this.logger = logging.getLogger("validation.grid")

    # Checks attributes presence
    attributes = this.MANDATORY_ATTRIBUTES_SKS['grid']
    for attribute in attributes:
        settings = validate_attribute_presence(settings, attribute)

    # Checks validity of attributes
    if is_state_valid('grid'):

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

    # Prints state of the validation process
    is_state_valid('grid', True)

    return settings


############
### MASK ###
############
def validate_mask_settings(settings):
    """
    TODO
    """
    this.logger = logging.getLogger("validation.mask")
    data = ''

    # Checks attributes presence
    this.optional_attributes = this.OPTIONAL_ATTRIBUTES_SKS['mask']
    for attribute in this.optional_attributes:
        settings = validate_attribute_presence(settings, attribute, 'optional')

    # Checks validity of attributes
    if is_state_valid('mask'):

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

    # Checks validity of data
    if is_state_valid('mask') and (not isinstance(data, (str))):

        # Checks if there is at least 3 points declared
        if len(data) < 3:
            msg = "Not enough vertices to create a mask delimitation (3 minimum)."
            this.logger.error(msg)
            this.ERRORS.append(ValueError(msg))

        # Checks if each vertex contains 2 coordinates
        for vertex in data:
            if len(vertex) != 2:
                msg = "The values of the 'data' attribute contains at least one invalid vertex. Format must be like : [[x0, y0], ..., [xn, yn]]."
                this.logger.error(msg)
                this.ERRORS.append(ValueError(msg))
                break
            
        # Checks type of coordinates
        for vertex in data:
            for coordinate in vertex:
                if np.isnan(coordinate) or (not isinstance(coordinate, (int, float))):
                    msg = "The values of the 'data' attribute contains at least one invalid vertex. Coordinates must be of type int or float."
                    this.logger.error(msg)
                    this.ERRORS.append(TypeError(msg))
                    break

    # Controls validity of attributes
    is_state_valid('mask', True)

    return settings


#########################
### GEOLOGIC FEATURES ###
#########################
# TOPOGRAPHY
# GEOLOGY
# FAULTS
# TODO - FRACTURES ?
# TODO - KARST
# TODO - FIELD

def validate_geologic_feature_settings(geologic_feature, settings, grid):
    """
    TODO
    """
    this.logger = logging.getLogger("validation.{}".format(geologic_feature))

    # Checks attributes presence
    this.optional_attributes = this.OPTIONAL_ATTRIBUTES_SKS[geologic_feature]
    for attribute in this.optional_attributes:
        settings = validate_attribute_presence(settings, attribute, 'optional')

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

            # Checks validity of data
            if is_state_valid(geologic_feature):

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
                                this.logger.error(msg)
                                this.ERRORS.append(ValueError(msg))
                            else:
                                msg = "The '{}' data dimensions do not match with the volume but with the surface of the grid. Data will be replicated on z-axis.".format(geologic_feature)
                                this.logger.warning(msg)

                    else:
                        if data.shape != (nx, ny, nz):
                            if data.shape[0:2] != (nx, ny):
                                msg = "The '{}' data dimensions do not match neither with the volume nor the surface of the grid (data : {}, grid : {}).".format(geologic_feature, data.shape, (nx, ny, nz))
                                this.logger.error(msg)
                                this.ERRORS.append(ValueError(msg))
                            else:
                                msg = "The '{}' data dimensions do not match with the volume but with the surface of the grid. Data will be replicated on z-axis.".format(geologic_feature)
                                this.logger.warning(msg)

    # Controls validity of attributes
    is_state_valid(geologic_feature, True)

    return settings


###############################################################################################
### SIM ###
###########

########################
### OUTLETS - INLETS ###
########################
def validate_points_feature_settings(points_feature, settings):
    """
    TODO
    """
    this.logger = logging.getLogger("validation.{}".format(points_feature))
    points = ''

    # Checks attributes presence
    this.optional_attributes = this.MANDATORY_ATTRIBUTES_SIM[points_feature]
    for attribute in this.optional_attributes:
        settings = validate_attribute_presence(settings, attribute)

    # Checks validity of attributes
    if is_state_valid(points_feature):
        
        # Checks if 'number' is of type int
        is_attribute_type_valid('number', settings['number'], (int))

        # Checks if 'number' attribute value is valid
        is_attribute_value_valid('number', settings['number'], '>', 0)

        # TODO 
        # 'shuffle'
        # 'importance'
        # 'per_outlet'

        # Handles '' in 'data' attributes
        if settings['data'] == '':
            settings['data'] == []

        # Checks data when provided
        if settings['data'] != []:

            # Checks if data type is valid:
            if is_attribute_type_valid('data', settings['data'], (str, list)):

                # Type is str:
                if isinstance(settings['data'], str):
                    path = settings['data']

                    # Checks if the datafile exist
                    if is_path_valid(path, 'data'):

                        # Tries to open file
                        points = read_file(path, 'data')

                # Type is list:
                if isinstance(settings['data'], list):
                    points = settings['data']

                # Checks validity of data
                if is_state_valid(points_feature):

                    # Checks if each points contains 2 or 3 coordinates
                    for point in points:
                        if len(point) not in [2, 3]:
                            msg = "The values of the 'data' attribute contains at least one invalid point. Format must be like : [[x0, y0], ..., [xn, yn]] or [[x0, y0, z0], ..., [xn, yn, zn]]."
                            this.logger.error(msg)
                            this.ERRORS.append(ValueError(msg))
                            break

                    # Checks type of coordinates
                    for point in points:
                        for coordinate in point:
                            if np.isnan(coordinate) or (not isinstance(coordinate, (int, float))):
                                msg = "The values of the 'data' attribute contains at least one invalid point. Coordinates must be of type int or float."
                                this.logger.error(msg)
                                this.ERRORS.append(TypeError(msg))
                                break

    # Prints state of the validation process
    is_state_valid(points_feature, True)

    return settings


# TODO 
def validate_tracers_settings(settings):
    """
    TODO
    """
    return settings


# TODO
def validate_fractures_settings(settings):
    """
    TODO
    """
    return settings