"""
Input validations functions
"""

import os
import logging
import numpy as np

# import _exceptions
# import core_settings

ERRORS   = []
WARNINGS = []

def validate_core_settings(core_settings:dict):
    """
    DOCSTRING

    Inputs
    ------

    Returns
    -------
    """
    global ERRORS

    logger = logging.getLogger("validations - core")
    msg = "The core settings file has been validated."
    logger.info(msg)
    return core_settings


def validate_sks_settings(sks_settings:dict):
    """
    DOCSTRING

    Inputs
    ------

    Returns
    -------
    """
    global ERRORS
    logger = logging.getLogger("validations - sks")

    ############
    ### GRID ###
    ############
    attribute = 'grid'

    # Checks presence of the 'grid' attribute
    if attribute not in sks_settings:
        msg = "The '{}' attribute is missing.".format(attribute)
        logger.error(msg)
        ERRORS.append(KeyError(msg))
    else:
        # Checks presence of the subattributes of the 'grid' attribute
        subattributes = ['x0', 'y0', 'z0', 'nx', 'ny', 'nz', 'dx', 'dy', 'dz']
        for subattribute in subattributes:
            if subattribute not in sks_settings[attribute]:
                msg = "The '{}' subattribute from the '{}' attribute is missing.".format(subattribute, attribute)
                logger.error(msg)
                ERRORS.append(KeyError(msg))

    if len(ERRORS) == 0:

        # Checks global size of the model
        grid_surface = sks_settings[attribute]['nx'] * sks_settings[attribute]['ny']
        grid_volume  = grid_surface * sks_settings[attribute]['nz']
        # if grid_volume > core_settings.grid_size_limit:
        #     pass
        #     TODO
    
        # Checks if subattributes are of type int or float
        subattributes = ['x0', 'y0', 'z0', 'dx', 'dy', 'dz']
        for subattribute in subattributes:
            if not isinstance(sks_settings[attribute][subattribute], (int, float)):
                msg = "The '{}' subattribute from the '{}' attribute must be of type int or float.".format(subattribute, attribute)
                logger.error(msg)
                ERRORS.append(TypeError(msg))

        # Checks if some subattributes are type int
        subattributes = ['nx', 'ny', 'nz']
        for subattribute in subattributes:
            if not isinstance(sks_settings[attribute][subattribute], int):
                msg = "The '{}' subattribute from the '{}' attribute must be of type int.".format(subattribute, attribute)
                logger.error(msg)
                ERRORS.append(TypeError(msg))

        # Checks if the values of some attributes are well upper 1
        subattributes = ['nx', 'ny', 'nz']
        for subattribute in subattributes:
            if sks_settings[attribute][subattribute] < 1:
                msg = "The '{}' subattribute from the '{}' attribute cannot be lower as 1.".format(subattribute, attribute)
                logger.error(msg)
                ERRORS.append(ValueError(msg))


    ############
    ### MASK ###
    ############
    attribute = 'mask'
    data = ''

    # Checks presence of the 'mask' attribute
    if attribute not in sks_settings:
        msg = "The '{}' attribute was missing.".format(attribute)
        logger.warning(msg)
        sks_settings[attribute] = {
            'data' : ''
        }

    # Checks presence of the subattributes of the 'mask' attribute
    subattributes = ['data']
    for subattribute in subattributes:
        if subattribute not in sks_settings[attribute]:
            msg = "The '{}' subattribute from the '{}' attribute was missing.".format(subattribute, attribute)
            logger.warning(msg)
            if subattribute == 'data':
                sks_settings[attribute][subattribute] = ''

    # If 'data' is not empty
    if sks_settings[attribute]['data'] != '':

        # If data type is valid:
        if isinstance(sks_settings[attribute]['data'], (str, list)):
            
            # If data are of type str:
            if isinstance(sks_settings[attribute]['data'], str):
                location = sks_settings[attribute]['data']
                
                # Checks if the data datafile exist
                if not os.path.exists(location):
                    msg = "The 'data' subattribute from the '{}' attribute is not valid. '{}' does not exist. The subattribute has been cleaned.".format(attribute, location)
                    logger.warning(msg)
                    sks_settings[attribute]['data'] = ''
                else:
                    try:
                        data = np.genfromtxt(location)
                    except Exception as err:
                        msg = "Impossible to read the file located by the 'data' subattribute from the '{}' attribute.".format(attribute)
                        logger.error(msg)
                        ERRORS.append(err)

            # If data are of type list:
            if isinstance(sks_settings[attribute]['data'], list):
                data = sks_settings[attribute]['data']

        # Data type is not recognized
        else:
            msg = "The 'data' subattribute type from the '{}' attribute is not valid. Type must be str or list. The subattribute has been cleaned.".format(attribute)
            logger.warning(msg)
            sks_settings[attribute]['data'] = ''
    

    # Checks if data is valid
    if data != '':

        # Checks if there is at least 3 data declared
        if len(data) < 3:
            msg = "Not enough vertices to create a mask delimitation (3 minimum). The 'data' subattribute has been cleaned."
            logger.warning(msg)
            sks_settings[attribute]['data'] = ''

        # Checks if each vertex contains 2 coordinates
        for vertex in data:
            if len(vertex) != 2:
                msg = "The 'data' subattribute from the '{}' attribute contains at least one invalid vertex. Format must be like : [[x0, y0], ..., [xn, yn]]. The subattribute has been cleaned.".format(attribute)
                logger.warning(msg)
                sks_settings[attribute]['data'] = ''
                break
            
        # Checks type of coordinates
        for vertex in data:
            for coordinate in vertex:
                if np.isnan(coordinate) or (not isinstance(coordinate, (int, float))):
                    msg = "The 'data' subattribute from the '{}' attribute contains at least one invalid vertex. Coordinates must be of type int or float. The subattribute has been cleaned.".format(attribute)
                    logger.warning(msg)
                    sks_settings[attribute]['data'] = ''
                    break

    ##################
    ### TOPOGRAPHY ###
    ##################
    attribute = 'topography'
    subattributes = ['data']

    # Checks presence of the 'topography' attribute
    if attribute not in sks_settings:
        msg = "The '{}' attribute was missing. Topography will be set from geological data.".format(attribute)
        logger.warning(msg)
        sks_settings[attribute] = {
            'data' : ''
        }

    # Checks presence of the subattributes of the 'topography' attribute
    for subattribute in subattributes:
        if subattribute not in sks_settings[attribute]:
            msg = "The '{}' subattribute from the '{}' attribute was missing.".format(subattribute, attribute)
            logger.warning(msg)
            if subattribute == 'data':
                sks_settings[attribute][subattribute] = ''

    # If 'data' is not empty
    if sks_settings[attribute]['data'] != '':

        # If data type is valid:
        if isinstance(sks_settings[attribute]['data'], str):
            location = sks_settings[attribute]['data']
                
            # Checks if the data datafile exist
            if not os.path.exists(location):
                msg = "The 'data' subattribute from the '{}' attribute is not valid. '{}' does not exist. The subattribute has been cleaned.".format(attribute, location)
                logger.warning(msg)
                sks_settings[attribute]['data'] = ''
            else:
                try:
                    data = np.genfromtxt(location, skip_header=3)
                except Exception as err:
                    msg = "Impossible to read the file located by the 'data' subattribute from the '{}' attribute.".format(attribute)
                    logger.error(msg)
                    ERRORS.append(err)

            

        # Data type is not recognized
        else:
            msg = "The 'data' subattribute type from the '{}' attribute is not valid. Type must be str. The subattribute has been cleaned.".format(attribute)
            logger.warning(msg)
            sks_settings[attribute]['data'] = ''


    #########################
    ### GEOLOGIC FEATURES ###
    #########################
    attributes = ['geology', 'karst', 'field']
    subattributes = ['data']
    valid_extensions = ['gslib', 'npy', 'png', 'jpg']
    nx, ny, nz = sks_settings['grid']['nx'], sks_settings['grid']['ny'], sks_settings['grid']['nz']

    # Checks presence of the geologic feature attributes
    for attribute in attributes:
        if attribute not in sks_settings:
            msg = "The '{}' attribute was missing.".format(attribute)
            logger.warning(msg)
            sks_settings[attribute] = {
                'data' : ''
            }

    # Checks presence of the subattributes
    for attribute in attributes:
        for subattribute in subattributes:
            if subattribute not in sks_settings[attribute]:
                msg = "The '{}' subattribute from the '{}' attribute was missing.".format(subattribute, attribute)
                logger.warning(msg)
                if subattribute == 'data':
                    sks_settings[attribute][subattribute] = ''


    # Checks if datafiles are valid
    for attribute in attributes:
        location = sks_settings[attribute]['data']
        if location != '':

            # Checks if datafile location type is valid:
            if not isinstance(location, str):
                msg = "The '{}' datafile location type is not valid. Type must be str. The attribute has been cleaned.".format(attribute) # TODO
                logger.warning(msg)
                sks_settings[attribute]['location'] = ''
            else:

                # Checks if location is valid
                if not os.path.exists(location):
                    msg = "The '{}' datafile location is not valid. '{}' does not exist. The attribute has been cleaned.".format(attribute, location)
                    logger.warning(msg)
                    sks_settings[attribute]['location'] = ''
                else:

                    # Checks if extension is valid
                    extension = location.split('.')[-1]
                    if extension not in valid_extensions:
                        msg = "The extension of the '{}' datafile location is not valid. '.{}' extension is not recognized. Valid extensions : {}. The attribute has been cleaned.".format(attribute, extension, valid_extensions)
                        logger.warning(msg)
                        sks_settings[attribute]['location'] = ''
                    else:

                        # Checks if datafile is readable and if size is valid
                        ### GSLIB
                        if extension == 'gslib':
                            try:
                                data = np.genfromtxt(location, skip_header=3, dtype=np.int8)
                            except Exception as err:
                                msg = "The '{}' datafile is not valid. The attribute has been cleaned.".format(attribute)
                                logger.error(msg)
                                ERRORS.append(err)
                                sks_settings[attribute]['location'] = ''
                            else:
                                # Checks if data dimensions are valid
                                data_size = len(data)
                                if data_size == grid_volume:
                                    pass
                                elif (data_size == grid_surface) and (grid_surface != grid_volume):
                                    msg = "The '{}' datafile dimensions match with grid surface dimensions. Data will be replicated on z-axis.".format(attribute)
                                    logger.warning(msg)
                                else:
                                    msg = "The '{}' datafile dimensions do not match with grid volume dimensions (data : {}, grid : {}). The attribute has been cleaned.".format(attribute, data_size, grid_volume)
                                    logger.warning(msg)
                                    sks_settings[attribute]['location'] = ''

                        ### Numpy_pickle
                        elif extension == 'npy':
                            try:
                                data = np.load(location)
                            except Exception as err:
                                msg = "The '{}' datafile is not valid. The attribute has been cleaned.".format(attribute)
                                logger.error(msg)
                                ERRORS.append(err)
                                sks_settings[attribute]['location'] = ''
                            else:
                                if data.shape != (nx, ny, nz):
                                    msg = "The '{}' datafile dimensions do not match with the grid dimensions (data : {}, grid : {}). The attribute has been cleaned.".format(attribute, data.shape, (nx, ny, nz))
                                    logger.warning(msg)
                                    sks_settings[attribute]['location'] = ''

    ###########
    ### END ###
    ###########

    print(sks_settings)

    # TODO
    # Count the warnings ?
    # Create subfunctions for each test

    if len(ERRORS) > 0:
        msg = "{} error(s) detected in the sks settings file".format(len(ERRORS))
        logger.error(msg)
        raise ERRORS[0]
    else:
        msg = "The sks settings file has been validated."
        logger.info(msg)

    return sks_settings


def validate_sim_settings(sim_settings:dict):
    """
    DOCSTRING

    Inputs
    ------

    Returns
    -------
    """
    global ERRORS
    logger = logging.getLogger("validations - simulations")

    
    ########################
    ### OUTLETS - INLETS ###
    ########################
    attributes = ['outlets', 'inlets']
    subattributes_o = ['data', 'number', 'shuffle', 'importance']
    subattributes_i = ['data', 'number', 'shuffle', 'per_outlet', 'importance']
    subattributes_  = [subattributes_o, subattributes_i]

    for attribute, subattributes in zip(attributes, subattributes_):

        # Checks presence of the attributes in the file
        if attribute not in sim_settings:
            msg = "The '{}' attribute is missing.".format(attribute)
            logger.error(msg)
            ERRORS.append(KeyError(msg))
        else:
            # Checks presence of the subattributes in the file
            for subattribute in subattributes:
                if subattribute not in sim_settings[attribute]:
                    msg = "The '{}' subattribute from the '{}' attribute is missing.".format(subattribute, attribute)
                    logger.error(msg)
                    ERRORS.append(KeyError(msg))

        if len(ERRORS) == 0:

            # Checks if 'number' is of type int
            if not isinstance(sim_settings[attribute]['number'], int):
                msg = "The 'number' subattribute from the '{}' attribute is not of type int.".format(attribute)
                logger.error(msg)
                ERRORS.append(TypeError(msg))

            # Checks if 'number' value is valid
            if sim_settings[attribute]['number'] < 1:
                msg = "The 'number' subattribute from the '{}' attribute cannot be lower than 1. At least one point must be declared.".format(attribute)
                logger.error(msg)
                ERRORS.append(ValueError(msg))

            # Checks data when provided
            if not (sim_settings[attribute]['data'] == '') or (sim_settings[attribute]['data'] == []):

                # Checks if data type is valid:
                if not isinstance(sim_settings[attribute]['data'], (str, list)):
                    msg = "The 'data' subbatribute type from the '{}' attribute is not valid. Type must be str or list.".format(attribute)
                    logger.error(msg)
                    ERRORS.append(KeyError(msg))
                else:

                    # If data is of type str:
                    if isinstance(sim_settings[attribute]['data'], str):
                        # Checks if data exist
                        location = sim_settings[attribute]['data']
                        if not os.path.exists(location):
                            msg = "The 'data' location from the '{}' attribute is not valid. '{}' does not exist.".format(attribute, location)
                            logger.error(msg)
                            ERRORS.append(KeyError(msg))
                        else:
                            try:
                                points = np.genfromtxt(location)
                            except Exception as err:
                                msg = "The 'data' datafile from the '{}' attribute is not valid.".format(attribute)
                                logger.error(msg)
                                ERRORS.append(err)

                    # If data is of type list:
                    if isinstance(sim_settings[attribute]['data'], list):
                        points = sim_settings[attribute]['data']

                    # Checks if each points contains 2 coordinates
                    for point in points:
                        if len(point) != 2:
                            msg = "The 'data' subattribute from the '{}' attribute contains at least one invalid point. Format must be like : [[x0, y0], ..., [xn, yn]].".format(attribute)
                            logger.error(msg)
                            ERRORS.append(ValueError(msg))
                            break
                    
                    # Checks type of coordinates
                    for point in points:
                        for coordinate in point:
                            if np.isnan(coordinate) or (not isinstance(coordinate, (int, float))):
                                msg = "The 'data' subattribute from the '{}' attribute contains at least one invalid point. Coordinates must be of type int or float.".format(attribute)
                                logger.error(msg)
                                ERRORS.append(ValueError(msg))
                                break

    # TODO ###################################################
    attributes = ['outlets', 'inlets']
    if len(sim_settings['outlets']['importance']) == 0:
        pass
    
        # Proceed to assert on some variables
                # features = [self.settings['outlets_importance'], self.settings['inlets_importance'], self.settings['inlets_per_outlet']]
                # for feature in features:
                #     assert isinstance(feature, list)
                # assert len(self.settings['outlets_importance']) == len(self.outlets)
    ###################################################


    ###############
    ### TRACERS ###
    ###############

    pass

    ###########
    ### END ###
    ###########

    print(sim_settings)

    if len(ERRORS) > 0:
        msg = "{} error(s) detected in the simulations settings file".format(len(ERRORS))
        logger.error(msg)
        raise ERRORS[0]
    else:
        msg = "The sks simulations file has been validated."
        logger.info(msg)
    return sim_settings