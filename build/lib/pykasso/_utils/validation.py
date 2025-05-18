"""
Module containing functions used for data flow validation.
"""

### Internal dependencies
import os
import logging

### External dependencies
import numpy as np

from ..core._namespaces import (
    SKS_VALID_ALGORITHM_VALUES,
    SKS_VALID_MODE_VALUES,
)

### Typing
from typing import Union


def is_filepath_valid(
    filepath: str,
) -> Union[FileNotFoundError, bool]:
    """
    Test if the filepath is valid.
    
    Parameters
    ----------
    filepath : str.
        String defining the filepath.
        
    Returns
    -------
    Union[FileNotFoundError, bool]
        Return ``True`` if test pass.
        Otherwise raise a ``FileNotFoundError`` exception.
        
    Raises
    ------
    FileNotFoundError
    """
    if not os.path.exists(filepath):
        msg = ("Filepath '{}' does not exist.".format(filepath))
        raise FileNotFoundError(msg)
    else:
        return True


def is_variable_type_valid(
    variable_name: str,
    variable_value,
    valid_types: tuple,
) -> Union[TypeError, bool]:
    """
    Test if the type of the variable is valid.

    Parameters
    ----------
    variable_name : str
        Name of the parameter.
    variable_value : any
        Value of the parameter.
    valid_types : tuple
        Accepted types.

    Returns
    -------
    Union[TypeError, bool]
        Return ``True`` if test pass.
        Otherwise raise a ``TypeError`` exception.

    Raises
    ------
    TypeError
    """
    if isinstance(variable_value, valid_types):
        return True
    else:
        msg = ("Parameter '{}' type is invalid. Valid type(s) : {}"
               .format(variable_name, valid_types))
        raise TypeError(msg)


def is_key_in_dict(
    dictionary: dict,
    dictionary_name: str,
    key_to_test: str,
) -> Union[KeyError, bool]:
    """
    Test key presence in the dictionary.

    Parameters
    ----------
    dictionary : dict
        Dictionary to test.
    dictionary_name : str
        Name of the dictionary.
    key_to_test : str
        Key to verify presence in dictionary.

    Returns
    -------
    Union[KeyError, bool]
        Return ``True`` if test pass.
        Otherwise raise a ``KeyError`` exception.

    Raises
    ------
    KeyError
    """
    if key_to_test in dictionary:
        return True
    else:
        msg = ("Key '{}' is missing in '{}' dictionary."
               .format(key_to_test, dictionary_name))
        raise KeyError(msg)


def is_variable_in_list(
    variable_name: str,
    variable_value,
    accepted_values: list,
) -> Union[ValueError, bool]:
    """
    TODO
    """
    if variable_value not in accepted_values:
        msg = ("Parameter '{}' value is invalid. Accepted values : {}"
               .format(variable_name, accepted_values))
        raise ValueError(msg)
    else:
        return True


def is_parameter_comparison_valid(
    parameter_name: str,
    parameter_value,
    logical_test: str,
    compared_to,
) -> Union[ValueError, bool]:
    """
    Test if the comparision returns true.

    Parameters
    ----------
    parameter_name : str
        Name of the parameter.
    parameter_value : any
        Value of the parameter.
    logical_test : str
        Logical test to use ('>', '>=', '<', '<=', '==', '!=').
    compared_to : _type_
        Value to compare.
        
    Returns
    -------
    Union[ValueError, bool]
        Return ``True`` if test pass.
        Otherwise raise a ``ValueError`` exception.

    Raises
    ------
    ValueError
    """
    logical_test_text = {
        '>': 'greater than',
        '>=': 'greater than or equal to',
        '<': 'less than',
        '<=': 'less than or equal to',
        '==': 'equal to',
        '!=': 'not equal to'
    }
    test = str(parameter_value) + logical_test + str(compared_to)
    if not eval(test):
        msg = ("The value of the '{}' parameter must be {} {}."
               .format(parameter_name,
                       logical_test_text[logical_test],
                       compared_to))
        raise ValueError(msg)
    else:
        return True

##########################
### Dictionary testing ###
##########################


def test_sks_settings(settings: dict) -> None:
    """
    """
    logger = logging.getLogger("sks.validation")
    
    ### Test 'seed' value
    try:
        is_variable_type_valid(
            variable_name='seed',
            variable_value=settings['seed'],
            valid_types=(int),
        )
    except TypeError as error:
        logger.error(error)
        raise
    
    ### Test 'algorithm' value
    try:
        is_variable_in_list(
            variable_name='algorithm',
            variable_value=settings['algorithm'],
            accepted_values=SKS_VALID_ALGORITHM_VALUES
        )
    except ValueError as error:
        logger.error(error)
        raise
    
    ### Test 'costs' value
    try:
        is_variable_type_valid(
            variable_name='costs',
            variable_value=settings['costs'],
            valid_types=(dict),
        )
    except TypeError as error:
        logger.error(error)
        raise
    
    ### Test 'factors' value
    try:
        is_variable_type_valid(
            variable_name='factors',
            variable_value=settings['factors'],
            valid_types=(dict),
        )
    except TypeError as error:
        logger.error(error)
        raise
    
    ### Test 'mode' value
    try:
        is_variable_in_list(
            variable_name='mode',
            variable_value=settings['mode'],
            accepted_values=SKS_VALID_MODE_VALUES
        )
    except ValueError as error:
        logger.error(error)
        raise
    
    return None


# def test_geologic_feature_settings(settings: dict) -> None:
#     """
#     """
#     return None

# def test_point_settings(kind: str, settings: dict) -> None:
#     """
#     """
#     logger = logging.getLogger("{}.validation".format(kind))
    
#     ### Test 'number' value
#     try:
#         is_variable_type_valid(
#             variable_name='number',
#             variable_value=settings['number'],
#             valid_types=(int),
#         )
#     except TypeError as error:
#         logger.error(error)
#         raise
    
#     ### Test 'data' value
    
    
#     ### Test 'shuffle' value
    
    
#     ### Test 'importance' value
    
    
#     ### Test 'subdomain' value
    
    
#     ### Test 'geology' value
    
    
#     ### Test 'seed' value
    
    
#     return None
