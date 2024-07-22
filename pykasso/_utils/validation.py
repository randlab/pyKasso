"""
Module containing functions used for data flow validation.
"""

### Internal dependencies
import os
import logging

### External dependencies
import numpy as np

### Typing
from typing import Union

# def is_list_length_valid(list_to_test: list, value: int, attribute: str) -> bool:
#     if len(data) < value:
#         msg = ("'{}' data length is too short ({} elements minimum)."
#                .format(attribute, value))
#         this.logger.critical(msg)
#         raise ValueError(msg)
#     else:
#         return True

def is_filepath_valid(filepath: str,
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
        Return ``true`` if test pass.
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


def is_variable_type_valid(variable_name: str,
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
        Return ``true`` if test pass.
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


def is_key_in_dict(dictionary: dict,
                   dictionary_name: str,
                   key_to_test: str
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
        Return ``true`` if test pass.
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


# def is_element_in_list(element_to_test: str,
#                        list_to_test: list,
#                        ) -> Union[KeyError, bool]:
#     """
#     Test element presence in the list.

#     Parameters
#     ----------
#     element_to_test : str
#         Element to verify presence in list.
#     list_to_test : list
#         List to test.

#     Returns
#     -------
#     Union[KeyError, bool]
#         Return ``true`` if test pass.
#         Otherwise raise a ``KeyError`` exception.

#     Raises
#     ------
#     KeyError
#     """
#     if element_to_test in list_to_test:
#         return True
#     else:
#         msg = ("Element '{}' is invalid. Valid element(s): '{}'."
#                .format(element_to_test, list_name))
#         raise KeyError(msg)
        

def is_parameter_comparison_valid(parameter_name: str,
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
        Return ``true`` if test pass.
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

# try:
#     do_something_in_app_that_breaks_easily()
# except AppError as error:
#     logger.error(error)
#     raise

