"""
Module containing functions used for data flow validation.
"""

### Internal dependencies


### External dependencies
import numpy as np

### Typing
from typing import Union

# class InvalidTypeError(Exception):
#     pass

# class variable(object):
#     def __init__(self, value=0):
#         self.__x = value

#     def __set__(self, obj, value):
#         if value < 0:
#             raise InvalidTypeError('x is less than zero')

#         self.__x  = value

#     def __get__(self, obj, objType):
#         return self.__x

# class MyClass(object):
#     x = variable()


def is_parameter_type_valid(parameter_name: str,
                            parameter_value,
                            valid_types: tuple,
                            ) -> Union[TypeError, bool]:
    """
    Test if the type of a parameter is valid.

    Parameters
    ----------
    parameter_name : str
        Name of the parameter.
    parameter_value : any
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
    if isinstance(parameter_value, valid_types):
        return True
    else:
        msg = ("Parameter '{}' type is invalid. Valid type(s) : {}"
               .format(parameter_name, valid_types))
        raise TypeError(msg)


def is_key_in_dict(dictionary: dict,
                   dictionary_name: str,
                   key_to_test: str
                   ) -> Union[KeyError, bool]:
    """
    Test if key is present in dictionary.

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
        


def is_parameter_value_valid(parameter, parameter_name,
                             logical_test, compared_to) -> bool:
    """Check validity of the parameter value."""
    logical_test_text = {
        '>': 'greater than',
        '>=': 'greater than or equal to',
        '<': 'less than',
        '<=': 'less than or equal to',
        '==': 'equal to',
        '!=': 'not equal to'
    }
    test = str(parameter) + logical_test + str(compared_to)
    if not eval(test):
        msg = ("The value of the '{}' parameter must be {} {}."
               .format(parameter_name,
                       logical_test_text[logical_test],
                       compared_to))
        raise ValueError(msg)
    else:
        return True
