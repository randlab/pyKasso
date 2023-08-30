"""
TODO
"""

### Internal dependencies


### External dependencies
import numpy as np


def is_parameter_type_valid(parameter, parameter_name, types):
    """Check value type."""
    # print(parameter)
    # print(type(parameter))
    # print(types)
    if not isinstance(parameter, types):
        msg = ("Parameter '{}' must be of type(s) : {}"
               .format(parameter_name, types))
        raise TypeError(msg)
    else:
        return True


def is_key_in_dict(dictionary, dictionary_name, key):
    """Check key presence in dictionary."""
    if key not in dictionary:
        msg = ("Key '{}' is missing in '{}' dictionary."
               .format(key, dictionary_name))
        raise KeyError(msg)
    else:
        return True


def is_parameter_value_valid(parameter, parameter_name,
                             logical_test, compared_to):
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
