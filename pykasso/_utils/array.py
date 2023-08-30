"""
TODO
"""

### External dependencies
import numpy as np


def normalize_array(array: np.ndarray) -> np.ndarray:
    min_value = array.min()
    max_value = array.max()
    normalized = (array - min_value) / (max_value - min_value)
    return normalized


def transform_array_where(array: np.ndarray, value: np.ndarray,
                          domain: np.ndarray) -> np.ndarray:
    test = (domain == 1)
    array = np.where(test, value, array)
    return array
