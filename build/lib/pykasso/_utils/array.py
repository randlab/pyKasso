"""
Module defining array manipulation functions.
"""

### External dependencies
import numpy as np


def normalize_array(array: np.ndarray) -> np.ndarray:
    """
    The function takes a numpy array and returns a new numpy array where each
    element is scaled to have values between 0 and 1.

    Parameters
    ----------
    array : np.ndarray
        Numpy array to normalize.

    Returns
    -------
    np.ndarray
        Normalized numpy array, with values ranging from 0 to 1.
    """
    min_value = array.min()
    max_value = array.max()
    normalized = (array - min_value) / (max_value - min_value)
    return normalized
