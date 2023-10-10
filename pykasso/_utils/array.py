"""
TODO
"""

### External dependencies
import numpy as np


def normalize_array(array: np.ndarray) -> np.ndarray:
    """TODO"""
    min_value = array.min()
    max_value = array.max()
    normalized = (array - min_value) / (max_value - min_value)
    return normalized
