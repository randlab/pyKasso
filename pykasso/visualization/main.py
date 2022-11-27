"""
TODO
"""

from pykasso.visualization import _matplotlib
from pykasso.visualization import _pyvista

def show(environment, feature='grid', engine='matplotlib', settings=None):
    """
    TODO
    """
    # According to engine, selects appropriates functions
    if engine == 'matplotlib':
        _show_grid = _matplotlib._show_grid
        _show_mask = _matplotlib._show_mask
        _show_geology = _matplotlib._show_geology
    elif engine == 'pyvista':
        _show_grid        = _pyvista._show_grid
        _show_mask        = _pyvista._show_mask
        _show_geology     = _pyvista._show_geology
        _show_numpy_array = _pyvista._show_numpy_array
        _show_topography  = _pyvista._show_topography
    else:
        pass

    
    # Displays the selected feature
    if feature == 'grid':
        _show_feature(environment, feature, _show_grid, settings)

    elif feature == 'mask':
        _show_feature(environment, feature, _show_mask, settings)

    elif feature == 'geology':
        _show_feature(environment, feature, _show_geology, settings)

    elif feature == 'topography':
        _show_feature(environment, feature, _show_topography, settings)

    # elif feature == 'fractures':
    #     _show_feature(environment, feature, _show_fractures, settings)

    elif feature == 'numpy_array':
        # if not isinstance TODO
        _show_numpy_array(environment, settings)

    else:
        print('ERROR : selected feature has been not recognized.')

    return None


def _show_feature(environment, feature, func, settings):
    """
    TODO
    """
    feature = feature.upper()
    if not hasattr(environment, feature):
        raise KeyError("ERROR : environment has no '{}' attribute.".format(feature))
    else:
        func(environment, settings)
    return None