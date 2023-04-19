"""
This module contains ...
"""

### Internal dependencies
import sys
import importlib

### External dependencies
import numpy as np

### Module variables
this = sys.modules[__name__]


def visualizer(project_directory: str, notebook: bool = True, *args, **kwargs):
    """
    TODO
    """
    # Tests if pyvista is installed
    find_pyvista = importlib.util.find_spec("pyvista")
    this.is_pyvista_installed = (find_pyvista is not None)
    if this.is_pyvista_installed:
        from ._pyvista import PyvistaVisualizer
        visualizer = __class_factory(PyvistaVisualizer)(project_directory,
                                                        notebook,
                                                        *args, **kwargs)
    else:
        from ._pyplot import PyplotVisualizer
        visualizer = __class_factory(PyplotVisualizer)(project_directory,
                                                       notebook,
                                                       *args, **kwargs)
    return visualizer


def __class_factory(BaseClass):
    class SpecificClass(BaseClass):
        def __init__(self, *args, **kwargs):
            super(SpecificClass, self).__init__(*args, **kwargs)
    return SpecificClass


def show_array(array: np.ndarray, ghost=False) -> None:
    """
    TODO
    """
    if this.is_pyvista_installed:
        from ._pyvista import _show_array
    else:
        from ._pyplot import _show_array
    _show_array(array, ghost)
    return None


# def show(environment, feature, settings={}):
#     """
#     TODO
#     """

#     ### Controls validity of settings
#     attributes = ['domain', 'inlets', 'outlets', 'tracers']
#     for attribute in attributes:
#         if attribute in settings:
#             if not hasattr(environment, attribute):
#                 print('WARNING : No {} available'.format(attribute))
#                 del settings[attribute]
#             else:
#                 if getattr(environment, attribute) is None:
#                     print('WARNING : No {} available'.format(attribute))
#                     del settings[attribute]

#     ### Plot feature
#     if _is_feature_valid(environment, feature):
#         _show_data(environment, feature, settings)

#     return None


    
# ###################
# ### DEBUG PLOTS ###
# ###################

# def debug(environment, step, engine='matplotlib', settings={}):
#     """
#     TODO
#     """
#     if engine == 'pyvista':
#         from pykasso.visualization import _pyvista

#         if step == 'model':
#             _pyvista._debug_plot_model(environment, settings)
#         elif step == 'fmm':
#             _pyvista._debug_plot_fmm(environment, settings)

#     return None
