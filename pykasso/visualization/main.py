"""
This module contains ...
"""

# TODO - automatic update when visualization class is already created

### Internal dependencies
import sys
import importlib

### External dependencies
import yaml
import numpy as np

### Internal dependencies
from ._pyplot import PyplotVisualizer

### Module variables
this = sys.modules[__name__]

GEOLOGICAL_FEATURES = [
    'grid',
    'domain',
    'geology',
    # 'beddings'
    'faults',
    'fractures',
]

SURFACES_FEATURES = [
    'topography',
    'bedrock',
    'water_level',
]

DOMAIN_FEATURES = [
    'delimitation',
    'topography',
    'bedrock',
    'water_level',
]

ANISOTROPIC_FEATURES = [
    'cost',
    'alpha',
    'beta',
    'time',
    'karst',
]

features = [
    GEOLOGICAL_FEATURES,
    DOMAIN_FEATURES,
    ANISOTROPIC_FEATURES
]
AUTHORIZED_FEATURES = [f for list_ in features for f in list_]

# Tests if pyvista is installed
find_pyvista = importlib.util.find_spec("pyvista")
is_pyvista_installed = find_pyvista is not None
if is_pyvista_installed:
    from ._pyvista import PyvistaVisualizer
    subclass = PyvistaVisualizer
else:
    subclass = PyplotVisualizer


class Visualizer(subclass):
    """
    TODO
    """
    # ajouter __method__ pour dire qu'elle existe pas
    # si yapa pyvista d'installé
    
    def __init__(self, project_directory: str, notebook: bool = True,
                 *args, **kwargs) -> None:
        """"""
        super().__init__(project_directory, *args, **kwargs)
        
    def _is_feature_name_valid(self, feature) -> bool:
        """
        TODO
        """
        if feature not in AUTHORIZED_FEATURES:
            msg = ("ERROR : selected feature has been not recognized "
                   "(authorized features : {})".format(AUTHORIZED_FEATURES))
            raise ValueError(msg)
        else:
            return True
        
    def _is_feature_data_valid(self, simulation, feature) -> bool:
        """
        TODO
        """
        if feature in GEOLOGICAL_FEATURES:
            if getattr(simulation, feature) is None:
                msg = "ERROR : '{}' attribute is None type.".format(feature)
                raise ValueError(msg)
        
        elif feature in DOMAIN_FEATURES:
            if getattr(simulation, 'domain') is None:
                msg = "ERROR : '{}' attribute is None type.".format(feature)
                raise ValueError(msg)
            if getattr(simulation.domain, feature) is None:
                msg = ("ERROR : '{}' domain's attribute"
                       "is of None type.".format(feature))
                raise ValueError(msg)
    
        elif feature in ANISOTROPIC_FEATURES:
            if not hasattr(simulation, 'maps'):  # TODO - à supprimer ?
                msg = "ERROR : environment has no 'maps' attribute."
                raise ValueError(msg)
            
        if feature in ['alpha', 'beta']:
            if simulation.fmm['algorithm'] != 'Riemann3':
                msg = "ERROR : environment has no '{}' data.".format(feature)
                raise ValueError(msg)
            
        return True
    
    
def show_array(array: np.ndarray, engine: str = 'pyvista',
               ghost=False) -> None:
    """
    TODO
    """
    if is_pyvista_installed:
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
