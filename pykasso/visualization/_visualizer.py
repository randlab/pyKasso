"""
This module contains ...
"""

### Local dependencies
from .._utils import ProjectReader

### Module variables
GEOLOGICAL_FEATURES = [
    'grid',
    'domain',
    'geology',
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


class Visualizer(ProjectReader):
    """
    TODO
    """
    
    def __init__(self, project_directory: str, *args, **kwargs) -> None:
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
