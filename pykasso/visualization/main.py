"""
TODO
"""

from pykasso.visualization import _matplotlib
from pykasso.visualization import _pyvista

def show(environment, feature='grid', engine='matplotlib', settings={}):
    """
    TODO
    """
    # According to engine, selects appropriates functions
    if engine == 'matplotlib':
        _show_data = _matplotlib._show_data
    elif engine == 'pyvista':
        _show_data  = _pyvista._show_data
    else:
        pass


    if 'mask' in settings:
        if getattr(environment, 'MASK') is None:
            print('Warning : No mask available.')
            del settings['mask']


    authorized_features = [
        'grid',
        'geology',
        'fractures',
        'topography'
    ]
    if feature in authorized_features:
        _show_feature(environment, feature, _show_data, settings)
    elif feature == 'array':
        pass
    else:
        print('ERROR : selected feature has been not recognized. Authorized features : {}'.format(authorized_features))


    return None


def _show_feature(environment, feature, func, settings):
    """
    TODO
    """
    feature = feature.upper()
    if not hasattr(environment, feature):
        raise KeyError("ERROR : environment has no '{}' attribute.".format(feature))
    else:
        if getattr(environment, feature) is None:
            raise KeyError("ERROR : '{}' attribute is type of None.".format(feature))
        else:
            func(environment, feature, settings)
    return None