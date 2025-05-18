"""pyKasso's wrappers functions."""

### Internal dependencies
import logging

### Local dependencies
from .._utils.validation import test_sks_settings
from ..model import _validations as val
from ..core._namespaces import DEFAULT_FMM_COSTS


DEFAULT_VALUES = {
    'sks': {
        'seed': 0,
        'algorithm': 'Isotropic3',
        'costs': DEFAULT_FMM_COSTS,
        'factors': {'F': 100, 'F1': 100, 'F2': 50}
    },
    'geology': {
        'data': None,
        'axis': 'z',
        'names': {},
        'costs': {},
        'model': {}
    },
    'faults': {
        'data': None,
        'axis': 'z',
        'names': {},
        'costs': {},
        'model': {}
    },
    'fractures': {
        'data': None,
        'axis': 'z',
        'names': {},
        # 'costs': {},
        'model': {},
        'seed': None,
    },
    'domain': {
        'delimitation': None,
        'topography': None,
        'bedrock': None,
        'water_table': None,
    },
    'outlets': {
        # 'number': ['required', ''],
        'data': [],
        'shuffle': False,
        'importance': [1],
        'subdomain': 'domain_surface',
        'geology': None,
        'seed': None,
    },
    'inlets': {
        # 'number': ['required', ''],
        'data': [],
        'shuffle': False,
        'importance': [1],
        # 'per_outlet': [1],
        'subdomain': 'domain_surface',
        'geology': None,
        'seed': None,
    },
}


def _parameters_validation(feature, kind):
    """
    This decorator validates input parameters before creatings modeling classes.
    """
    def _(function):
        def _wrapper(*args, **kwargs):
            logger = logging.getLogger("validation.{}".format(feature))
            model = args[0]
            
            # Add feature dictionary if value is missing
            if feature not in model.model_parameters:
                if kind == 'required':
                    msg = "The '{}' key is missing.".format(feature)
                    logger.error(msg)
                    raise KeyError(msg)
                else:
                    model.model_parameters[feature] = {}
            
            # Add default feature values
            user_params = model.model_parameters[feature].copy()
            default_params = DEFAULT_VALUES[feature].copy()
            for (key, value) in default_params.items():
                if key not in user_params:
                    msg = ("The '{}' attribute is missing. Set to default"
                           " value.").format(key)
                    logger.warning(msg)
            default_params.update(user_params)
            
            # Test special key presences
            if feature == 'sks':
                # Travel cost
                costs = default_params['costs'].copy()
                default_costs = DEFAULT_FMM_COSTS.copy()
                default_costs.update(costs)
                default_params['costs'] = default_costs
                # Mode
                if default_params['algorithm'] == 'Isotropic3':
                    default_params['mode'] = 'A'
                else:
                    default_params.setdefault('mode', 'D')
            if feature in ['outlets', 'inlets']:
                for key in ['number']:
                    if key not in default_params:
                        msg = ("The mandatory '{}' attribute is missing."
                               ).format(key)
                        logger.error(msg)
                        raise KeyError(msg)
            
            # Control values
            if feature == 'sks':
                test_sks_settings(default_params)
            elif feature in ['geology', 'faults', 'fractures']:
                pass
            elif feature == 'domain':
                pass
            elif feature in ['inlets', 'outlets']:
                val.validate_settings_points(default_params, feature)
                # if isinstance(default_params['data'], str)
                
                pass
            
            # Update dictionary
            model.model_parameters[feature] = default_params
            msg = "'{}' parameters have been validated.".format(feature)
            logger.info(msg)
            result = function(*args, **kwargs)
            return model
        return _wrapper
    return _


def _memoize(feature):
    """
    This decorator caches the results of function calls, preventing the need
    to recompute results for the same inputs.
    """
    def _(function):
        def _wrapper(*args, **kwargs):
            logger = logging.getLogger("construction.{}".format(feature))
            model = args[0]
            memoization = model.project._memoization
            if model.model_parameters[feature] is not memoization['settings'][feature]:
                # print('is not') # TODO
                result = function(*args, **kwargs)
                memoization['settings'][feature] = model.model_parameters[feature]
                memoization['model'][feature] = getattr(model, feature)
                logger.info("'{}' has been constructed".format(feature))
            else:
                setattr(model, feature, memoization['model'][feature])
                msg = "'{}' has been reused from previous simulation".format(feature)
                logger.info(msg)
            return model
        return _wrapper
    return _


def _logging(feature=None, step=None):
    """
    This decorator records messages to a log fileand tracks events, errors,
    and informational messages.
    """
    def _(function):
        def _wrapper(*args, **kwargs):
            if feature is not None:
                logger = logging.getLogger("{}.{}".format(feature, step))
            else:
                logger = logging.getLogger(".")
            try:
                result = function(*args, **kwargs)
            except Exception as err:
                msg = "Critical error during '{}'".format(function.__name__)
                logger.critical(msg)
                raise err
            else:
                logger.debug("'{}' went well".format(function.__name__))
                return result
        return _wrapper
    return _
