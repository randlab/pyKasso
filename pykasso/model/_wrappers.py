"""pyKasso's wrappers functions."""

### Internal dependencies
import sys
import time
import logging

### Local dependencies
import pykasso.model._validations as val
from pykasso.core._namespaces import DEFAULT_FMM_COSTS


DEFAULT_VALUES = {
    'sks': {
        'seed': 0,
        'algorithm': 'Isotropic3',
        'costs': DEFAULT_FMM_COSTS,
        'mode': 'A',
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
        'costs': {},
        'model': {},
        'seed': 0,
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
        # 'importance': ['required', []],
        'subdomain': 'domain_surface',
        'geology': None,
        'seed': 0,
    },
    'inlets': {
        # 'number': ['required', ''],
        'data': [],
        'shuffle': False,
        # 'importance': ['required', []],
        # 'per_outlet': ['required', []],
        'subdomain': 'domain_surface',
        'geology': None,
        'seed': 0,
    },
}


def _parameters_validation(feature, kind):
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
                costs = default_params['costs'].copy()
                default_costs = DEFAULT_FMM_COSTS.copy()
                default_costs.update(costs)
                default_params['costs'] = default_costs
            if feature == 'outlets':
                for key in ['number', 'importance']:
                    if key not in default_params:
                        msg = ("The mandatory '{}' attribute is missing."
                               ).format(key)
                        logger.error(msg)
                        raise KeyError(msg)
            if feature == 'inlets':
                for key in ['number', 'importance', 'per_outlet']:
                    if key not in default_params:
                        msg = ("The mandatory '{}' attribute is missing."
                               ).format(key)
                        logger.error(msg)
                        raise KeyError(msg)
            
            # Control values
            # TODO : finish the validating functions
            if feature == 'sks':
                val.validate_settings_sks(default_params)
            elif feature in ['geology', 'faults', 'fractures']:
                pass
            elif feature == 'domain':
                # val.validate_settings_domain(default_params)
                pass
            elif feature in ['inlets', 'outlets']:
                val.validate_settings_points(default_params, feature)
                # if isinstance(default_params['data'], str)
                
                # val.validate_settings_io(default_params)
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
    TODO
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
    TODO
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
