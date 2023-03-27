"""pyKasso's wrappers functions."""

### Internal dependencies
import sys
import time
import logging

### Local dependencies
import pykasso.core._validations as val

### Modules variables
sks  = sys.modules['pykasso.core.sks']
this = sys.modules[__name__]
this.time_functions = {}


def _debug_level(level:float, iteration_mode:bool=False):
    """Debugging system.

    Parameters
    ----------
    level : float
        Debug level.
    iteration_mode : bool, optional
        _description_, by default False
    """
    def _(function):
        def _wrapper(*args, **kwargs):
            model = args[0] # retrieves the model
            result = ''
            if not model.debug_mode:
                result = function(*args, **kwargs)
            else:
                if not iteration_mode:
                    if model.debug_level['model'] >= level:
                        result = function(*args, **kwargs)
                else:
                    if (model.debug_level['simulation'] >= level) and (model.iteration <= model.debug_level['iteration']):
                        result = function(*args, **kwargs)
            if result == '':
                msg = "DEBUG MODE".format(level) # TODO
                sks.logger.critical(msg)
                sys.exit(msg)
            else:
                return result
        return _wrapper
    return _

def _parameters_validation(feature, kind):
    def _(function):
        def _wrapper(*args, **kwargs):
            sks.logger = logging.getLogger("{}.validation".format(feature))
            model = args[0]
            model.SKS_SETTINGS = val.validate_attribute_presence(model.SKS_SETTINGS, feature, kind, default_value={})
            if feature in ['sks', 'grid']:
                model.SKS_SETTINGS[feature] = val.validate_settings(feature, model.SKS_SETTINGS[feature])
            else:
                model.SKS_SETTINGS[feature] = val.validate_settings(feature, model.SKS_SETTINGS[feature], model.grid)
            sks.logger.info("'{}' settings have been validated".format(feature))
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
            sks.logger = logging.getLogger("{}.construction".format(feature))
            model = args[0]
            if model.SKS_SETTINGS[feature] is not sks.ACTIVE_PROJECT['settings'][feature]:
                # print('is not') # TODO
                result = function(*args, **kwargs)
                sks.ACTIVE_PROJECT['settings'][feature] = model.SKS_SETTINGS[feature]
                sks.ACTIVE_PROJECT['model'][feature]    = getattr(model, feature)
                sks.logger.info("'{}' has been constructed".format(feature))
            else:
                setattr(model, feature, sks.ACTIVE_PROJECT['model'][feature])
                sks.logger.info("'{}' has been reused from previous simulation".format(feature))
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
                sks.logger = logging.getLogger("{}.{}".format(feature, step))
            try:
                t0 = time.perf_counter()
                result = function(*args, **kwargs)
                t1 = time.perf_counter()
            except Exception as err:
                sks.logger.critical("Critical error during '{}'".format(function.__name__))
                raise err
            else:
                sks.logger.debug("'{}' went well".format(function.__name__))
                this.time_functions[function.__name__] = t1 - t0
                return result
        return _wrapper
    return _