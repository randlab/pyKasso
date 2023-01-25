""" TODO"""

import sys
import logging
import pykasso.core._validations as val

this = sys.modules['pykasso.core.sks']

def _debug_level(level, iteration_mode=False):
    """
    TODO
    """
    def _(function):
        def _wrapper(*args, **kwargs):
            model = args[0]
            result = ''
            if not model.debug_mode:
                result = function(*args, **kwargs)
            else:
                if not iteration_mode:
                    if model.debug_level['model'] >= level:
                        result = function(*args, **kwargs)
                else:
                    if (model.debug_level['model'] >= level) and (model.iteration <= model.debug_level['iteration']):
                        result = function(*args, **kwargs)
            if result == '':
                msg = "DEBUG MODE".format(level) # TODO
                this.logger.critical(msg)
                sys.exit(msg)
            else:
                return result
        return _wrapper
    return _

def _parameters_validation(feature, kind):
    def _(function):
        def _wrapper(*args, **kwargs):
            this.logger = logging.getLogger("{}.validation".format(feature))
            model = args[0]
            model.SKS_SETTINGS = val.validate_attribute_presence(model.SKS_SETTINGS, feature, kind)
            if feature == 'grid':
                model.SKS_SETTINGS[feature] = val.validate_settings(feature, model.SKS_SETTINGS[feature])
            elif feature == 'sks':
                model.SKS_SETTINGS['sks'] = val.validate_settings(feature, model.SKS_SETTINGS[feature])
            elif feature == 'fmm':
                model.SKS_SETTINGS = val.validate_settings(feature, model.SKS_SETTINGS)
            else:
                model.SKS_SETTINGS[feature] = val.validate_settings(feature, model.SKS_SETTINGS[feature], model.grid)
            this.logger.info("'{}' settings have been validated".format(feature))
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
            this.logger = logging.getLogger("{}.construction".format(feature))
            model = args[0]
            if model.SKS_SETTINGS[feature] != this.ACTIVE_PROJECT['settings'][feature]:
                result = function(*args, **kwargs)
                this.ACTIVE_PROJECT['settings'][feature] = model.SKS_SETTINGS[feature]
                this.ACTIVE_PROJECT['model'][feature]    = getattr(model, feature)
                this.logger.info("'{}' has been constructed".format(feature))
            else:
                setattr(model, feature, this.ACTIVE_PROJECT['model'][feature])
                this.logger.info("'{}' has been reused from previous simulation".format(feature))
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
                this.logger = logging.getLogger("{}.{}".format(feature, step))
            try:
                result = function(*args, **kwargs)
            except Exception as err:
                this.logger.critical("Critical error during '{}'".format(function.__name__))
                raise err
            else:
                this.logger.debug("'{}' went well".format(function.__name__))
                return result
        return _wrapper
    return _