"""pyKasso's wrappers functions."""

### Internal dependencies
import sys
import time
import logging

### Local dependencies
import pykasso.model._validations as val

### Modules variables
# TODO - save time exec somewhere and remove the 'this' var
this = sys.modules[__name__]
this.time_functions = {}


# def _debug_level(level: float, iteration_mode: bool = False):
#     """Debugging system.

#     Parameters
#     ----------
#     level : float
#         Debug level.
#     iteration_mode : bool, optional
#         _description_, by default False
#     """
#     def _(function):
#         def _wrapper(*args, **kwargs):
#             model = args[0] # retrieves the model
#             result = ''
#             if not model.debug_mode:
#                 result = function(*args, **kwargs)
#             else:
#                 if not iteration_mode:
#                     if model.debug_level['model'] >= level:
#                         result = function(*args, **kwargs)
#                 else:
#                     if (model.debug_level['simulation'] >= level) and (model.iteration <= model.debug_level['iteration']):
#                         result = function(*args, **kwargs)
#             if result == '':
#                 msg = "DEBUG MODE : {}".format(level)  # TODO
#                 sks.logger.critical(msg)
#                 sys.exit(msg)
#             else:
#                 return result
#         return _wrapper
#     return _


def _parameters_validation(feature, kind):
    def _(function):
        def _wrapper(*args, **kwargs):
            logger = logging.getLogger("{}.validation".format(feature))
            model = args[0]
            model.model_parameters = val.validate_attribute_presence(
                model.model_parameters, feature, kind, default_value={})
            if feature in ['sks', 'grid']:
                model.model_parameters[feature] = val.validate_settings(
                    feature, model.model_parameters[feature])
            else:
                model.model_parameters[feature] = val.validate_settings(
                    feature, model.model_parameters[feature], model.grid)
            msg = "'{}' settings have been validated".format(feature)
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
            logger = logging.getLogger("{}.construction".format(feature))
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
                t0 = time.perf_counter()
                result = function(*args, **kwargs)
                t1 = time.perf_counter()
            except Exception as err:
                msg = "Critical error during '{}'".format(function.__name__)
                logger.critical(msg)
                raise err
            else:
                logger.debug("'{}' went well".format(function.__name__))
                this.time_functions[function.__name__] = t1 - t0
                return result
        return _wrapper
    return _
