""" TODO"""

import sys
import logging

this = sys.modules['pykasso.core.sks']

def _decorator_logging(stage, feature_name):
    def _(function):
        # this.logger.critical(function.__name__)
        # print(function.__name__)
        def _wrapper_logging(*args, **kwargs):
            words = {
                'validation'   : [
                    "Validating the '{}' settings",
                    "'{}' settings have been validated",
                ],
                'construction' : [
                    "Constructing '{}'",
                    "'{}' has been constructed",
                ],
                'initialization' : [
                    "Initializating '{}'",
                    "'{}' has been initialized",
                ],
            }
            if feature_name in ['points', 'geology'] :
                feature = this.feature_kind
            else:
                feature = feature_name
            
            logger = logging.getLogger("{}.{}".format(feature, stage))
            logger.debug(words[stage][0].format(feature))
            try:
                result = function(*args, **kwargs)
            except Exception as err:
                logger.critical("Critical error during '{}' {}".format(feature, stage))
                raise err
            else:
                logger.debug(words[stage][1].format(feature))
                return result
                    
        return _wrapper_logging
    return _