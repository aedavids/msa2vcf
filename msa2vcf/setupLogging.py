'''
Created on Oct 30, 2018

@author: Andrew Davidson aedavids@ucsc.edu
'''

import os
import json
import logging.config

# # https://stackoverflow.com/a/20885799
# # https://docs.python.org/3.8/library/importlib.html?highlight=importlib#module-importlib.resources
# try:
#     import importlib.resources as pkg_resources
# except ImportError:
#     # Try backported to PY<37 `importlib_resources`.
#     import importlib_resources as pkg_resources

# https://stackoverflow.com/a/47972610
import pkg_resources


#from . import templates  # relative-import the *package* containing the templates


_alreadyInitialized = False
_default_path = 'logging.ini.json'  #'logging.json',
def setupLogging (
        default_path=_default_path,
        default_level=logging.WARN,
        env_key='LOG_CFG'):
    """
    configure python logging facilities. Be default uses config file logging.ini.json
    
    returns path to config file
    """
    path = default_path
    value = os.getenv(env_key, None)
    if value:
        path = value
      
    ret = None  
    #print("AEDWIP path:{}".format(path))
    if os.path.exists(path):
#         logger = logging.getLogger('root')
#         print("AEDWIP loading log configuration:{}".format(default_path))
#         with open(path, 'rt') as f:
#             config = json.load(f)
#         logging.config.dictConfig(config)
        ret = path
        _loadConfig(path)
        
    elif path == _default_path:
        # load config from package
#         package = pkg_resources.Package
#         print("AEDWP package {}".format(package))
#         foo = pkg_resources.open_text(package, _default_path)
#         print("foo:{}".format(foo))
        filepath = pkg_resources.resource_filename(__name__, path)
#         print("AEDWIP filepath:{}".format(filepath))
        ret = filepath
        _loadConfig(filepath)
    
    else:
        ret = "logging.basisConfig()"
        logging.basicConfig(level=default_level)
        #print("AEDWIP loading log loading basicConfig:")
        
    return ret        

def _loadConfig(path):
    logger = logging.getLogger('root')
#     print("AEDWIP loading log configuration:{}".format(path))
    with open(path, 'rt') as f:
        config = json.load(f)
    logging.config.dictConfig(config)
    #logger.warning("loading log configuration:{}".format(path))

