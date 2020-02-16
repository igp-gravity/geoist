# -*- coding: utf-8 -*-
"""
Created on Sat Sep 15 18:00:46 2018

@author: chens
"""

import logging
import pathlib
import geoist

def setlogger(modname=__name__, level = logging.INFO, logfile = 'geoist.log'):
    logger = logging.getLogger(modname)
    logger.setLevel(level=level)
    filepath = geoist.USER_DATA_PATH
    #filepath = pathlib.Path(__file__).parent
    handler = logging.FileHandler(pathlib.Path(filepath,logfile))
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logger

logger = setlogger()

def setlevel(level = 20):
    logger.setLevel(level=level)

def setname(logname = __name__):
    logger.name = logname

def debug(str='Debugging'):  #10
    logger.debug(str)
      
def info(str='INFO'):   #20
    logger.info(str)

def warning(str='Warning'):  #30
    logger.warning(str)
    
def error(str='ERROR'):  #40
    logger.error(str)

def critical(str='CRITICAL'):  #50
    logger.critical(str)
