# -*- coding: utf-8 -*-
"""
Created on Sat Dec 22 10:38:53 2018

@author: chens
"""

import geoist as gi
import importlib

importlib.reload(gi)
gi.log.info('tttt')

#gi.log.setname('base_log')
#gi.log.setlevel(10)
gi.log.info('info tt')
gi.log.debug('debug')