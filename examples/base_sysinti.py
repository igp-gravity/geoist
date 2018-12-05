# -*- coding: utf-8 -*-
"""
Created on Sat Dec  1 10:46:03 2018
测试导入模块后的init内容与作用
@author: chens
"""

import geoist.catqc as cat
print(cat.dirname(__file__))
print(cat.basename(__file__))
print(cat.isfile(__file__))
modules = cat.glob.glob(cat.dirname(__file__)+r"//*.py")
print(modules)
