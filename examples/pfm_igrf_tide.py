# -*- coding: utf-8 -*-
"""
测试igrf和tide的调用
"""
#stdlib imports

#third party imports
# local imports
import datetime
import geoist.pfm as gip

m1=gip.IGRF()
m1.pnt(45,105,100)

g1=gip.TideModel()
gdate = datetime.datetime(2018, 11, 1, 12, 00, 00)    
gm, gs, g =g1.solve_longman(45,105,100,gdate)
print("Earthtide value in this time is:",g)