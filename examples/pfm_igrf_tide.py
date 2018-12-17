# -*- coding: utf-8 -*-
"""
测试igrf和tide的调用
"""
#stdlib imports

#third party imports
# local imports
from datetime import datetime
import geoist.pfm as gip

m1=gip.IGRF()
aa=m1.pnt(45,105,100)
print(aa)
g1=gip.TideModel()
gdate = datetime(2018, 11, 1, 12, 00, 00)  
gm, gs, g =g1.solve_longman(45,105,100,gdate)  

value='2018-11-01 12:00:00'
gdd=datetime.strptime(value,'%Y-%m-%d %H:%M:%S')
gm, gs, g =g1.solve_longman(45,105,100,gdd)
print("Earthtide value in this time is:",g)