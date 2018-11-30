# -*- coding: utf-8 -*-
"""
 Name        : tri.py
 Created on  : 2018/11/03 17:00
 Author      : Steve Chen <chenshi@cea-igp.ac.cn>
 Affiliation : Institute of Geophysics, CEA.
 Version     : 0.1.0
 Copyright   : Copyright (C) 2018-2020 GEOIST Development Team. All Rights Reserved.
 License     : Distributed under the MIT License. See LICENSE.txt for more info.
 Github      : https://igp-gravity.github.io/
 Description : Application for ***.
"""

"""
Code function

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