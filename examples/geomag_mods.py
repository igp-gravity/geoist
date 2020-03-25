#!/usr/bin/env python
# coding: utf-8

# #### MAGNETIC FIELD MODEL COLLETCTION
# * Auther: Shi Chen at CEA-IGP
# * Date: 2019-09-29

import datetime as dt

from geoist.magmod.magnetic_model.loader_igrf import load_model_igrf
from geoist.magmod.magnetic_model.loader_wmm import load_model_wmm
from geoist.magmod.magnetic_model.loader_emm import load_model_emm
from geoist.magmod.magnetic_model.loader_shc import (
    load_model_shc, load_model_shc_combined,
)
from geoist.magmod.magnetic_model.loader_mio import (
    load_model_swarm_mio_internal,
    load_model_swarm_mio_external,
)

from geoist.magmod.data import (
    EMM_2010_STATIC, EMM_2010_SECVAR, WMM_2015, WMM_2020, IGRF13, 
    CHAOS6_CORE_LATEST, CHAOS6_STATIC,
    IGRF11, IGRF12, SIFM,
)
from geoist.magmod.time_util import (
    decimal_year_to_mjd2000, decimal_year_to_mjd2000_simple,mjd2000_to_decimal_year,mjd2000_to_year_fraction
)
from geoist.magmod.util import datetime_to_decimal_year, vnorm

from geoist.magmod.magnetic_model.parser_mio import parse_swarm_mio_file
from geoist.magmod.magnetic_model.tests.data import SWARM_MIO_SHA_2_TEST_DATA
import geoist.magmod._pymm as pymm

print(SWARM_MIO_SHA_2_TEST_DATA) #DIFI4 is a type of MIO SHA model
print(mjd2000_to_decimal_year(7305))
print(mjd2000_to_decimal_year([5479., 7305., 6392.0]))
d1 = dt.datetime(2015,1,1)   # import time , location(lat, lon)
d11 = datetime_to_decimal_year(d1) # datetime to decimal year
loc = (30.0, 40.0, 1000.0)
wmm2015 = load_model_wmm(WMM_2015)  #load wmm2015 model
igrf11 = load_model_igrf(IGRF11)    #load igrf11 model
igrf12 = load_model_shc(IGRF12, interpolate_in_decimal_years=True)    #load igrf12 model

igrf13 = load_model_igrf(IGRF13)
wmm2020 = load_model_wmm(WMM_2020)

emm = load_model_emm(EMM_2010_STATIC, EMM_2010_SECVAR)  #load emm model
options = {"scale": [1, 1, -1]}   #-1 is Z direction

## renew for the IGRF13 and WMM2020 models on 2020-03-22
d2 = 2022.5
loc1 = (80.0, 0.0, 100.0)
wmm2020.eval(decimal_year_to_mjd2000(d2), loc1, 0, 0, **options) 
igrf13.eval(decimal_year_to_mjd2000(d2), loc1, 0, 0, **options)
## the result has been checked with WMM2020testvalues.pdf 

wmm2015.eval(decimal_year_to_mjd2000(d11), loc, 0, 0, **options) # 0,0 mean input,output using GEODETIC_ABOVE_WGS84
igrf11.eval(decimal_year_to_mjd2000(d11), loc, 0, 0, **options)
igrf12.eval(decimal_year_to_mjd2000(d11), loc, 0, 0, **options) #North-X, East-Y, Vertical-Z
vnorm(igrf12.eval(decimal_year_to_mjd2000(d11), (30.0, 40.0, 1000.0), 0, 0)) #Total intensity
mjd2000_to_decimal_year([5479., 7305., 6392.0])

emm.eval(decimal_year_to_mjd2000(d11), loc, 0, 0, **options)

# load mio model
mio = load_model_swarm_mio_internal(SWARM_MIO_SHA_2_TEST_DATA)
options = {"f107": 70, "scale": [1, 1, -1]}
mio.eval(decimal_year_to_mjd2000(d11), (30.0, 40.0, 1000.0), 0, 0, **options)

mio1 = load_model_swarm_mio_external(SWARM_MIO_SHA_2_TEST_DATA)
mio1.eval(5661.87, [(30.0, 40.0, 6400.0), (30.0, 40.0, 8000.0),], **options)
print('mio1 has finished')
# load DIFI4 model
#import sys

path = "D:\\MyProjects\\2018-科学院地球所\\models\\"

difi3 = load_model_swarm_mio_internal(path+'SW_OPER_MIO_SHAi2D_20131201T000000_20170129T235959_0301.txt')
difi4 = load_model_swarm_mio_external(path+'SW_OPER_MIO_SHAi2D_20131201T000000_20171231T235959_0401.txt')
difi42 = load_model_swarm_mio_internal(path+'SW_OPER_MIO_SHAi2D_20131201T000000_20171231T235959_0401.txt')

difi3.eval(decimal_year_to_mjd2000(d11), loc, 0, 0, **options)
difi4.eval(decimal_year_to_mjd2000(d11), loc, 0, 0, **options)
loc = (45.0, 105.0, 1.0)
difi4.eval(decimal_year_to_mjd2000(d11), loc, 0, 0, **options)
difi4.eval(decimal_year_to_mjd2000(datetime_to_decimal_year(dt.datetime(2019,1,1,12,0,30))), loc, 0, 0, **options)
difi4.eval(decimal_year_to_mjd2000(datetime_to_decimal_year(dt.datetime(2019,1,1,0,0,30))), loc, 0, 0, **options)
difi4.eval(decimal_year_to_mjd2000(datetime_to_decimal_year(dt.datetime(2019,1,1,13,0,30))), loc, 0, 0, **options)

#get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
import numpy as np
magdifi = np.zeros(3*24).reshape(24,3)
magdifi2 = np.zeros(3*24).reshape(24,3)

for i in range(24):
   t1 = decimal_year_to_mjd2000(datetime_to_decimal_year(dt.datetime(2019,1,1,i,0,30)))
   magdifi[i] = difi4.eval(t1, loc, 0, 0, **options)
   magdifi2[i] = difi42.eval(t1, loc, 0, 0, **options)

plt.title("DIFI-4 IONOSPHERE MAGNETIC FIELD MODEL")
plt.xlabel("UTC time/h")
plt.ylabel("Sq intensity/nT")
plt.plot(magdifi[:,0],'bo', label = "North-X")
plt.plot(magdifi[:,1],'ro', label = "East-Y")
plt.plot(magdifi[:,2],'go', label = "Vertical-Z")
plt.plot(magdifi[:,0] + magdifi2[:,0],'b', label = "North-X2")
plt.plot(magdifi[:,1] + magdifi2[:,1],'r', label = "East-Y2")
plt.plot(magdifi[:,2] + magdifi2[:,2],'g', label = "Vertical-Z2")
plt.legend()
plt.show()
vnorm(magdifi)



