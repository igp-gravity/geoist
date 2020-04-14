#!/usr/bin/env python
# coding: utf-8

# #### MAGNETIC FIELD MODEL COLLETCTION
# * Auther: Shi Chen at CEA-IGP
# * Date: 2019-09-29

import datetime as dt

from geoist.magmod.data import DIFI4
from geoist.magmod.magnetic_model.loader_mio import (
    load_model_swarm_mio_internal,
    load_model_swarm_mio_external)

from geoist.magmod.time_util import decimal_year_to_mjd2000
from geoist.magmod.util import datetime_to_decimal_year

loc = (45.0, 105.0, 1.0)
options = {"f107": 70, "scale": [1, 1, -1]}  #-1 is Z direction

# load DIFI4 model
difi4 = load_model_swarm_mio_external(DIFI4)
difi42 = load_model_swarm_mio_internal(DIFI4)

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


