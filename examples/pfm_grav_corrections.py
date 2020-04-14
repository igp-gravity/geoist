# -*- coding: utf-8 -*-
"""
Created on Sun Apr 12 22:46:16 2020

@author: chens
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
from geoist.pfm.hawaii_gravity import fetch_hawaii_gravity
from geoist.pfm import normgra, GRS80

lat = 45.0
height = 1000.0
grav = 976783.05562266 #mGal

datapath = Path(normgra.__file__).parent
filename = Path(datapath, 'data', 'yngravity.csv')
gradata = pd.read_csv(filename)
print(gradata.keys())

gamma0 = normgra.gamma_closed_form(lat, 0)
gamma = normgra.gamma_closed_form(lat, height)
gamma1 = normgra.gamma_closed_form(lat, height, ellipsoid = GRS80)

somi0 = normgra.gamma_somigliana_free_air(lat, 0)
somi = normgra.gamma_somigliana_free_air(lat, height)
somi1 = normgra.gamma_somigliana_free_air(lat, height, ellipsoid = GRS80)
bug = normgra.bouguer_plate(height)
print(gamma, gamma1)
print(somi, somi1)
print(bug)
print(gamma, somi, gamma- gamma0)
print(gamma1, somi1, somi - somi0)
# Use a Bouguer plate to remove the effect of topography
bouguer = grav -gamma - bug
print(bouguer, bug)


data = fetch_hawaii_gravity()

gamma = normgra.gamma_closed_form(data['lat'], data['height'])
disturbance = data['gravity'] - gamma
bouguer = disturbance - normgra.bouguer_plate(data['topography'])
shape = data['shape']
x, y = data['x'].reshape(shape), data['y'].reshape(shape)