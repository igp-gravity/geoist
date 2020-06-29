# -*- coding: utf-8 -*-
# Name        : giconstants.py
# Created on  : 2018/11/24 17:00
# Author      : Steve Chen <chenshi@cea-igp.ac.cn>
# Affiliation : Institute of Geophysics, CEA.
# Version     : 0.1.0
# Copyright   : Copyright (C) 2018-2020 GEOIST Development Team. All Rights Reserved.
# License     : Distributed under the MIT License. See LICENSE.txt for more info.
# Github      : https://igp-gravity.github.io/
"""
   Holds all physical constants and unit conversions used in GEOIST.
   All modules should import the constants from here!
   All constants should be in SI, unless otherwise stated!
"""


#: The default thermal diffusivity in :math:`m^2/s`
THERMAL_DIFFUSIVITY = 0.000001

#: The default thermal diffusivity but in :math:`m^2/year`
THERMAL_DIFFUSIVITY_YEAR = 31.5576

#: Conversion factor from SI units to Eotvos: :math:`1/s^2 = 10^9\ Eotvos`
SI2EOTVOS = 1000000000.0

#: Conversion factor from SI units to mGal: :math:`1\ m/s^2 = 10^5\ mGal`
SI2MGAL = 100000.0

#: The gravitational constant in :math:`m^3 kg^{-1} s^{-1}`
G = 0.00000000006673

#: Proportionality constant used in the magnetic method in henry/m (SI)
CM = 10. ** (-7)

#: Conversion factor from tesla to nanotesla
T2NT = 10. ** (9)

#: The mean earth radius in meters
MEAN_EARTH_RADIUS = 6378137.0

#: Permeability of free space in :math:`N A^{-2}`
PERM_FREE_SPACE = 4 * \
    3.141592653589793115997963468544185161590576171875 * (10 ** -7)
