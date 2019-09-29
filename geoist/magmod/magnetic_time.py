#-------------------------------------------------------------------------------
#
#  Magnetic time calculations
#
# Author: Martin Paces <martin.paces@eox.at>
#
#-------------------------------------------------------------------------------
# Copyright (C) 2018 EOX IT Services GmbH
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies of this Software or works derived from this Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#-------------------------------------------------------------------------------

from math import pi
from numpy import sin, cos, arctan2
from .solar_position import sunpos

DEG2RAD = pi/180.0
RAD2HOUR = 12.0/pi


def mjd2000_to_magnetic_universal_time(mjd2000, lat_ngp, lon_ngp,
                                       lat_sol=None, lon_sol=None):
    """ Evaluate magnetic universal time for the given MJD2000 and
    coordinates of the North Geomagnetic Pole.

    The magnetic universal time is geometrically equivalent
    to the magnetic dipole longitude of the sub-solar point.

    Function allows specification of user defined sub-solar
    latitude and longitude lat_sol and lon_sol in degrees
    overriding the default Sun model.
    """
    if lat_sol is None or lon_sol is None:
        lat_sol, lon_sol = get_subsol(mjd2000)

    return _mjd2000_to_magnetic_universal_time(
        mjd2000, lat_ngp, lon_ngp, lat_sol, lon_sol
    )


def get_subsol(mjd2000):
    """ Calculate sub-solar point coordinates for the given MJD2000 time. """
    declination, _, hour_angle, _, _ = sunpos(mjd2000, 0, 0, rad=0)
    return declination, -hour_angle


def _mjd2000_to_magnetic_universal_time(mjd2000, lat_ngp, lon_ngp,
                                        lat_sol, lon_sol):
    latr_sol = DEG2RAD * lat_sol
    latr_ngp = DEG2RAD * lat_ngp
    dif_lonr = DEG2RAD * (lon_sol - lon_ngp)

    tmp_y = cos(latr_sol)
    tmp_x = sin(latr_ngp)*tmp_y*cos(dif_lonr) - cos(latr_ngp)*sin(latr_sol)
    tmp_y *= sin(dif_lonr)

    return RAD2HOUR * (pi - arctan2(tmp_y, tmp_x))
