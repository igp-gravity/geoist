#-------------------------------------------------------------------------------
#
#  Solar position calculation
#
#  Deprecated old interface.
#
# Project: VirES
# Author: Martin Paces <martin.paces@eox.at>
#
#-------------------------------------------------------------------------------
# Copyright (C) 2017 EOX IT Services GmbH
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all
# copies of this Software or works derived from this Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#-------------------------------------------------------------------------------

from . import _pysunpos


def sunpos(time_mjd2k, lat, lon, rad=6371.2, dtt=0):
    """ Calculate solar equatorial and horizontal coordinates
    for given MJD2000 times and geocentric coordinates (lat, lon, rad).

    Outputs: arrays of the Sun equatorial and horizontal coordinates:
      declination
      right ascension
      local hour angle (global hour angle for 0, 0 latitude and longitude)
      local azimuth
      local zenith

      Parameters:
        time_mjd2k - array of MJD2000 times (up to 15 dimensions).
        lat - array of latitudes [deg]
        lon - array of longitudes [deg]
        rad - array of radii [km] (set to 0 to disable the parallax correction)
        dtt - array of offsets to TT [sec]
    """
    result = _pysunpos.sunpos(time_mjd2k, lat, lon, rad, dtt)
    return tuple(result[..., idx] for idx in range(result.shape[-1]))


def sunpos_original(time_mjd2k, lat, lon, rad=6371.2, dtt=0, pres=1.0, temp=20.0):
    """ Calculate solar equatorial and horizontal coordinates
    for given MJD2000 times and geocentric coordinates (lat, lon, rad).
    This is the original implementation.

    Outputs: arrays of the Sun equatorial and horizontal coordinates:
      declination
      right ascension
      local hour angle (global hour angle for 0, 0 latitude and longitude)
      local azimuth
      local zenith

      Parameters:
        time_mjd2k - array of MJD2000 times (up to 15 dimensions).
        lat - array of latitudes [deg]
        lon - array of longitudes [deg]
        rad - array of radii [km] (set to 0 to disable the parallax correction)
        dtt - array of offsets to TT [sec]
        pres - array of offsets of pressures [atm] (refraction correction)
        temp - array of offsets to temperatures [dgC] (refraction coorection)
    """
    result = _pysunpos.sunpos_original(
        time_mjd2k, lat, lon, rad, dtt, pres, temp
    )
    return tuple(result[..., idx] for idx in range(result.shape[-1]))
