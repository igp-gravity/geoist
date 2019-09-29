#-------------------------------------------------------------------------------
#
#  utilities utilities
#
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
# The above copyright notice and this permission notice shall be included in all
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
# pylint: disable=no-name-in-module

from datetime import datetime
from numpy import sqrt, asarray
from ._pymm import (
    GEODETIC_ABOVE_WGS84, GEOCENTRIC_SPHERICAL, GEOCENTRIC_CARTESIAN,
    convert, vrot_sph2geod, vrot_sph2cart, vrot_cart2sph,
)

SPHERICAL_COORD_TYPES = (GEODETIC_ABOVE_WGS84, GEOCENTRIC_SPHERICAL)

def vrotate(arr, coord_in, coord_out, coord_type_in, coord_type_out):
    """ Rotate vectors from one coordinate system to another.
        Input:
            arr - array of the source vectors
            coord_in - source coordinates
            coord_out - destination coordinates
            coord_type_in - source coordinate system type
            coord_type_out - destination coordinate system type
        Output:
            arr_out - array of the rotated vectors
    """
    # pylint: disable=too-many-return-statements
    if coord_type_in == coord_type_out:
        return arr

    coord_out = None if coord_out is None else asarray(coord_out)
    coord_in = None if coord_in is None else asarray(coord_in)

    if coord_type_in == GEODETIC_ABOVE_WGS84:
        if coord_type_out == GEODETIC_ABOVE_WGS84:
            return arr
        elif coord_type_out == GEOCENTRIC_SPHERICAL:
            return vrot_sph2geod(arr, coord_out[..., 0] - coord_in[..., 0])
        elif coord_type_out == GEOCENTRIC_CARTESIAN:
            return vrot_sph2cart(arr, coord_in[..., 0], coord_in[..., 1])

    elif coord_type_in == GEOCENTRIC_SPHERICAL:
        if coord_type_out == GEODETIC_ABOVE_WGS84:
            return vrot_sph2geod(arr, coord_out[..., 0] - coord_in[..., 0])
        elif coord_type_out == GEOCENTRIC_CARTESIAN:
            return vrot_sph2cart(arr, coord_in[..., 0], coord_in[..., 1])

    elif coord_type_in == GEOCENTRIC_CARTESIAN:
        if coord_type_out in SPHERICAL_COORD_TYPES:
            return vrot_cart2sph(arr, coord_out[..., 0], coord_out[..., 1])

    raise ValueError("Unsupported coordinate system type!")


def vnorm(arr):
    """ Calculate norm for each vector form an input array of vectors:
        |v| = sqrt(dot(v, v))
    """
    arr = asarray(arr)
    return sqrt((arr*arr).sum(axis=arr.ndim-1))


def vincdecnorm(arr):
    """ Convenience function converting magnetic field vector in the local
    northing, easting, elevation (up-pointing) geodetic frame to inclination
    in degrees (from -90 to 90, positive values point down), declination
    in degrees (from -180 to 180, 0 points north, clockwise), and vector
    intensity (vector norm).
    This conversion is equivalent to convert of a Cartesian vector to
    spherical coordinates.
    The function accepts either one vector or an array of vectors and returns
    a tuple of the inclination, declination and norm.
    """
    tmp = convert(arr, GEOCENTRIC_CARTESIAN, GEOCENTRIC_SPHERICAL)
    return -tmp[..., 0], tmp[..., 1], tmp[..., 2]


def datetime_to_decimal_year(time):
    """ Convert time given by a `datetime.datetime` object to a decimal year
    value.
    """
    if not isinstance(time, datetime):
        raise TypeError("The input must be a datetime object.")

    year_start = datetime(year=time.year, month=1, day=1)
    next_year_start = datetime(year=time.year+1, month=1, day=1)

    year_elapsed = (time - year_start).total_seconds()
    year_total = (next_year_start - year_start).total_seconds()

    return time.year + year_elapsed / year_total
