#-------------------------------------------------------------------------------
#
#  Geomagnetic model - spherical harmonic evaluation in the dipole coordinates.
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
# pylint: disable=no-name-in-module

from ._pymm import (
    GEODETIC_ABOVE_WGS84, GEOCENTRIC_CARTESIAN, GEOCENTRIC_SPHERICAL,
    POTENTIAL, GRADIENT, convert, sheval,
)
from .dipole_coords import convert_to_dipole, vrot_from_dipole


def sheval_dipole(arr_in, degree, coef_g, coef_h, lat_ngp, lon_ngp,
                  coord_type_in=GEODETIC_ABOVE_WGS84,
                  coord_type_out=GEODETIC_ABOVE_WGS84,
                  mode=GRADIENT, is_internal=True,
                  scale_potential=1.0, scale_gradient=1.0):
    """
    Evaluate spherical harmonic model in the dipole coordinate frame given
    by the provided North Geomagnetic Pole coordinates.

    Parameters:
       arr_in - array of 3D coordinates (up to 16 dimensions).
       degree - degree of the spherical harmonic model.
       coef_g - vector of spherical harmonic model coefficients.
       coef_h - vector of spherical harmonic model coefficients.
       lat_ngp - North Geomagnetic Pole latitude.
       lon_ngp - North Geomagnetic Pole longitude.
       coord_type_in - type of the input coordinates.
       mode - quantity to be evaluated:
                  POTENTIAL
                  GRADIENT (default)
                  POTENTIAL_AND_GRADIENT
       rad_ref - reference (Earth) radius
       is_internal - boolean flag set to True by default. When set to False
                     external field evaluation is used.
       scale_potential - scalar value multiplied with the result potentials.
       scale_gradient - scalar or 3 element array multiplied with the result
                        gradient components.
    """
    arr_in_dipole = convert_to_dipole(arr_in, lat_ngp, lon_ngp, coord_type_in)
    result = sheval(
        arr_in_dipole, mode=mode, is_internal=is_internal,
        degree=degree, coef_g=coef_g, coef_h=coef_h,
        coord_type_in=GEOCENTRIC_SPHERICAL,
        coord_type_out=GEOCENTRIC_SPHERICAL,
        scale_potential=scale_potential,
    )

    if mode == GRADIENT:
        return rotate_vectors_from_dipole(
            result, lat_ngp, lon_ngp, arr_in_dipole, arr_in,
            coord_type_in, coord_type_out,
        ) * scale_gradient
    elif mode == POTENTIAL:
        return result
    else: #mode == POTENTIAL_AND_GRADIENT
        potential, gradient = result
        return potential, rotate_vectors_from_dipole(
            gradient, lat_ngp, lon_ngp, arr_in_dipole, arr_in,
            coord_type_in, coord_type_out,
        ) * scale_gradient


def rotate_vectors_from_dipole(vectors, lat_ngp, lon_ngp,
                               coords_dipole, coord_in,
                               coord_type_in=GEODETIC_ABOVE_WGS84,
                               coord_type_out=GEODETIC_ABOVE_WGS84):
    """ Rotate vectors from the Dipole spherical coordinate frame
    to the requested output coordinate frame.
    """
    lat_dipole, lon_dipole = coords_dipole[..., 0], coords_dipole[..., 1]
    if coord_type_out == GEOCENTRIC_CARTESIAN:
        lat_out, lon_out = None, None
    else:
        coords_out = convert(coord_in, coord_type_in, coord_type_out)
        lat_out, lon_out = coords_out[..., 0], coords_out[..., 1]
        del coords_out

    return vrot_from_dipole(
        vectors, lat_ngp, lon_ngp, lat_dipole, lon_dipole, lat_out, lon_out,
        coord_type_out
    )
