#-------------------------------------------------------------------------------
#
#  Magnetic Model
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

from numpy import asarray, empty, full, nan, nditer, ndim
from .._pymm import GRADIENT, GEOCENTRIC_SPHERICAL, sheval
from ..sheval_dipole import rotate_vectors_from_dipole
from ..dipole_coords import convert_to_dipole


class GeomagneticModel(object):
    """ Abstract base class of the Earth magnetic field model. """
    # list of the required model parameters
    parameters = ("time", "location")

    def eval(self, time, location,
             input_coordinate_system=GEOCENTRIC_SPHERICAL,
             output_coordinate_system=GEOCENTRIC_SPHERICAL,
             **options):
        """ Evaluate magnetic field for the given MJD2000 times and coordinates.
        """
        raise NotImplementedError

    @property
    def validity(self):
        """ Get model's validity range as a tuple of two MJD2000 times.
        In case of an unconstrained validity rage (-inf, +inf) tuple is
        returned.
        """
        raise NotImplementedError


class SphericalHarmomicGeomagneticModel(GeomagneticModel):
    """ Earth magnetic field model calculated by the Spherical Harmonic
    Expansion.
    """

    def __init__(self, coefficients):
        self.coefficients = coefficients

    @property
    def validity(self):
        return self.coefficients.validity

    @property
    def degree(self):
        """ Get maximum degree of the model. """
        return self.coefficients.degree

    @property
    def min_degree(self):
        """ Get minimum degree of the model. """
        return self.coefficients.min_degree

    def eval(self, time, location,
             input_coordinate_system=GEOCENTRIC_SPHERICAL,
             output_coordinate_system=GEOCENTRIC_SPHERICAL,
             **options):
        time = asarray(time)
        location = asarray(location)

        mask = self.coefficients.is_valid(time)

        return self._eval_masked(
            mask, time, location, input_coordinate_system,
            output_coordinate_system, **options
        )

    def _eval_masked(self, mask, time, location, input_coordinate_system,
                     output_coordinate_system, **options):
        if ndim(time):
            eval_model = self._eval_multi_time
            time = time[mask]
        else:
            eval_model = self._eval_single_time

        result = full(location.shape, nan)
        result[mask] = eval_model(
            time, location[mask], input_coordinate_system,
            output_coordinate_system, **options
        )
        return result

    def _eval_multi_time(self, time, coords, input_coordinate_system,
                         output_coordinate_system, **options):
        """ Evaluate spherical harmonic for multiple times. """
        result = empty(coords.shape)
        if result.size > 0:
            iterator = nditer(
                [
                    time, coords[..., 0], coords[..., 1], coords[..., 2],
                    result[..., 0], result[..., 1], result[..., 2],
                ],
                op_flags=[
                    ['readonly'], ['readonly'], ['readonly'], ['readonly'],
                    ['writeonly'], ['writeonly'], ['writeonly'],
                ],
            )
            for time_, coord0, coord1, coord2, vect0, vect1, vect2 in iterator:
                vect0[...], vect1[...], vect2[...] = self._eval_single_time(
                    time_, [coord0, coord1, coord2], input_coordinate_system,
                    output_coordinate_system, **options
                )
        return result

    def _eval_single_time(self, time, coords, input_coordinate_system,
                          output_coordinate_system, **options):
        """ Evaluate spherical harmonic for a single time."""
        coefficients = self.coefficients
        coeff, degree = coefficients(time, **options)
        return sheval(
            coords, degree, coeff[..., 0], coeff[..., 1],
            is_internal=coefficients.is_internal, mode=GRADIENT,
            coord_type_in=input_coordinate_system,
            coord_type_out=output_coordinate_system,
            scale_gradient=-asarray(options.get("scale", 1.0))
        )


class DipoleSphericalHarmomicGeomagneticModel(SphericalHarmomicGeomagneticModel):
    """ Earth magnetic field model calculated by the Spherical Harmonic
    Expansion in the dipole coordinates.
    The dipole coordinates are defined north latitude and longitude of
    the dipole axis.
    """

    def __init__(self, coefficients, north_pole):
        SphericalHarmomicGeomagneticModel.__init__(self, coefficients)
        self.north_pole = north_pole

    def _eval_masked(self, mask, time, location, input_coordinate_system,
                     output_coordinate_system, **options):
        lat_ngp, lon_ngp = self.north_pole
        scale = options.pop('scale', None)

        result = full(location.shape, nan)

        if ndim(time):
            eval_model = self._eval_multi_time
            time = time[mask]
        else:
            eval_model = self._eval_single_time

        location = location[mask]
        location_dipole = convert_to_dipole(
            location, lat_ngp, lon_ngp, input_coordinate_system
        )

        result[mask] = rotate_vectors_from_dipole(
            eval_model(
                time, location_dipole,
                GEOCENTRIC_SPHERICAL, GEOCENTRIC_SPHERICAL,
                **options
            ),
            lat_ngp, lon_ngp, location_dipole, location,
            input_coordinate_system, output_coordinate_system
        )

        if scale is not None:
            result[mask] *= scale

        return result
