#-------------------------------------------------------------------------------
#
#  MIO Magnetic Model
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
# pylint: disable=too-many-arguments,too-many-locals

from numpy import asarray, empty, nditer, ndim, isnan, full, newaxis, ones
from .._pymm import GEOCENTRIC_SPHERICAL, convert
from ..magnetic_time import mjd2000_to_magnetic_universal_time
from .model import GeomagneticModel, DipoleSphericalHarmomicGeomagneticModel

MIO_HEIGHT = 110.0 # km
MIO_EARTH_RADIUS = 6371.2 # km
MIO_WOLF_RATIO = 0.014850


class DipoleMIOPrimaryGeomagneticModel(GeomagneticModel):
    """ Composed model switching between MIO primary fields evaluation
    above and below the Ionosphere.
    """
    parameters = ("time", "location", "f107", "subsolar_point")

    @property
    def degree(self):
        """ Get maximum degree of the model. """
        return max(
            self.model_below_ionosphere.degree,
            self.model_above_ionosphere.degree
        )

    @property
    def min_degree(self):
        """ Get minimum degree of the model. """
        return max(
            self.model_below_ionosphere.min_degree,
            self.model_above_ionosphere.min_degree
        )

    def __init__(self, model_below_ionosphere, model_above_ionosphere,
                 height=MIO_HEIGHT, earth_radius=MIO_EARTH_RADIUS):
        self.model_below_ionosphere = model_below_ionosphere
        self.model_above_ionosphere = model_above_ionosphere
        self.height = height
        self.earth_radius = earth_radius

    def eval(self, time, location,
             input_coordinate_system=GEOCENTRIC_SPHERICAL,
             output_coordinate_system=GEOCENTRIC_SPHERICAL,
             **options):
        time = asarray(time)
        location = convert(
            location, input_coordinate_system, GEOCENTRIC_SPHERICAL
        )

        radius_ionosphere = self.height + self.earth_radius
        radius = location[..., 2]
        mask_below_ionosphere = radius <= radius_ionosphere

        result = empty(location.shape)

        def _eval_model(model, mask, **options):

            def _mask_array(data):
                return data[mask] if data.ndim else data

            def _mask_option(key):
                data = options.get(key)
                if data is not None:
                    options[key] = _mask_array(asarray(data))

            _mask_option('f107')
            _mask_option('lat_sol')
            _mask_option('lon_sol')

            result[mask] = model.eval(
                _mask_array(time), location[mask], GEOCENTRIC_SPHERICAL,
                output_coordinate_system, **options
            )

        _eval_model(
            self.model_below_ionosphere, mask_below_ionosphere, **options
        )
        _eval_model(
            self.model_above_ionosphere, ~mask_below_ionosphere, **options
        )

        return result

    @property
    def validity(self):
        start_above, stop_above = self.model_above_ionosphere.validity
        start_below, stop_below = self.model_below_ionosphere.validity
        return (max(start_above, start_below), min(stop_above, stop_below))


class DipoleMIOGeomagneticModel(DipoleSphericalHarmomicGeomagneticModel):
    """ Swarm MIO model calculated by the Spherical Harmonic
    Expansion in the dipole coordinates.
    The dipole coordinates are defined by the north pole latitude and longitude.
    The model requires time dependent F10.7 index values.

    Options:
        coefficients - spherical harmonic coefficients
        north_pole - a tuple of dipole north pole lat/lon coordinates
                     or callable returning a tuple of dipole north pole
                     lat/lon coordinates for a given MJD2000 time.
        wolf_ratio - Wolf ratio (F10.7 scale)
        height - radius of the ionosphere (a + h)
        earth_radius - mean Earth radius used by the MIO model.
    """
    # list of the required model parameters
    parameters = ("time", "location", "f107", "subsolar_point")

    def __init__(self, coefficients, north_pole, wolf_ratio=MIO_WOLF_RATIO,
                 height=MIO_HEIGHT, earth_radius=MIO_EARTH_RADIUS):
        DipoleSphericalHarmomicGeomagneticModel.__init__(
            self, coefficients, north_pole
        )
        self.wolf_ratio = wolf_ratio

    def eval(self, time, location,
             input_coordinate_system=GEOCENTRIC_SPHERICAL,
             output_coordinate_system=GEOCENTRIC_SPHERICAL,
             **options):
        time = asarray(time)
        location = asarray(location)

        lat_ngp, lon_ngp = self.north_pole
        mut = mjd2000_to_magnetic_universal_time(
            time, lat_ngp, lon_ngp,
            lat_sol=options.pop('lat_sol', None),
            lon_sol=options.pop('lon_sol', None),
        )

        mio_scale = 1.0 + self.wolf_ratio * asarray(options.pop('f107'))

        mask = self.coefficients.is_valid(time) & ~isnan(mio_scale)

        if time.ndim:
            if not mio_scale.ndim:
                mio_scale = full(time.shape, mio_scale)
            mio_scale = mio_scale[mask, newaxis]
            mut = mut[mask]

        scale = (options.pop('scale', 1.0) * ones(3)) * mio_scale

        return self._eval_masked(
            mask, time, location,
            input_coordinate_system, output_coordinate_system,
            scale=scale, mut=mut, **options
        )

    def _eval_multi_time(self, time, coords, input_coordinate_system,
                         output_coordinate_system, mut, **options):
        result = empty(coords.shape)
        if result.size > 0:
            iterator = nditer(
                [
                    result[..., 0], result[..., 1], result[..., 2],
                    time, coords[..., 0], coords[..., 1], coords[..., 2],
                    mut,
                ],
                op_flags=[
                    ['writeonly'], ['writeonly'], ['writeonly'],
                    ['readonly'], ['readonly'], ['readonly'], ['readonly'],
                    ['readonly'],
                ],
            )
            for item in iterator:
                (
                    vect0, vect1, vect2,
                    time_, coord0, coord1, coord2, mut_,
                ) = item
                vect0[...], vect1[...], vect2[...] = self._eval_single_time(
                    time_, [coord0, coord1, coord2],
                    input_coordinate_system, output_coordinate_system,
                    mut=mut_, **options
                )
        return result
