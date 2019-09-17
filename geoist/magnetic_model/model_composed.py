#-------------------------------------------------------------------------------
#
#  Aggregated Magnetic Model
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

from collections import namedtuple
from numpy import inf, zeros, asarray
from .._pymm import (
    GEOCENTRIC_SPHERICAL, GEODETIC_ABOVE_WGS84, GEOCENTRIC_CARTESIAN,
    convert, vrot_sph2geod, vrot_sph2cart,
)
from .model import GeomagneticModel


Component = namedtuple("_Component", ["model", "scale", "parameters"])


def _validity_overlap(validity1, validity2):
    start1, end1 = validity1
    start2, end2 = validity2
    start = max(start1, start2)
    end = min(end1, end2)
    if end < start:
        return -inf, -inf
    return start, end


class ComposedGeomagneticModel(GeomagneticModel):
    """ Composed Earth magnetic field model aggregating multiple models
    into one.
    """

    def __init__(self, *models):
        self._parameters = set()
        self._components = []
        self._validity = (-inf, inf)
        for model in models:
            self.push(model)

    def push(self, model, scale=1.0, **parameters):
        """ Add model. """
        self._parameters.update(model.parameters)
        self._validity = _validity_overlap(self.validity, model.validity)
        self._components.append(Component(model, scale, parameters))

    @property
    def validity(self):
        return self._validity

    @property
    def parameters(self):
        """ required parameters. """
        return tuple(self._parameters)

    def eval(self, time, location,
             input_coordinate_system=GEOCENTRIC_SPHERICAL,
             output_coordinate_system=GEOCENTRIC_SPHERICAL,
             **options):

        # convert input coordinates to spherical coordinates
        coord_sph = convert(
            location, input_coordinate_system, GEOCENTRIC_SPHERICAL
        )

        # get output dimension
        time = asarray(time)
        location = asarray(location)
        if time.ndim > (location.ndim - 1):
            shape = time.shape
        else:
            shape = location.shape[:-1]

        result = zeros(shape + (3,))
        final_scale = options.pop("scale", None)
        for model, scale, params in self._components:
            args = options.copy()
            args.update(params)
            result += model.eval(time, coord_sph, scale=scale, **args)

        # rotate result to the desired coordinate frame
        if output_coordinate_system == GEODETIC_ABOVE_WGS84:
            if input_coordinate_system == GEODETIC_ABOVE_WGS84:
                coord_out = location
            else:
                coord_out = convert(
                    coord_sph, GEOCENTRIC_SPHERICAL, GEODETIC_ABOVE_WGS84
                )
            result = vrot_sph2geod(result, coord_out[..., 0] - coord_sph[..., 0])
        elif output_coordinate_system == GEOCENTRIC_CARTESIAN:
            result = vrot_sph2cart(result, coord_sph[..., 0], coord_sph[..., 1])

        # apply the final scale
        if final_scale is not None:
            result *= final_scale

        return result
