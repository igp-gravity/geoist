#-------------------------------------------------------------------------------
#
#  Swarm MIO_SHA_2* coefficients loaders
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

from io import open
from numpy import arange
from .model_mio import (
    DipoleMIOPrimaryGeomagneticModel, DipoleMIOGeomagneticModel,
    MIO_EARTH_RADIUS,
)
from .coefficients_mio import SparseSHCoefficientsMIO
from .parser_mio import parse_swarm_mio_file


def load_model_swarm_mio_internal(path):
    """ Load internal (secondary field) model from a Swarm MIO_SHA_2* product.
    """
    coefficients, params = load_coeff_swarm_mio_internal(path)
    return _create_mio_model(coefficients, params)


def load_model_swarm_mio_external(path, above_ionosphere=None):
    """ Load external (primary field) model from a Swarm MIO_SHA_2* product.
    """
    with open(path, encoding="ascii") as file_in:
        params = parse_swarm_mio_file(file_in)

    if above_ionosphere is None:
        return _create_composed_mio_model(
            _get_coeff_swarm_mio_external(params, False),
            _get_coeff_swarm_mio_external(params, True),
            params
        )
    else:
        return _create_mio_model(
            _get_coeff_swarm_mio_external(params, above_ionosphere), params
        )


def _create_composed_mio_model(coefficients_below_ionosphere,
                               coefficients_above_ionosphere, params):
    return DipoleMIOPrimaryGeomagneticModel(
        _create_mio_model(coefficients_below_ionosphere, params),
        _create_mio_model(coefficients_above_ionosphere, params),
        height=params["height"],
    )


def _create_mio_model(coefficients, params):
    return DipoleMIOGeomagneticModel(
        coefficients, north_pole=(params["lat_NGP"], params["lon_NGP"]),
        wolf_ratio=params["wolf_ratio"], height=params["height"],
    )


def load_coeff_swarm_mio_internal(path):
    """ Load internal model coefficients and other parameters
    from a Swarm MIO_SHA_2* product file.
    """
    with open(path, encoding="ascii") as file_in:
        data = parse_swarm_mio_file(file_in)

    return SparseSHCoefficientsMIO(
        data["nm"], data["gh"],
        ps_extent=(data["pmin"], data["pmax"], data["smin"], data["smax"]),
        is_internal=True,
    ), data


def load_coeff_swarm_mio_external(path, above_ionosphere=True):
    """ Load external model coefficients from a Swarm MIO_SHA_2* product file.
    Use the `above_ionosphere` to pick the right model variant.
    """
    with open(path, encoding="ascii") as file_in:
        data = parse_swarm_mio_file(file_in)

    return _get_coeff_swarm_mio_external(data, above_ionosphere), data


def _get_coeff_swarm_mio_external(data, above_ionosphere):
    """ Create coefficient object for the given source data. """
    indices = data["nm"]
    coefficients = data["qs"]
    if above_ionosphere:
        is_internal = True
        coefficients = convert_external_mio_coeff(
            data["degree_max"], indices, coefficients, data["height"]
        )
    else:
        is_internal = False

    return SparseSHCoefficientsMIO(
        indices, coefficients,
        ps_extent=(data["pmin"], data["pmax"], data["smin"], data["smax"]),
        is_internal=is_internal,
    )


def convert_external_mio_coeff(degree, indices, coefficients, height):
    """ Convert external coefficients to internal ones. """
    nrrad = -(1.0 + height/MIO_EARTH_RADIUS)
    order = arange(degree + 1, dtype='float')
    scale = (order/(order + 1)) * nrrad**(2*order + 1)
    return (scale[indices[:, 0]] * coefficients.transpose()).transpose()
