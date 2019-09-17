#-------------------------------------------------------------------------------
#
#  IGRF file format model loader
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
from .model import SphericalHarmomicGeomagneticModel
from .coefficients import SparseSHCoefficientsTimeDependentDecimalYear
from .parser_igrf import parse_igrf_file


def load_model_igrf(path):
    """ Load model from an IGRF coefficient file.

    Note: IGRF file format was discontinued after IGRF11. Starting by IGRF12
    the IRFF models uses SHC format which is loaded by the SHC file loader.
    """
    return SphericalHarmomicGeomagneticModel(load_coeff_igrf(path))


def load_coeff_igrf(path):
    """ Load coefficients from an IGRF file. """
    with open(path, encoding="ascii") as file_in:
        data = parse_igrf_file(file_in)

    return SparseSHCoefficientsTimeDependentDecimalYear(
        data["nm"], data["gh"], data["t"]
    )
