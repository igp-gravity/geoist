#-------------------------------------------------------------------------------
#
#  WMM model loader
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
from .parser_wmm import parse_wmm_file


def load_model_wmm(path):
    """ Load model from a WMM COF file. """
    return SphericalHarmomicGeomagneticModel(load_coeff_wmm(path))


def load_coeff_wmm(path):
    """ Load coefficients from a WMM COF file. """
    with open(path, encoding="ascii") as file_in:
        data = parse_wmm_file(file_in)

    return SparseSHCoefficientsTimeDependentDecimalYear(
        data["nm"], data["gh"], data["t"]
    )
