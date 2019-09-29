#-------------------------------------------------------------------------------
#
#  EMM model loader
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
from .coefficients import (
    SparseSHCoefficientsTimeDependentDecimalYear,
    SparseSHCoefficientsConstant,
    CombinedSHCoefficients,
)
from .parser_emm import combine_emm_coefficients, parse_emm_file


def load_model_emm(path_static, path_secvar):
    """ Load model from a EMM coefficient files. """
    return SphericalHarmomicGeomagneticModel(
        load_coeff_emm(path_static, path_secvar)
    )


def load_coeff_emm(path_static, path_secvar):
    """ Load coefficients from a EMM coefficient files. """

    with open(path_static, encoding="ascii") as file_static:
        with open(path_secvar, encoding="ascii") as file_secvar:
            data_variable, data_constant = combine_emm_coefficients(
                parse_emm_file(file_static),
                parse_emm_file(file_secvar),
            )

    return CombinedSHCoefficients(
        SparseSHCoefficientsTimeDependentDecimalYear(
            data_variable["nm"], data_variable["gh"], data_variable["t"],
        ),
        SparseSHCoefficientsConstant(
            data_constant["nm"], data_constant["gh"][:, 0]
        )
    )
