#-------------------------------------------------------------------------------
#
#  EMM model loader
#
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
