#-------------------------------------------------------------------------------
#
#  Swarm MMA_SHA_2* product file format parser
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

from numpy import array
from spacepy import pycdf

CDF_EPOCH_TYPE = pycdf.const.CDF_EPOCH.value
CDF_EPOCH_2000 = 63113904000000.0
CDF_EPOCH_TO_DAYS = 1.0/86400000.0


def read_swarm_mma_2f_geo_internal(cdf):
    """ Read Swarm MMA_SHA_2F product CDF file and returns the internal
    geographic frame model coefficients.
    """
    return read_swarm_mma_2f_coefficients(cdf, "gh", "geo")


def read_swarm_mma_2f_geo_external(cdf):
    """ Read Swarm MMA_SHA_2F product CDF file and returns the external
    geographic frame model coefficients.
    """
    return read_swarm_mma_2f_coefficients(cdf, "qs", "geo")


def read_swarm_mma_2f_sm_internal(cdf):
    """ Read Swarm MMA_SHA_2F product CDF file and returns the internal
    solar-magnetic frame model coefficients.
    """
    return read_swarm_mma_2f_coefficients(cdf, "gh", "sm")


def read_swarm_mma_2f_sm_external(cdf):
    """ Read Swarm MMA_SHA_2F product CDF file and returns the external
    solar-magnetic frame model coefficients.
    """
    return read_swarm_mma_2f_coefficients(cdf, "qs", "sm")


def read_swarm_mma_2c_internal(cdf):
    """ Read Swarm MMA_SHA_2C product CDF file and returns a tuple
    containing the internal model coefficients.
    """
    return (
        read_swarm_mma_2c_coefficients(cdf, "gh", "1"),
        read_swarm_mma_2c_coefficients(cdf, "gh", "2"),
    )


def read_swarm_mma_2c_external(cdf):
    """ Read Swarm MMA_SHA_2C product CDF file and returns a tuple
    containing the internal model coefficients.
    """
    return (
        read_swarm_mma_2c_coefficients(cdf, "qs", "1"),
        read_swarm_mma_2c_coefficients(cdf, "qs", "2"),
    )


def read_swarm_mma_2f_coefficients(cdf, variable, frame):
    """ Read a set of Swarm MMA_SHA_2F coefficients.

    The function expect a spacepy.pycdf.CDF object.
    """
    return read_swarm_mma_coefficients(
        cdf, "t_" + variable, "nm_" + variable, "%s_%s" % (variable, frame),
        variable
    )


def read_swarm_mma_2c_coefficients(cdf, variable, subset):
    """ Read a single set of Swarm MMA_SHA_2C coefficients.

    The function expect a spacepy.pycdf.CDF object.
    """
    source_variable = "%s_%s" % (variable, subset)
    return read_swarm_mma_coefficients(
        cdf, "t_" + source_variable, "nm_" + source_variable, source_variable,
        variable
    )


def read_swarm_mma_coefficients(cdf, time_variable, nm_variable, coeff_variable,
                                target_variable):
    """ Read a single set of Swarm MMA_SHA_2* coefficients.

    The function expect a spacepy.pycdf.CDF object.
    """
    time = cdf.raw_var(time_variable)
    nm_idx = array(cdf[nm_variable][...])
    if nm_idx.ndim == 3  and nm_idx.shape[0] == 1:
        nm_idx = nm_idx[0]
    return {
        "degree_min": nm_idx[:, 0].min(),
        "degree_max": nm_idx[:, 0].max(),
        "t": _cdf_rawtime_to_mjd2000(time[0], time.type()),
        "nm": nm_idx,
        target_variable: array(cdf[coeff_variable][0].transpose()),
    }


def _cdf_rawtime_to_mjd2000(raw_time, cdf_type):
    """ Convert an array of CDF raw time values to array of MJD2000 values.
    """
    if cdf_type == CDF_EPOCH_TYPE:
        return (raw_time - CDF_EPOCH_2000) * CDF_EPOCH_TO_DAYS
    else:
        raise TypeError("Unsupported CDF time type %r !" % cdf_type)
