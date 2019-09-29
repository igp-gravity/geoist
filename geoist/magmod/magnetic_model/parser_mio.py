#-------------------------------------------------------------------------------
#
#  Swarm MIO_SHA_2* product file format parser
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
from .parser_shc import _strip_shc_comments


def parse_swarm_mio_file(file_in):
    """ Parse Swarm MIO_SHA_2* product file format and return a dictionary
    containing the parsed model data.
    """
    lines = _strip_shc_comments(file_in)
    data = parse_swarm_mio_header(next(lines))
    data["nm"], data["qs"], data["gh"] = parse_swarm_mio_coefficients(lines, data)
    data["degree_min"] = data["nm"][:, 0].min()
    data["degree_max"] = data["nm"][:, 0].max()
    return data

def parse_swarm_mio_coefficients(lines, data):
    """ Parse the Swarm MIO_SHA_2* coefficients. """
    nm_idx = []
    coeff = []
    for line in lines:
        fields = line.split()
        nm_idx.append([int(v) for v in fields[:2]])
        coeff.append([float(v) for v in fields[2:]])
    nm_idx = array(nm_idx)
    coeff = array(coeff)

    # shape of the 2D Fourier series coefficients arrays
    coeff_shape = (
        len(coeff), data["smax"] - data["smin"] + 1,
        data["pmax"] - data["pmin"] + 1, 2
    )
    coeff = coeff.reshape(coeff_shape)

    # split internal and external coefficients
    size = nm_idx.shape[0] >> 1
    if (nm_idx[:size] != nm_idx[size:]).any():
        raise ValueError("Mismatch between the QS and HG coefficients!")

    return nm_idx[:size], coeff[:size], coeff[size:]


def parse_swarm_mio_header(line):
    """ Parse the Swarm MIO_SHA_2* file header. """
    fields = line.split()
    colat_ngp = float(fields[6])
    lat_ngp = float(fields[7])
    return {
        "nmax": int(fields[0]),
        "mmax": int(fields[1]),
        "pmin": int(fields[2]),
        "pmax": int(fields[3]),
        "smin": int(fields[4]),
        "smax": int(fields[5]),
        "theta_NGP": colat_ngp,
        "phi_NGP": lat_ngp,
        "lat_NGP": 90.0 - colat_ngp,
        "lon_NGP": 180 - ((180 - lat_ngp) % 360),
        "height": float(fields[8]),
        "wolf_ratio": float(fields[9]),
    }
