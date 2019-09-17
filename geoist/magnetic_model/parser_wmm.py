#-------------------------------------------------------------------------------
#
#  World Magnetic Model file format parser
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

import re
from numpy import array

RE_TERMINATING_LINE = re.compile(r'^9+$')
WMM_VALIDITY_PERIOD = 5.0 # years


def parse_wmm_file(file_in):
    """ Parse World Magnetic Model COF file format and return a dictionary
    containing the parsed model data.
    """
    lines = file_in
    data = parse_wmm_header(next(lines))
    data["nm"], data["gh"], data["t"] = parse_wmm_coefficients(lines, data)
    data["degree_min"] = data["nm"][:, 0].min()
    data["degree_max"] = data["nm"][:, 0].max()
    return data


def parse_wmm_coefficients(lines, data):
    """ Parse the WMM COF coefficients. """
    nm_index = []
    coeff = []

    # parse coefficients and convert them to a sparse structure
    for line in lines:
        fields = line.split()
        if len(fields) == 1 and RE_TERMINATING_LINE.match(fields[0]):
            break
        n_idx, m_idx = int(fields[0]), int(fields[1])
        coef_g, coef_h, coef_dg, coef_dh = [float(v) for v in fields[2:6]]
        if coef_g != 0 or coef_dg != 0:
            nm_index.append((n_idx, m_idx))
            coeff.append((coef_g, coef_dg))
        if m_idx > 0 and (coef_h != 0 or coef_dh != 0):
            nm_index.append((n_idx, -m_idx))
            coeff.append((coef_h, coef_dh))

    # convert (t0, dt) to (t0, t1)
    epoch = data["epoch"]
    times = array([epoch, epoch + WMM_VALIDITY_PERIOD])
    nm_index = array(nm_index)
    coeff = array(coeff)
    coeff[:, 1] = coeff[:, 0] + coeff[:, 1]*WMM_VALIDITY_PERIOD

    return nm_index, coeff, times


def parse_wmm_header(line):
    """ Parse the WMM COF file header. """
    fields = line.split()
    return {
        "epoch": float(fields[0]),
        "name": fields[1],
        "version": fields[2],
    }
