#-------------------------------------------------------------------------------
#
#  IGRF format parser
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
# The above copyright notice and this permission notice shall be included in
# all
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

IGRF_EXTRAPOLATION_PERIOD = 5.0 # years


def parse_igrf_file(file_in):
    """ Parse IGRF file format and return a dictionary containing the parsed
    model data.
    """
    lines = strip_igrf_comments(file_in)
    data = {}
    data["labels"] = parse_igrf_header(next(lines))
    data["t"] = parse_igrf_times(next(lines))
    data["nm"], data["gh"] = parse_igrf_coefficients(lines)
    data["degree_min"] = data["nm"][:, 0].min()
    data["degree_max"] = data["nm"][:, 0].max()
    return data


def parse_igrf_coefficients(lines):
    """ Parse IGRF coefficients. """
    nm_idx = []
    coeff = []
    for line in lines:
        fields = line.split()
        label, n_idx, m_idx = fields[0], int(fields[1]), int(fields[2])
        if label == "h":
            m_idx = -m_idx
        nm_idx.append((n_idx, m_idx))
        coeff.append([float(v) for v in fields[3:]])
    coeff = array(coeff)
    coeff[:, -1] = coeff[:, -2] + coeff[:, -1] * IGRF_EXTRAPOLATION_PERIOD
    return array(nm_idx), coeff


def parse_igrf_times(line):
    """ Parse SHC times. """
    times = [float(v) for v in line.split()[3:-1]]
    times.append(times[-1] + IGRF_EXTRAPOLATION_PERIOD)
    return array(times)


def parse_igrf_header(line):
    """ Parse IGRF header with the column labels. """
    return line.split()


def strip_igrf_comments(file_in):
    """ Strip initial comments and empty lines from a text file stream. """
    for line in file_in:
        line = line.partition("#")[0].strip()
        if line:
            yield line
            break
    for line in file_in:
        line = line.strip()
        yield line
