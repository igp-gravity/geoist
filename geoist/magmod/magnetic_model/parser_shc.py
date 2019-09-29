#-------------------------------------------------------------------------------
#
#  SHC format parser
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

from numpy import inf, array


def parse_shc_file(file_in):
    """ Parse SHC file and return a dictionary containing the parsed model data.
    """
    lines = _strip_shc_comments(file_in)
    data = _parse_shc_header(lines)
    data["nm"], data["gh"] = _parse_shc_coefficients(lines)
    return data


def parse_shc_header(file_in):
    """ Parse SHC file header and return a dictionary containing the parsed
    model info.
    """
    lines = _strip_shc_comments(file_in)
    return _parse_shc_header(lines)


def _parse_shc_header(lines):
    data = _parse_shc_header_line(next(lines))
    data["t"] = times = _parse_shc_times(next(lines))
    if "validity_start" not in data:
        data["validity_start"] = times.min() if times.size > 1 else -inf
    if "validity_end" not in data:
        data["validity_end"] = data["t"].max() if times.size > 1 else +inf
    return data


def _parse_shc_coefficients(lines):
    """ Parse SHC coefficients. """
    nm_index = []
    coefficients = []
    for line in lines:
        fields = line.split()
        nm_index.append([int(v) for v in fields[:2]])
        coefficients.append([float(v) for v in fields[2:]])
    return array(nm_index), array(coefficients)


def _parse_shc_times(line):
    """ Parse SHC times. """
    return array([float(v) for v in line.split()])


def _parse_shc_header_line(line):
    """ Parse the SHC file header. """
    fields = line.split()
    header = {
        "degree_min": int(fields[0]),
        "degree_max": int(fields[1]),
        "ntime": int(fields[2]),
        "spline_order": int(fields[3]),
        "nstep": int(fields[4]),
    }
    if fields[5:7]:
        header["validity_start"] = float(fields[5])
        header["validity_end"] = float(fields[6])
    return header


def _strip_shc_comments(file_in):
    """ Strip initial comments and empty lines from a text file stream. """
    for line in file_in:
        line = line.partition("#")[0].strip()
        if line:
            yield line
            break
    for line in file_in:
        line = line.strip()
        yield line
