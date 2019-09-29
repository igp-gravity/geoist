#-------------------------------------------------------------------------------
#
#  EMM format parser
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
from numpy import array, stack


RE_HEADER_LINE = re.compile(r'^-+$')
EMM_VALIDITY_PERIOD = 5.0 # years


def combine_emm_coefficients(data_static, data_secvar):
    """ Combine coefficient read from the separate EMM files. """
    return (
        combine_emm_variable(data_static, data_secvar),
        extract_emm_constant(data_static, data_secvar),
    )


def combine_emm_variable(data_static, data_secvar):
    """ Combine EMM variable coefficient. """
    epoch_static = data_static["epoch"]
    epoch_secvar = data_secvar["epoch"]
    nm_static, gh_static = data_static["nm"], data_static["gh"]
    nm_secvar, gh_secvar = data_secvar["nm"], data_secvar["gh"]
    size = nm_secvar.shape[0]

    if epoch_static != epoch_secvar or (nm_secvar != nm_static[:size]).any():
        ValueError("Incompatible sets of the static and variation coefficients!")

    return {
        "epoch": epoch_secvar,
        "degree_min": data_secvar["degree_min"],
        "degree_max": data_secvar["degree_max"],
        "t": array([epoch_secvar, epoch_secvar + EMM_VALIDITY_PERIOD]),
        "nm": nm_secvar,
        "gh": stack((
            gh_static[:size, 0],
            gh_static[:size, 0] + gh_secvar[:, 0]*EMM_VALIDITY_PERIOD
        ), axis=-1)
    }


def extract_emm_constant(data_static, data_secvar):
    """ Extract EMM constant coefficient. """
    epoch = data_secvar["epoch"]
    size = data_secvar["nm"].shape[0]
    nm_static, gh_static = data_static["nm"], data_static["gh"]
    nm_static, gh_static = nm_static[size:], gh_static[size:]
    nonzeros, = gh_static[:, 0].nonzero()
    nm_static, gh_static = nm_static[nonzeros], gh_static[nonzeros]
    return {
        "epoch": epoch,
        "degree_min": nm_static[:, 0].min(),
        "degree_max": nm_static[:, 0].max(),
        "t": array([epoch]),
        "nm": nm_static,
        "gh": gh_static,
    }


def parse_emm_file(file_in):
    """ Parse EMM file format and return a dictionary containing the parsed
    model data.
    """
    lines = file_in
    data = parse_emm_header(lines)
    data["nm"], data["gh"], data["t"] = parse_emm_coefficients(lines, data)
    data["degree_min"] = data["nm"][:, 0].min()
    data["degree_max"] = data["nm"][:, 0].max()
    return data


def parse_emm_header(lines):
    """ Parse EMM file header. """
    data = {}
    for line in lines:
        fields = line.split()
        if len(fields) == 1 and RE_HEADER_LINE.match(fields[0]):
            break
        if len(fields) == 2:
            if fields[0] == "Epoch:":
                data["epoch"] = float(fields[1])

    return data


def parse_emm_coefficients(lines, data):
    """ Parse EMM coefficients. """
    nm_index = []
    coeff = []
    for line in lines:
        fields = line.split()
        n_idx, m_idx = int(fields[0]), int(fields[1])
        coef_g, coef_h = float(fields[2]), float(fields[3])
        nm_index.append((n_idx, m_idx))
        coeff.append(coef_g)
        if m_idx > 0:
            nm_index.append((n_idx, -m_idx))
            coeff.append(coef_h)
    coeff = array(coeff).reshape((len(coeff), 1))
    return array(nm_index), coeff, array([data["epoch"]])
