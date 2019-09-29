#-------------------------------------------------------------------------------
#
#  Spherical Harmonic Expansion - Geomagnetic Model - test dataset
#
#  External Spherical Harmonic Coefficients - SWARM MMA
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
# all copies of this Software or works derived from this Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#-------------------------------------------------------------------------------

from numpy import zeros

DEGREE = 3
TIME = 2016.99965847

DATA = [
    (1, 0, 6.37661120235457446),
    (1, 1, -0.295167002837779691),
    (1, -1, -3.5069934604844577),
    (2, 0, -0.666513682337298352),
    (2, 1, 1.11417106500435548),
    (2, -1, 2.65809245405974481),
    (2, 2, -0.413685898767694404),
    (2, -2, 0.638934950391676026),
    (3, 0, 2.26438999698106125),
    (3, 1, 0.931464083813372312),
    (3, -1, 0.473701177575287014),
    (3, 2, 0.382569838536052365),
    (3, -2, 0.822668017809991881),
    (3, 3, -0.233587111584525969),
    (3, -3, -0.0149047101935679896)
]

COEF_G = zeros(((DEGREE + 1)*(DEGREE + 2))//2)
COEF_H = zeros(((DEGREE + 1)*(DEGREE + 2))//2)

for n, m, coeff in DATA:
    if m >= 0:
        COEF_G[(n*(n + 1))//2 + m] = coeff
    else:
        COEF_H[(n*(n + 1))//2 - m] = coeff
