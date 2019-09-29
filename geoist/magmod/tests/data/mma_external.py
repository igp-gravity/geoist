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

DEGREE = 2
TIME = 2016.99965847

DATA = [
    (1, 0, 30.7117594582906754),
    (1, 1, 0.741232780798051327),
    (1, -1, -8.95006273357913607),
    (2, 0, 0.968979330373007097),
    (2, 1, 1.42731572806433937),
    (2, -1, 7.0650510135243545),
]

COEF_Q = zeros(((DEGREE + 1)*(DEGREE + 2))//2)
COEF_S = zeros(((DEGREE + 1)*(DEGREE + 2))//2)

for n, m, coeff in DATA:
    if m >= 0:
        COEF_Q[(n*(n + 1))//2 + m] = coeff
    else:
        COEF_S[(n*(n + 1))//2 - m] = coeff
