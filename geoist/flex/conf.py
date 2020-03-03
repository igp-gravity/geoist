# Copyright 2019 Pascal Audet
#
# This file is part of PlateFlex.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""

plateflex.cpwt.conf_cpwt
------------------------

Global variable that controls the spatio-spectral
localization of the Morlet wavelet. Default value 
in parentheses.

.. rubric:: Internal wavenumber of Morlet wavelet

``k0``: float
    Internal Morlet wavenumber (5.336 or higher)

plateflex.flex.conf_flex
------------------------

Global variables that define model parameters for
the flexural model. Default values in parentheses.

.. rubric:: Earth parameters - fixed

``E`` : float
    Young's modulus (100 GPa)
``nu`` : float
    Poisson's ratio (0.25)
``g`` : float
    Gravitational acceleration (9.81 m/s^2)
``G`` : float)
    Gravitational constant (6.67e-11*1.e5 mGal)

.. rubric:: Earth parameters - variable

``zc`` : float
    Crustal thickness (35.e3 m)
``rhom``: float
    Uppermost mantle density (3200. kg/m^3)
``rhoc`` : float
    Crustal density (2700. kg/m^3)
``rhoa`` : float
    Air density (0. kg/m^3)
``rhow`` : float
    Water density (1030. kg/m^3)
``rhof`` : float
    Fluid density at topo/fluid interface (==rhoa *or* ==rhow, depending on water depth)
``wd`` : float
    Water depth (0.e3 m) -  Automatically determined from 
    :class:`~plateflex.classes.TopoGrid` object

plateflex.conf
--------------

Global variables that control sampling of the posterior
from MCMC chains using :class:`pymc3`. Default values in parentheses.

.. rubric:: Bayes sampling

``draws`` : int 
    Number of draws (i.e., samples) in single MCMC chain (500)
``tunes`` : int
    Number of tuning (i.e., burn-in) samples (500)
``cores`` : int
    Number of cores (i.e., parallel MCMC chains) (4)

"""

# bayes parameters
global samples, tunes
draws = 500
tunes = 500
cores = 4