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
This :mod:`~plateflex` module contains the following functions: 

- :func:`~plateflex.estimate.bayes_estimate_cell`: Set up :mod:`~pymc` model and estimate the parameters
  of the elastic plate model using a probabilistic Bayesian inference method.  
- :func:`~plateflex.estimate.get_bayes_estimates`: Explore the output of sampling the :mod:`~pymc` model
- :func:`~plateflex.estimate.L2_estimate_cell`: Set up non-linear curve fitting to estimate the parameters
  of the elastic plate model using non-linear least-squares from the function :func:`scipy.optimize.curve_fit`.  
- :func:`~plateflex.estimate.get_L2_estimates`: Explore the output the non-linear inversion
- :func:`~plateflex.estimate.real_xspec_functions`: Calculate the analytical admittance and coherence functions. 

Internal functions are also available to define predicted admittance and coherence data
with ``theano`` decorators to be incorporated as pymc variables. These functions are
used within :class:`~plateflex.classes.Project` methods as with :mod:`~plateflex.plotting`
functions.

.. warning::

    If you plan to estimate model parameters over entire grids, the non-linear least-squares method
    is orders of magnitude faster than the probabilistic method and should be preferred. The probabilistic
    method, on the other hand, gives useful estimation statistics from sampling the posterior distribution
    and can provide better insight on the trade-off between parameters. Estimates from the two methods are
    otherwise asymptotically identical (i.e., given infinite sampling of the posterior)

"""

# -*- coding: utf-8 -*-
import numpy as np
import pymc3 as pm
from geoist.flex import flex
from geoist.flex import conf as cf
from theano.compile.ops import as_op
import theano.tensor as tt
from scipy.optimize import curve_fit
import pandas as pd


def bayes_estimate_cell(k, adm, eadm, coh, ecoh, alph=False, atype='joint'):
    """
    Function to estimate the parameters of the flexural model at a single cell location
    of the input grids. 

    :type k: :class:`~numpy.ndarray`
    :param k: 1D array of wavenumbers
    :type adm: :class:`~numpy.ndarray`
    :param adm: 1D array of wavelet admittance
    :type eadm: :class:`~numpy.ndarray`
    :param eadm: 1D array of error on wavelet admittance
    :type coh: :class:`~numpy.ndarray`
    :param coh: 1D array of wavelet coherence
    :type ecoh: :class:`~numpy.ndarray`
    :param ecoh: 1D array of error on wavelet coherence
    :type alph: bool, optional
    :param alph: Whether or not to estimate parameter ``alpha``
    :type atype: str, optional
    :param atype: Whether to use the admittance (`'admit'`), coherence (`'coh'`) or both (`'joint'`)

    :return:
        (tuple): Tuple containing:
            * ``trace`` : :class:`~pymc3.backends.base.MultiTrace`
                Posterior samples from the MCMC chains
            * ``summary`` : :class:`~pandas.core.frame.DataFrame`
                Summary statistics from Posterior distributions
            * ``map_estimate`` : dict
                Container for Maximum a Posteriori (MAP) estimates

    """

    with pm.Model() as model:

        # k is an array - needs to be passed as distribution
        k_obs = pm.Normal('k', mu=k, sigma=1., observed=k)

        # Prior distributions
        Te = pm.Uniform('Te', lower=1., upper=250.)
        F = pm.Uniform('F', lower=0., upper=0.9999)

        if alph:

            # Prior distribution of `alpha`
            alpha = pm.Uniform('alpha', lower=0., upper=np.pi)
            admit_exp, coh_exp = real_xspec_functions_alpha(k_obs, Te, F, alpha)

        else:
            admit_exp, coh_exp = real_xspec_functions_noalpha(k_obs, Te, F)

        # Select type of analysis to perform
        if atype=='admit':

            # Uncertainty as observed distribution
            sigma = pm.Normal('sigma', mu=eadm, sigma=1., \
                observed=eadm)

            # Likelihood of observations
            admit_obs = pm.Normal('admit_obs', mu=admit_exp, \
                sigma=sigma, observed=adm)

        elif atype=='coh':

            # Uncertainty as observed distribution
            sigma = pm.Normal('sigma', mu=ecoh, sigma=1., \
                observed=ecoh)

            # Likelihood of observations
            coh_obs = pm.Normal('coh_obs', mu=coh_exp, \
                sigma=sigma, observed=coh)

        elif atype=='joint':

            # Define uncertainty as concatenated arrays
            ejoint = np.array([eadm, ecoh]).flatten()

            # Define array of observations and expected values as concatenated arrays
            joint = np.array([adm, coh]).flatten()
            joint_exp = tt.flatten(tt.concatenate([admit_exp, coh_exp]))

            # Uncertainty as observed distribution
            sigma = pm.Normal('sigma', mu=ejoint, sigma=1., \
                observed=ejoint)

            # Likelihood of observations
            joint_obs = pm.Normal('admit_coh_obs', mu=joint_exp, \
                sigma=sigma, observed=joint)

        # Sample the Posterior distribution
        trace = pm.sample(cf.draws, tune=cf.tunes, cores=cf.cores)

        # Get Max a porteriori estimate
        map_estimate = pm.find_MAP()

        # Get Summary
        summary = pm.summary(trace)

    return trace, summary, map_estimate

def get_bayes_estimates(summary, map_estimate):
    """
    Returns digestible estimates from the Posterior distributions.

    :type summary: :class:`~pandas.core.frame.DataFrame`
    :param summary: Summary statistics from Posterior distributions
    :type map_estimate: dict
    :param map_estimate: Container for Maximum a Posteriori (MAP) estimates

    :return: 
        (tuple): tuple containing:
            * mean_te (float) : Mean value of elastic thickness ``Te`` from posterior (km)
            * std_te (float)  : Standard deviation of elastic thickness ``Te`` from posterior (km)
            * C2_5_te (float) : Lower limit of 95% confidence interval on ``Te`` (km)
            * C97_5_te (float) : Upper limit of 95% confidence interval on ``Te`` (km)
            * MAP_te (float) : Maximum a Posteriori ``Te`` (km)
            * mean_F (float)  : Mean value of load ratio ``F`` from posterior
            * std_F (float)   : Standard deviation of load ratio ``F`` from posterior
            * C2_5_F (float) : Lower limit of 95% confidence interval on ``F``
            * C97_5_F (float) : Upper limit of 95% confidence interval on ``F``
            * MAP_F (float)  : Maximum a Posteriori load ratio ``F``
            * mean_a (float, optional)  : Mean value of initial phase difference ``alpha`` from posterior
            * std_a (float, optional)   : Standard deviation of initial phase difference `alpha`` from posterior
            * C2_5_a (float, optional) : Lower limit of 95% confidence interval on ``alpha``
            * C97_5_a (float, optional) : Upper limit of 95% confidence interval on ``alpha``
            * MAP_a (float, optional)  : Maximum a Posteriori initial phase difference ``alpha``

    """

    mean_a = None

    # Go through all estimates
    for index, row in summary.iterrows():
        if index=='Te':
            mean_te = row['mean']
            std_te = row['sd']
            C2_5_te = row['hpd_2.5']
            C97_5_te = row['hpd_97.5']
            MAP_te = np.float(map_estimate['Te'])
        elif index=='F':
            mean_F = row['mean']
            std_F = row['sd']
            C2_5_F = row['hpd_2.5']
            C97_5_F = row['hpd_97.5']
            MAP_F = np.float(map_estimate['F'])
        elif index=='alpha':
            mean_a = row['mean']
            std_a = row['sd']
            C2_5_a = row['hpd_2.5']
            C97_5_a = row['hpd_97.5']
            MAP_a = np.float(map_estimate['alpha'])

    if mean_a is not None:
        return mean_te, std_te, C2_5_te, C97_5_te, MAP_te, \
            mean_F, std_F, C2_5_F, C97_5_F, MAP_F, \
            mean_a, std_a, C2_5_a, C97_5_a, MAP_a
    else:
        return mean_te, std_te, C2_5_te, C97_5_te, MAP_te, \
            mean_F, std_F, C2_5_F, C97_5_F, MAP_F


def L2_estimate_cell(k, adm, eadm, coh, ecoh, alph=False, atype='joint'):
    """
    Function to estimate the parameters of the flexural model at a single cell location
    of the input grids. 

    :type k: :class:`~numpy.ndarray`
    :param k: 1D array of wavenumbers
    :type adm: :class:`~numpy.ndarray`
    :param adm: 1D array of wavelet admittance
    :type eadm: :class:`~numpy.ndarray`
    :param eadm: 1D array of error on wavelet admittance
    :type coh: :class:`~numpy.ndarray`
    :param coh: 1D array of wavelet coherence
    :type ecoh: :class:`~numpy.ndarray`
    :param ecoh: 1D array of error on wavelet coherence
    :type alph: bool, optional
    :param alph: Whether or not to estimate parameter ``alpha``
    :type atype: str, optional
    :param atype: Whether to use the admittance (`'admit'`), coherence (`'coh'`) or both (`'joint'`)

    :return:
        (tuple): Tuple containing:
            * ``summary`` : :class:`~pandas.core.frame.DataFrame`
                Summary statistics from L2 estimation

    """
    
    def pred_admit(k, Te, F, alpha):

        return flex.real_xspec_functions(k, Te, F, alpha)[0]

    def pred_coh(k, Te, F, alpha):

        return flex.real_xspec_functions(k, Te, F, alpha)[1]

    def pred_joint(k, Te, F, alpha):

        admittance, coherence = flex.real_xspec_functions(k, Te, F, alpha)
        return np.array([admittance, coherence]).flatten()

    if atype=='admit':
        y_obs = adm
        y_err = eadm
        if alph:
            theta0 = np.array([20., 0.5, np.pi/2.])
            function = lambda k, Te, F, alpha: pred_admit(k, Te, F, alpha)
            p1fit, p1cov = curve_fit(function, k, y_obs, p0=theta0, \
                sigma=y_err, absolute_sigma=True, max_nfev=1000, \
                bounds=([2., 0.0001, 0.0001], [250., 0.9999, np.pi-0.001]))

            # Calculate best fit function
            pred = pred_admit(k, p1fit[0], p1fit[1], p1fit[2])
            
            # calculate reduced chi-square
            rchi2 = np.sum((pred - y_obs)**2\
                    /y_err**2)/(len(pred)-len(p1fit))

        else:
            theta0 = np.array([20., 0.5])
            function = lambda k, Te, F: pred_admit(k, Te, F, alpha=np.pi/2.)
            p1fit, p1cov = curve_fit(function, k, y_obs, p0=theta0, \
                sigma=y_err, absolute_sigma=True, max_nfev=1000, \
                bounds=([2., 0.0001], [250., 0.9999]))

            # Calculate best fit function
            pred = pred_admit(k, p1fit[0], p1fit[1], np.pi/2.)
            
            # calculate reduced chi-square
            rchi2 = np.sum((pred - y_obs)**2\
                    /y_err**2)/(len(pred)-len(p1fit))

    elif atype=='coh':
        y_obs = coh
        y_err = ecoh
        if alph:
            theta0 = np.array([20., 0.5, np.pi/2.])
            function = lambda k, Te, F, alpha: pred_coh(k, Te, F, alpha)
            p1fit, p1cov = curve_fit(function, k, y_obs, p0=theta0, \
                sigma=y_err, absolute_sigma=True, max_nfev=1000, \
                bounds=([2., 0.0001, 0.0001], [250., 0.9999, np.pi-0.001]))

            # Calculate best fit function
            pred = pred_coh(k, p1fit[0], p1fit[1], p1fit[2])
            
            # calculate reduced chi-square
            rchi2 = np.sum((pred - y_obs)**2\
                    /y_err**2)/(len(pred)-len(p1fit))

        else:
            theta0 = np.array([20., 0.5])
            function = lambda k, Te, F: pred_coh(k, Te, F, alpha=np.pi/2.)
            p1fit, p1cov = curve_fit(function, k, y_obs, p0=theta0, \
                sigma=y_err, absolute_sigma=True, max_nfev=1000, \
                bounds=([2., 0.0001], [250., 0.9999]))

            # Calculate best fit function
            pred = pred_coh(k, p1fit[0], p1fit[1], np.pi/2.)
            
            # calculate reduced chi-square
            rchi2 = np.sum((pred - y_obs)**2\
                    /y_err**2)/(len(pred)-len(p1fit))

    elif atype=='joint':
        y_obs = np.array([adm, coh]).flatten()
        y_err = np.array([eadm, ecoh]).flatten()
        if alph:
            theta0 = np.array([20., 0.5, np.pi/2.])
            function = lambda k, Te, F, alpha: pred_joint(k, Te, F, alpha)
            p1fit, p1cov = curve_fit(function, k, y_obs, p0=theta0, \
                sigma=y_err, absolute_sigma=True, max_nfev=1000, \
                bounds=([2., 0.0001, 0.0001], [250., 0.9999, np.pi-0.001]))

            # Calculate best fit function
            pred = pred_joint(k, p1fit[0], p1fit[1], p1fit[2])
            
            # calculate reduced chi-square
            rchi2 = np.sum((pred - y_obs)**2\
                    /y_err**2)/(len(pred)-len(p1fit))

        else:
            theta0 = np.array([20., 0.5])
            function = lambda k, Te, F: pred_joint(k, Te, F, alpha=np.pi/2.)
            p1fit, p1cov = curve_fit(function, k, y_obs, p0=theta0, \
                sigma=y_err, absolute_sigma=True, max_nfev=1000, \
                bounds=([2., 0.0001], [250., 0.9999]))

            # Calculate best fit function
            pred = pred_joint(k, p1fit[0], p1fit[1], np.pi/2.)
            
            # calculate reduced chi-square
            rchi2 = np.sum((pred - y_obs)**2\
                    /y_err**2)/(len(pred)-len(p1fit))

    p1err = np.sqrt(np.diag(p1cov))

    if alph:

        # Store summary
        summary = pd.DataFrame(data={'mean':[p1fit[0], p1fit[1], p1fit[2]], \
            'std':[p1err[0], p1err[1], p1err[2]], 'chi2':[rchi2, rchi2, rchi2]}, \
            index=['Te', 'F', 'alpha'])

    else:

        summary = pd.DataFrame(data={'mean':[p1fit[0], p1fit[1]], \
            'std':[p1err[0], p1err[1]], 'chi2':[rchi2, rchi2]}, \
            index=['Te', 'F'])

    return summary


def get_L2_estimates(summary):
    """
    Returns digestible estimates from the L2 estimates.

    :type summary: :class:`~pandas.core.frame.DataFrame`
    :param summary: Summary statistics from Posterior distributions

    :return: 
        (tuple): tuple containing:
            * mean_te (float) : Mean value of elastic thickness from posterior (km)
            * std_te (float)  : Standard deviation of elastic thickness from posterior (km)
            * mean_F (float)  : Mean value of load ratio from posterior
            * std_F (float)   : Standard deviation of load ratio from posterior
            * mean_a (float, optional)  : Mean value of phase difference between initial loads from posterior
            * std_a (float, optional)   : Standard deviation of phase difference between initial loads from posterior
            * rchi2 (float)   : Reduced chi-squared value
    """

    mean_a = None

    # Go through all estimates
    for index, row in summary.iterrows():
        if index=='Te':
            mean_te = row['mean']
            std_te = row['std']
            rchi2 = row['chi2']
        elif index=='F':
            mean_F = row['mean']
            std_F = row['std']
        elif index=='alpha':
            mean_a = row['mean']
            std_a = row['std']

    if mean_a is not None:
        return mean_te, std_te, mean_F, std_F, mean_a, std_a, rchi2
    else:
        return mean_te, std_te, mean_F, std_F, rchi2


def real_xspec_functions(k, Te, F, alpha=np.pi/2.):
    """
    Calculate analytical expressions for the real component of admittance
    and coherence functions. 

    :type k: np.ndarray
    :param k: Wavenumbers (rad/m)
    :type Te: float
    :param Te: Effective elastic thickness (km)
    :type F: float
    :param F: Subsurface-to-surface load ratio [0, 1[
    :type alpha: float, optional
    :param alpha: Phase difference between initial applied loads (rad)

    :return:  
        (tuple): tuple containing:
            * admittance (:class:`~numpy.ndarray`): Real admittance function (shape: ``len(k)``)
            * coherence (:class:`~numpy.ndarray`): Coherence functions (shape: ``len(k)``)

    """

    admittance, coherence = flex.real_xspec_functions(k, Te, F, alpha)

    return admittance, coherence


@as_op(itypes=[tt.dvector, tt.dscalar, tt.dscalar], 
    otypes=[tt.dvector, tt.dvector])
def real_xspec_functions_noalpha(k, Te, F):
    """
    Calculate analytical expressions for the real component of admittance, 
    coherency and coherence functions. 
    """

    adm, coh = real_xspec_functions(k, Te, F)

    return adm, coh

@as_op(itypes=[tt.dvector, tt.dscalar, tt.dscalar, tt.dscalar], 
    otypes=[tt.dvector, tt.dvector])
def real_xspec_functions_alpha(k, Te, F, alpha):
    """
    Calculate analytical expressions for the real component of admittance, 
    coherency and coherence functions. 
    """

    adm, coh = real_xspec_functions(k, Te, F, alpha)

    return adm, coh


