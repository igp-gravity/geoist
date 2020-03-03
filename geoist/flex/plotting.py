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
This :mod:`~plateflex`` module contains the following functions for plotting:

- :func:`~plateflex.plotting.plot_real_grid`
- :func:`~plateflex.plotting.plot_bayes_stats`
- :func:`~plateflex.plotting.plot_functions`

"""

# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde as kde
import seaborn as sns
sns.set()


def plot_real_grid(grid, log=False, mask=None, title=None, save=None, clabel=None, contours=None, **kwargs):
    """
    Plot 2D image of any real-valued 2D array, used in several context throughout
    :mod:`~plateflex`. For example, it can be used to plot the input grids of topography
    or gravity anomalies, the real or imaginary values of the wavelet transform at a
    giben wavenumber index, the wavelet scalograms at a given wavenumber index, the
    wavelet admittance or coherence at a given wavenumber index, or the final grids of
    results.

    :type grid: :class:`~numpy.ndarray` 
    :param grid: Array of real-valued data
    :type log: bool, optional 
    :param log: Whether or not to take the log of the array values (useful in scalogram)
    :type mask: np.ndarray, optional 
    :param mask: Array of booleans for masking data points 
    :type title: str, optional
    :param title: Title of plot
    :type save: str, optional
    :param save: Name of file for to save figure
    :type clabel: str, optional
    :param clabel: Label for colorbar
    :type contours: List
    :param contours: Contours to overlay on maps (e.g., useful for plotting outline of land areas)

    """

    # Take log of real values
    if log:
        if not np.all(grid==np.absolute(grid)):
            raise(Exception('cannot plot log of grid containing \
                negative values'))
        grid = np.log(grid)

    # Apply mask
    if mask is not None:
        grid = np.ma.masked_where(mask, grid)

    # Plot figure and add colorbar
    plt.figure()
    plt.imshow(grid, origin='lower', **kwargs)
    cbar = plt.colorbar()

    # Add units on colorbar label
    if clabel is not None:
        cbar.set_label(clabel)

    # Add contours
    if contours is not None:
        try: 
            for n, contour in enumerate(contours):
                plt.plot(contour[:,1], contour[:,0], lw=1.25, c='k')
        except:
            print("No contours exist for map. Passing.")

    # Plot title if requested
    if title is not None:
        plt.title(title)
    
    # Save figure
    if save is not None:
        plt.savefig(save+'.png')

    # Show
    plt.show()


def plot_bayes_stats(trace, summary, map_estimate, title=None, save=None):
    """
    Extract results from variables ``trace``, ``summary`` and ``map_estimate`` to 
    plot marginal and joint posterior distributions. Automatically determines
    how to plot results from those variables.

    :type trace: :class:`~pymc3.backends.base.MultiTrace`
    :param trace: Posterior samples from the MCMC chains
    :type summary: :class:`~pandas.core.frame.DataFrame`
    :param summary: Summary statistics from Posterior distributions
    :type map_estimate: dict
    :param map_estimate: Container for Maximum a Posteriori (MAP) estimates
    :type title: str, optional 
    :param title: Title of plot
    :type save: str, optional
    :param save: Name of file for to save figure

    """

    from plateflex import estimate

    # Extract results from summary and map_estimate
    results = estimate.get_bayes_estimates(summary, map_estimate)

    # Collect keys in trace object
    keys = []
    for var in trace.varnames:
        if var[-1]=='_':
            continue
        keys.append(var)

    # This means we searched for Te and F only
    if len(keys)==2:

        # Collect pymc chains as ``pandas.DataFrame`` object
        data = np.array([trace['Te'], trace['F']]).transpose()
        data = pd.DataFrame(data, columns=['Te (km)', 'F'])

        # Plot marginal and joint distributions as histograms and kernel density functions
        g = sns.PairGrid(data)
        g.map_diag(plt.hist, lw=1)
        g.map_lower(sns.kdeplot)

        # Set unused plot axes to invisible
        ax = g.axes[0][1]
        ax.set_visible(False)

        # Text for Te statistics
        tetext = '\n'.join((
            r'$\mu$ = {0:.0f} km'.format(results[0]),
            r'$\sigma$ = {0:.0f} km'.format(results[1]),
            r'$95\%$ CI = [{0:.0f}, {1:.0f}] km'.format(results[2], results[3]),
            r'MAP = {0:.0f} km'.format(results[4])))

        # Insert text as box
        ax1 = g.axes[0][0]
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax1.text(1.05, 0.9, tetext, transform=ax1.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)

        # Text for F statistics
        Ftext = '\n'.join((
            r'$\mu$ = {0:.2f}'.format(results[5]),
            r'$\sigma$ = {0:.2f}'.format(results[6]),
            r'$95\%$ CI = [{0:.2f}, {1:.2f}]'.format(results[7], results[8]),
            r'MAP = {0:.2f}'.format(results[9])))

        # Insert text as box
        ax2 = g.axes[1][1]
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax2.text(0.135, 1.4, Ftext, transform=ax2.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)

    # This means we searched for Te, F and alpha
    elif len(keys)==3:

        # Collect pymc chains as ``pandas.DataFrame`` object
        data = np.array([trace['Te'], trace['F'], trace['alpha']]).transpose()
        data = pd.DataFrame(data, columns=['Te (km)', 'F', r'$\alpha$'])

        # Plot marginal and joint distributions as histograms and kernel density functions
        g = sns.PairGrid(data)
        g.map_diag(plt.hist, lw=1)
        g.map_lower(sns.kdeplot)

        # Set unused plot axes to invisible
        ax = g.axes[0][1]
        ax.set_visible(False)
        ax = g.axes[0][2]
        ax.set_visible(False)
        ax = g.axes[1][2]
        ax.set_visible(False)

        # Text for Te statistics
        tetext = '\n'.join((
            r'$\mu$ = {0:.0f} km'.format(results[0]),
            r'$\sigma$ = {0:.0f} km'.format(results[1]),
            r'$95\%$ CI = [{0:.0f}, {1:.0f}] km'.format(results[2], results[3]),
            r'MAP = {0:.0f} km'.format(results[4])))

        # Insert text as box
        ax1 = g.axes[0][0]
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax1.text(1.05, 0.9, tetext, transform=ax1.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)

        # Text for F statistics
        Ftext = '\n'.join((
            r'$\mu$ = {0:.2f}'.format(results[5]),
            r'$\sigma$ = {0:.2f}'.format(results[6]),
            r'$95\%$ CI = [{0:.2f}, {1:.2f}]'.format(results[7], results[8]),
            r'MAP = {0:.2f}'.format(results[9])))

        # Insert text as box
        ax2 = g.axes[1][1]
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax2.text(0.135, 1.4, Ftext, transform=ax2.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)

        # Text for alpha statistics
        atext = '\n'.join((
            r'$\mu$ = {0:.2f}'.format(results[10]),
            r'$\sigma$ = {0:.2f}'.format(results[11]),
            r'$95\%$ CI = [{0:.2f}, {1:.2f}]'.format(results[12], results[13]),
            r'MAP = {0:.2f}'.format(results[14])))

        # Insert text as box
        ax3 = g.axes[2][2]
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax3.text(0.135, 1.4, atext, transform=ax3.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)

    else:
        raise(Exception('there are less than 2 or more than 3 variables in pymc3 chains'))

    # Plot title if requested
    if title is not None:
        plt.suptitle(title)

    # Save figure
    if save is not None:
        plt.savefig(save+'.png')

    plt.show()


def plot_functions(k, adm, eadm, coh, ecoh, padm=None, pcoh=None, title=None, save=None):
    """
    Function to plot observed and predicted (``None`` by default) admittance and coherence functions. 
    Both admittance and coherence are plotted regardless of method to estimate 
    the model paramters.

    :type k: :class:`~numpy.ndarray`
    :param k: 1D array of wavenumbers
    :type adm: :class:`~numpy.ndarray`
    :param adm: 1D array of observed wavelet admittance
    :type eadm: :class:`~numpy.ndarray`
    :param eadm: 1D array of error on observed wavelet admittance
    :type coh: :class:`~numpy.ndarray`
    :param coh: 1D array of observed wavelet coherence
    :type ecoh: :class:`~numpy.ndarray`
    :param ecoh: 1D array of error on observed wavelet coherence
    :type padm: :class:`~numpy.ndarray`
    :param padm: 1D array of predicted wavelet admittance
    :type pcoh: :class:`~numpy.ndarray`
    :param pcoh: 1D array of predicted wavelet coherence
    :type title: str, optional 
    :param title: Title of plot
    :type save: str, optional
    :param save: Name of file for to save figure

    """

    from plateflex import estimate

    # Plot as 2 subplots
    f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

    # Plot observed admittance with error bars
    ax1.errorbar(k*1.e3,adm,yerr=eadm, marker='*')

    if padm is not None:
        # Plot predicted admittance
        ax1.plot(k*1.e3,padm)

    # Plot observed coherence with error bars
    ax2.errorbar(k*1.e3,coh,yerr=ecoh, marker='*')

    if pcoh is not None:
        # Plot predicted coherence
        ax2.plot(k*1.e3,pcoh)

    # Add all labels
    ax1.set_ylabel('Admittance (mGal/m)')
    ax1.set_xscale('log')

    ax2.set_ylabel('Coherence')
    ax2.set_xlabel('Wavenumber (rad/km)')
    ax2.set_xscale('log')

    # Plot title if requested
    if title is not None:
        plt.suptitle(title)

    # Save figure
    if save is not None:
        plt.savefig(save+'.png')

    plt.show()
