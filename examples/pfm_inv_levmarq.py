# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 12:41:42 2020

@author: chens

GravMag: 2D gravity inversion for the relief of a basin
"""
from geoist.inversion import Smoothness1D
from geoist.pfm.basin2d import PolygonalBasinGravity
from geoist.pfm import talwani
from geoist.inversion.geometry import Polygon
from geoist.vis import giplt 
from geoist.pfm import giutils as utils
import numpy as np
import matplotlib.pyplot as plt

# Make some synthetic data to test the inversion
# The model will be a polygon.
# Reverse x because vertices must be clockwise.
xs = np.linspace(0, 100000, 100)[::-1]
depths = (-1e-15*(xs - 50000)**4 + 8000 -
          3000*np.exp(-(xs - 70000)**2/(10000**2)))
depths -= depths.min()  # Reduce depths to zero
props = {'density': -300}
model = Polygon(np.transpose([xs, depths]), props)
x = np.linspace(0, 100000, 100)
z = -100*np.ones_like(x)
data = utils.contaminate(talwani.gz(x, z, [model]), 0.5, seed=0)

# Make the solver using smoothness regularization and run the inversion
misfit = PolygonalBasinGravity(x, z, data, 50, props, top=0)
regul = Smoothness1D(misfit.nparams)
solver = misfit + 1e-4*regul
# This is a non-linear problem so we need to pick an initial estimate
initial = 3000*np.ones(misfit.nparams)
solver.config('levmarq', initial=initial).fit()
plt.figure()
plt.subplot(2, 1, 1)
plt.plot(x, data, 'ok', label='observed')
plt.plot(x, solver[0].predicted(), '-r', linewidth=2, label='predicted')
plt.legend()
ax = plt.subplot(2, 1, 2)
giplt.polygon(model, fill='gray', alpha=0.5, label='True')
# The estimate_ property of our solver gives us the estimate basin as a polygon
# So we can directly pass it to plotting and forward modeling functions
giplt.polygon(solver.estimate_, style='o-r', label='Estimated')
ax.invert_yaxis()
plt.legend()
plt.show()
