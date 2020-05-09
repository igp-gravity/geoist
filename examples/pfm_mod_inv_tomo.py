# -*- coding: utf-8 -*-
"""
Created on Fri May  8 12:35:57 2020

@author: chens
"""

import numpy as np
from geoist.inversion.mesh import SquareMesh
from geoist.pfm.inv3d import SRTomo
from geoist.inversion import Smoothness2D, LCurve, ttime2d
from geoist.pfm import giutils
from geoist.gridder import scatter, circular_scatter

area = (0, 2, 0, 2)
shape = (10, 10)
model = SquareMesh(area, shape)
vp = 4*np.ones(shape)
vp[3:7,3:7] = 10
print(vp)

model.addprop('vp', vp.ravel())
src_loc_x, src_loc_y = scatter(area, 30, seed=0)
src_loc = np.transpose([src_loc_x, src_loc_y])
rec_loc_x, rec_loc_y = circular_scatter(area, 20, random=True, seed=0)
rec_loc = np.transpose([rec_loc_x, rec_loc_y])
srcs = [src for src in src_loc for _ in rec_loc]
recs = [rec for _ in src_loc for rec in rec_loc]
tts = ttime2d.straight(model, 'vp', srcs, recs)
tts = giutils.contaminate(tts, 0.01, percent=True, seed=0)

# Now we can setup a tomography by creating the necessary data misfit
# (``SRTomo``) and regularization (``Smoothness2D``) objects. We'll normalize
# the data misfit by the number of data points to make the scale of the
# regularization parameter more tractable.

mesh = SquareMesh(area, shape)
datamisfit = (1./tts.size)*SRTomo(tts, srcs, recs, mesh)
regul = Smoothness2D(mesh.shape)

# The tomography solver will be the ``LCurve`` solver. It works by calling
# ``fit()`` and accessing ``estimate_``, exactly like any other solver:

regul_params = [10**i for i in range(-10, -2, 1)]
tomo = LCurve(datamisfit, regul, regul_params)
_ = tomo.fit()
print(np.array_repr(tomo.estimate_.reshape(shape), precision=0))

# When ``fit()`` is called, the ``LCurve``  will run the inversion for each
# value of the regularization parameter, build an l-curve, and find the
# best solution (i.e., the corner value of the l-curve).

# The ``LCurve`` object behaves like a normal multi-objective function.
# In fact, it will try to mirror the objective function that resulted in the
# best solution.
# You can index it to access the data misfit and regularization parts.
# For example, to get the residuals vector or the predicted data:

predicted = tomo[0].predicted()
residuals = tomo[0].residuals()
print('%.4f %.4f' % (residuals.mean(), residuals.std()))

# The estimated regularization parameter is stored in ``regul_param_``:

print('Hyperparameter Lambda value is {}'.format(tomo.regul_param_))
tomo.plot_lcurve()

# You can view the optimization information for the run corresponding to the
# best estimate using the ``stats_`` attribute:

list(sorted(tomo.stats_))
print(tomo.stats_['method'])
