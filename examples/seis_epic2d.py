# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 09:19:19 2020

@author: chens

Seismic: 2D epicenter estimation assuming a homogeneous and flat Earth
"""

import numpy
from geoist.pfm import giutils
from geoist.inversion.geometry import Square
import matplotlib.pyplot as plt
from geoist.vis import giplt
from geoist.inversion import ttime2d
from geoist.inversion.epic2d import Homogeneous

# Make a velocity model to calculate traveltimes
area = (0, 10, 0, 10)
vp, vs = 2, 1
model = [Square(area, props={'vp': vp, 'vs': vs})]

src = (5, 5)
srcs = [src, src, src, src]
recs = [(1, 2),(3,6),(4,7),(2,8)]

#giutils.connect_points(src, rec_points)
ptime = ttime2d.straight(model, 'vp', srcs, recs)
stime = ttime2d.straight(model, 'vs', srcs, recs)
# Calculate the residual time (S - P) with added noise
traveltime, error = giutils.contaminate(stime - ptime, 0.05, percent=True,
                                      return_stddev=True)
solver = Homogeneous(traveltime, recs, vp, vs)
    
initial = (1, 1)    
estimate = solver.config('levmarq', initial=initial).fit().estimate_

plt.figure(figsize=(10, 4))
plt.subplot(1, 2, 1)
plt.title('Epicenter + %d recording stations' % (len(recs)))
plt.axis('scaled')
giplt.points(src, '*y', label="True")
giplt.points(recs, '^r', label="Stations")
giplt.points(initial, '*b', label="Initial")
giplt.points([estimate], '*g', label="Estimate")
giplt.set_area(area)
plt.legend(loc='lower right', shadow=True, numpoints=1, prop={'size': 12})
plt.xlabel("X")
plt.ylabel("Y")
ax = plt.subplot(1, 2, 2)
plt.title('Travel-time residuals + error bars')
s = numpy.arange(len(traveltime)) + 1
width = 0.3
plt.bar(s - width, traveltime, width, color='g', label="Observed",
        yerr=error)
plt.bar(s, solver.predicted(), width, color='r', label="Predicted")
ax.set_xticks(s)
plt.legend(loc='upper right', shadow=True, prop={'size': 12})
plt.xlabel("Station number")
plt.ylabel("Travel-time residual")
plt.show()