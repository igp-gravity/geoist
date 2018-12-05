"""
GravMag: Use the DipoleMagDir class to estimate the magnetization direction
of dipoles with known centers
"""
# 3rd imports
import numpy as np
import matplotlib.pyplot as plt
# local imports
from geoist import gridder
from geoist.inversion import geometry
from geoist.pfm.giutils import ang2vec, contaminate
from geoist.pfm import sphere
from geoist.pfm.magdir import DipoleMagDir
from geoist.vis import giplt

# Make noise-corrupted synthetic data
inc, dec = -10.0, -15.0  # inclination and declination of the Geomagnetic Field
model = [geometry.Sphere(3000, 3000, 1000, 1000,
                       {'magnetization': ang2vec(6.0, -20.0, -10.0)}),
         geometry.Sphere(7000, 7000, 1000, 1000,
                       {'magnetization': ang2vec(10.0, 3.0, -67.0)})]
area = (0, 10000, 0, 10000)
x, y, z = gridder.scatter(area, 1000, z=-150, seed=0)
tf = contaminate(sphere.tf(x, y, z, model, inc, dec), 5.0, seed=0)

# Give the centers of the dipoles
centers = [[3000, 3000, 1000], [7000, 7000, 1000]]

# Estimate the magnetization vectors
solver = DipoleMagDir(x, y, z, tf, inc, dec, centers).fit()

# Print the estimated and true dipole monents, inclinations and declinations
print('Estimated magnetization (intensity, inclination, declination)')
for e in solver.estimate_:
    print(e)

# Plot the fit and the normalized histogram of the residuals
plt.figure(figsize=(14, 5))
plt.subplot(1, 2, 1)
plt.title("Total Field Anomaly (nT)", fontsize=14)
plt.axis('scaled')
nlevels = giplt.contour(y, x, tf, (50, 50), 15, interp=True, color='r',
                      label='Observed', linewidth=2.0)
giplt.contour(y, x, solver.predicted(), (50, 50), nlevels, interp=True,
            color='b', label='Predicted', style='dashed', linewidth=2.0)
plt.legend(loc='upper left', shadow=True, prop={'size': 13})
plt.xlabel('East y (m)', fontsize=14)
plt.ylabel('North x (m)', fontsize=14)
plt.subplot(1, 2, 2)
residuals_mean = np.mean(solver.residuals())
residuals_std = np.std(solver.residuals())
# Each residual is subtracted from the mean and the resulting
# difference is divided by the standard deviation
s = (solver.residuals() - residuals_mean) / residuals_std
plt.hist(s, bins=21, range=None, normed=True, weights=None,
         cumulative=False, bottom=None, histtype='bar', align='mid',
         orientation='vertical', rwidth=None, log=False,
         color=None, label=None)
plt.xlim(-4, 4)
plt.title("mean = %.3f    std = %.3f" % (residuals_mean, residuals_std),
          fontsize=14)
plt.ylabel("P(z)", fontsize=14)
plt.xlabel("z", fontsize=14)
plt.show()
