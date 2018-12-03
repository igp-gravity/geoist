"""
GravMag: 3D forward modeling of total-field magnetic anomaly using rectangular
prisms (model with induced and remanent magnetization)
"""
import matplotlib.pyplot as plt
from geoist import gridder
from geoist.inversion import geometry
from geoist.pfm import prism, giutils
from geoist.vis import giplt

# The regional field
inc, dec = 30, -15
bounds = [-5000, 5000, -5000, 5000, 0, 5000]
model = [
    geometry.Prism(-4000, -3000, -4000, -3000, 0, 2000,
                 {'magnetization': giutils.ang2vec(1, inc, dec)}),
    geometry.Prism(-1000, 1000, -1000, 1000, 0, 2000,
                 {'magnetization': giutils.ang2vec(1, inc, dec)}),
    # This prism will have magnetization in a different direction
    geometry.Prism(2000, 4000, 3000, 4000, 0, 2000,
                 {'magnetization': giutils.ang2vec(3, -10, 45)})]
# Create a regular grid at 100m height
shape = (200, 200)
area = bounds[:4]
xp, yp, zp = gridder.regular(area, shape, z=-500)
# Calculate the anomaly for a given regional field

tf = prism.tf(xp, yp, zp, model, inc, dec)
# Plot
plt.figure()
plt.title("Total-field anomaly (nT)")
plt.axis('scaled')
giplt.contourf(yp, xp, tf, shape, 15)
plt.colorbar()
plt.xlabel('East y (km)')
plt.ylabel('North x (km)')
giplt.m2km()
plt.show()
