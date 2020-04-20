"""
GravMag: Forward modeling of the gravitational potential and its derivatives
using 3D model
"""
# 3rd imports
import matplotlib.pyplot as plt
# local imports
from geoist import gridder
from geoist.inversion import geometry
from geoist.pfm import prism
from geoist.vis import giplt

model = [geometry.Prism(-4000, -3000, -4000, -3000, 0, 2000, {'density': 1000}),
         geometry.Prism(-1000, 1000, -1000, 1000, 0, 2000, {'density': -900}),
         geometry.Prism(2000, 4000, 3000, 4000, 0, 2000, {'density': 1300})]
shape = (100, 100)
xp, yp, zp = gridder.regular((-5000, 5000, -5000, 5000), shape, z=-150)

field0 = prism.potential(xp, yp, zp, model)

from geoist.pfm import giutils
field0 = giutils.contaminate(field0, 0.05, percent = True)

fields = [prism.gx(xp, yp, zp, model),
          prism.gy(xp, yp, zp, model),
          prism.gz(xp, yp, zp, model),
          prism.gxx(xp, yp, zp, model),
          prism.gxy(xp, yp, zp, model),
          prism.gxz(xp, yp, zp, model),
          prism.gyy(xp, yp, zp, model),
          prism.gyz(xp, yp, zp, model),
          prism.gzz(xp, yp, zp, model)]
titles = ['potential', 'gx', 'gy', 'gz',
          'gxx', 'gxy', 'gxz', 'gyy', 'gyz', 'gzz']
plt.figure(figsize=(8, 8))
plt.axis('scaled')
plt.title(titles[0])
levels = giplt.contourf(yp * 0.001, xp * 0.001, field0, shape, 15)
cb = plt.colorbar()
giplt.contour(yp * 0.001, xp * 0.001, field0, shape,
            levels, clabel=False, linewidth=0.1)
plt.show()

plt.figure(figsize=(8, 8))
plt.subplots_adjust(left=0.03, right=0.95, bottom=0.05, top=0.92, hspace=0.3)
plt.suptitle("Potential fields produced by a 3 prism model")
for i, field in enumerate(fields):
    plt.subplot(3, 3, i + 1)
    plt.axis('scaled')
    plt.title(titles[i+1])
    levels = giplt.contourf(yp * 0.001, xp * 0.001, field, shape, 15)
    cb = plt.colorbar()
    giplt.contour(yp * 0.001, xp * 0.001, field, shape,
                levels, clabel=False, linewidth=0.1)
plt.show()
