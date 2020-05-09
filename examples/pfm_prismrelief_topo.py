"""
Meshing: Generate a 3D prism model of the topography
"""
from geoist import gridder
from geoist.pfm import giutils
from geoist.inversion import mesh

area = (-150, 150, -300, 300)
shape = (100, 50)
x, y = gridder.regular(area, shape)
height = (-80 * giutils.gaussian2d(x, y, 100, 200, x0=-50, y0=-100, angle=-60) +
          100 * giutils.gaussian2d(x, y, 50, 100, x0=80, y0=170))

nodes = (x, y, -1 * height)  # -1 is to convert height to z coordinate
reference = 0  # z coordinate of the reference surface
relief = mesh.PrismRelief(reference, shape, nodes)
relief.addprop('density', (2670 for i in range(relief.size)))

prop = 'density'
for prism in relief:
    if prism is None or (prop is not None and prop not in prism.props):
        continue
    x1, x2, y1, y2, z1, z2 = prism.get_bounds()
    

# Tesseroid
w, e = -2, 2
s, n = -2, 2
bounds = (w, e, s, n, 500000, 0)

x, y = gridder.regular((w, e, s, n), (50, 50))
height = (250000 +
          -100000 * giutils.gaussian2d(x, y, 1, 5, x0=-1, y0=-1, angle=-60) +
          250000 * giutils.gaussian2d(x, y, 1, 1, x0=0.8, y0=1.7))

mesh = mesh.TesseroidMesh(bounds, (20, 50, 50))
mesh.addprop('density', (2670 for i in range(mesh.size)))
mesh.carvetopo(x, y, height)    

# from geoist.vis import myv
# myv.figure()
# myv.prisms(relief, prop='density', edges=False)
# axes = myv.axes(myv.outline())
# myv.wall_bottom(axes.axes.bounds, opacity=0.2)
# myv.wall_north(axes.axes.bounds)
# myv.show()
