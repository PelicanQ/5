import numpy as np
import calfem.core as cfc
import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.utils as cfu
from constants import a,b,h,t,d,c,L
import matplotlib as mpl
import calfem.vis_mpl as cfv

right_sym = 5
top_sym = 15
fixed_wall = 10
left_wall_line = 20
dirichlet_wall = 25
copper = 30
nylon = 40

geom = cfg.Geometry() # copper

geom.point([0.0, 0.5*L-a-b])
geom.point([a, 0.5*L-a-b])
geom.point([a, 0.5*L-a-b-h])
geom.point([a+t, 0.5*L-a-b-h])
geom.point([a+t, 0.5*L-a-b])
geom.point([c+d, 0.5*L-a-b])
geom.point([c+d, 0.0])
geom.point([a+c+d, 0.0])
geom.point([L-2*d, 0.3*L-d])
geom.point([L, 0.3*L-d])
geom.point([L, 0.3*L])
geom.point([L-2*d, 0.3*L])
geom.point([a+c+d, d])
geom.point([a+c+d, 0.5*L-b-d])
geom.point([a+c, 0.5*L-b])
geom.point([a, 0.5*L-b])
geom.point([a, 0.5*L])
geom.point([0.0, 0.5*L])
geom.point([0.0, 0.5*L-b])
geom.point([0.0, 0.0]) # tis one has to be last


# Copper

geom.spline([0, 1])
geom.spline([1, 2])
geom.spline([2, 3])
geom.spline([3, 4])
geom.spline([4, 5])
geom.spline([5, 6])
geom.spline([6, 7])
geom.spline([7, 8], marker=dirichlet_wall)
geom.spline([8, 9], marker=dirichlet_wall)
geom.spline([9, 10], marker=right_sym)
geom.spline([10, 11], marker=dirichlet_wall)
geom.spline([11, 12], marker=dirichlet_wall)
geom.spline([12, 13], marker=dirichlet_wall)
geom.spline([13, 14], marker=dirichlet_wall)
geom.spline([14, 15], marker=dirichlet_wall)
geom.spline([15, 16], marker=dirichlet_wall)
geom.spline([16, 17], marker=top_sym)
geom.spline([17, 18], marker=left_wall_line)
geom.spline([18, 0], marker=fixed_wall)

# Nylon
geom.spline([19, 0], marker=fixed_wall)
geom.spline([6, 19])

# Create two Surfaces. Loop all point numbers except origin, hence minus one
lines = [i for i in range(len(geom.points)-1)]
geom.surface(lines, marker=copper)
geom.surface([0,1,2,3,4,5,20,19], marker=nylon)
# cfv.draw_geometry(geom)
# cfv.show_and_wait()