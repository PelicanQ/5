import numpy as np
import calfem.core as cfc
import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.utils as cfu
from constants import a,b,h,t,d,c,L
import matplotlib as mpl
import calfem.vis_mpl as cfv

left_wall = 10
left_wall_line = 20
dirichlet_wall = 25
copper = 30
nylon = 40

g = cfg.Geometry() # copper

g.point([0.0, 0.5*L-a-b])
g.point([a, 0.5*L-a-b])
g.point([a, 0.5*L-a-b-h])
g.point([a+t, 0.5*L-a-b-h])
g.point([a+t, 0.5*L-a-b])
g.point([c+d, 0.5*L-a-b])
g.point([c+d, 0.0])
g.point([a+c+d, 0.0])
g.point([L-2*d, 0.3*L-d])
g.point([L, 0.3*L-d])
g.point([L, 0.3*L])
g.point([L-2*d, 0.3*L])
g.point([a+c+d, d])
g.point([a+c+d, 0.5*L-b-d])
g.point([a+c, 0.5*L-b])
g.point([a, 0.5*L-b])
g.point([a, 0.5*L])
g.point([0.0, 0.5*L])
g.point([0.0, 0.5*L-b])
g.point([0.0, 0.0]) # tis one has to be last


# Copper

g.spline([0, 1])
g.spline([1, 2])
g.spline([2, 3])
g.spline([3, 4])
g.spline([4, 5])
g.spline([5, 6])
g.spline([6, 7])
g.spline([7, 8], marker=dirichlet_wall)
g.spline([8, 9], marker=dirichlet_wall)
g.spline([9, 10])
g.spline([10, 11], marker=dirichlet_wall)
g.spline([11, 12], marker=dirichlet_wall)
g.spline([12, 13], marker=dirichlet_wall)
g.spline([13, 14], marker=dirichlet_wall)
g.spline([14, 15], marker=dirichlet_wall)
g.spline([15, 16], marker=dirichlet_wall)
g.spline([16, 17])
g.spline([17, 18], marker=left_wall_line)
g.spline([18, 0])

# Nylon
g.spline([19, 0])
g.spline([6, 19])

# Create two Surfaces. Loop all point numbers except origin, hence minus one
lines = [i for i in range(len(g.points)-1)]
g.surface(lines, marker=copper)
g.surface([0,1,2,3,4,5,20,19], marker=nylon)
# cfv.draw_geometry(g)
# cfv.show_and_wait()