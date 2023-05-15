import calfem.core as cfc
import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.utils as cfu
import calfem.vis_mpl as cfv
import numpy as np
from constants import (D1, D2, T_inf, a_c, cp1, cp2, h_flow, ny2, rho1, rho2, thickness, T_0, E1, E2, ny1, ny2, ptype, ep)
from geom import copper, dirichlet_wall, geom, left_wall_line, fixed_wall
from lines_along_spline import extract_lines, line_length

mesh = cfm.GmshMesh(geom)

mesh.el_type = 2
mesh.dofs_per_node = 2 # three cuz temp and poss'es
mesh.el_size_factor = 0.15
mesh.return_boundary_elements = True
coords, edof, dofs, bdofs, elementmarkers, belms = mesh.create()

# cfv.draw_geometry(geom, True, True, True)
# cfv.show()

# This makes it easier to accees eg. the coordinate of dof 45
dof_to_index = {}
for idx, dof in enumerate(dofs):
  dof_to_index[dof[0]] = idx

nDofs = np.size(dofs)
ex, ey = cfc.coordxtr(edof, coords, dofs)




K = np.empty((nDofs, nDofs))

for el_x, el_y, el_edof, marker in zip(ex, ey, edof, elementmarkers):
  E = E1 if marker == copper else E2
  ny = ny1 if marker == copper else ny2
  D = cfc.hooke(ptype, E, ny)[np.ix_([0, 1, 3], [0, 1, 3])]
  K_e = cfc.plante(el_x, el_y, ep, D)
  cfc.assem(el_edof, K, K_e)

bc = np.array([],'i')
bcVal = np.array([],'i')

bc, bcVal = cfu.applybc(bdofs, bc, bcVal, fixed_wall, 0.0)
bc, bcVal = cfu.applybc(bdofs, bc, bcVal, left_wall_line, 0.0)

F = np.zeros([nDofs,1])

# SAKNAR KOPPLING TILL TEMPERATUR
# cfu.applyforce(bdofs, f, mark_load, value = -10e5, dimension=2)

a, f = cfc.solveq(K, F, bc, bcVal)

a = np.random.random(a.shape)*0.1

ed = cfc.extractEldisp(edof, a)
von_mises = []

for el_x, el_y, eltopo, marker, disp in zip(ex, ey, edof, elementmarkers, ed):
  E = E1 if marker == copper else E2
  ny = ny1 if marker == copper else ny2
  D = cfc.hooke(ptype, E, ny)[np.ix_([0, 1, 3], [0, 1, 3])]

  es, et = cfc.plants(el_x, el_y, ep, D, disp)
  von_mises.append( np.math.sqrt( pow(es[0,0],2) - es[0,0]*es[0,1] + pow(es[0,1], 2) + 3*pow(es[0,2],2) ) )


cfv.draw_element_values(von_mises, coords, edof, mesh.dofs_per_node, mesh.el_type)
cfv.colorbar()
cfv.figure()
cfv.draw_displacements(a, coords, edof, mesh.dofs_per_node, mesh.el_type, True)
cfv.show()
