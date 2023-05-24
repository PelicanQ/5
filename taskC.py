import calfem.core as cfc
import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.utils as cfu
import calfem.vis_mpl as cfv
import numpy as np
from constants import (D1, D2, T_inf, a_c, cp1, cp2, h_flow, ny2, rho1, rho2, element_size_factor, T_0, E1, E2, ny1, ny2, ptype, ep, alfa1, alfa2)
from geom import copper, dirichlet_wall, geom, left_wall_line, fixed_wall, right_sym, top_sym
from lines_along_spline import extract_lines, line_length
import main 

# We gotta make a new mesh
mesh = cfm.GmshMesh(geom)
mesh.el_type = 2
mesh.dofs_per_node = 2 # three cuz temp and poss'es
mesh.el_size_factor = element_size_factor 
mesh.return_boundary_elements = True
coords, edof, dofs, bdofs, elementmarkers, belms = mesh.create()

nDofs = np.size(dofs)
ex, ey = cfc.coordxtr(edof, coords, dofs)

#Now make the K
K = np.empty((nDofs, nDofs))

for el_x, el_y, el_edof, marker in zip(ex, ey, edof, elementmarkers):
  E = E1 if marker == copper else E2
  ny = ny1 if marker == copper else ny2
  D = cfc.hooke(ptype, E, ny)[np.ix_([0, 1, 3], [0, 1, 3])]
  K_e = cfc.plante(el_x, el_y, ep, D)
  cfc.assem(el_edof, K, K_e)

bc = np.array([],'i')
bcVal = np.array([],'i')

# Fix one wall in both dimensions and constrain the symmerty lines in the perpendicular dimension 
bc, bcVal = cfu.applybc(bdofs, bc, bcVal, fixed_wall, 0.0)
bc, bcVal = cfu.applybc(bdofs, bc, bcVal, left_wall_line, 0.0)
bc, bcVal = cfu.applybc(bdofs, bc, bcVal, right_sym, 0, 1)
bc, bcVal = cfu.applybc(bdofs, bc, bcVal, top_sym, 0, 2)

# Lets organize temperatures according to element dofs
e_T = np.empty(main.edof.shape)

for idx, (el_x, el_y, el_topo, marker) in enumerate(zip(main.ex, main.ey, main.edof, main.elementmarkers)):
  e_T[idx, 0] = main.temps_stat[main.dof_to_index[el_topo[0]]]
  e_T[idx, 1] = main.temps_stat[main.dof_to_index[el_topo[1]]]
  e_T[idx, 2] = main.temps_stat[main.dof_to_index[el_topo[2]]]

delta_T = np.mean(e_T, 1) - T_inf

F = np.zeros([nDofs,1])

for el_x, el_y, eltopo, marker, dT in zip(ex, ey, edof, elementmarkers, delta_T):
  E = E1 if marker == copper else E2
  ny = ny1 if marker == copper else ny2
  alfa = alfa1 if marker == copper else alfa2

  D = cfc.hooke(ptype, E, ny)[np.ix_([0, 1, 3], [0, 1, 3])]
  ep_0 = (1+ny)*alfa*dT*np.array([[1], [1], [0]])
  es = D * ep_0
  f_e = cfc.plantf(el_x, el_y, ep, np.transpose(es))

  f_indexes = eltopo - 1
  F[f_indexes] +=  np.reshape(f_e, (6,1))

a, f = cfc.solveq(K, F, bc, bcVal)

# We now have displacements. Lets compute von Mises strains

element_displacements = cfc.extract_ed(edof, a)
von_mises_el = np.empty(edof.shape[0])
von_mises_node = np.empty(coords.shape[0])

for idx, (el_x, el_y, eltopo, marker, ed, dT) in enumerate(zip(ex, ey, edof, elementmarkers, element_displacements, delta_T)):
  E = E1 if marker == copper else E2
  ny = ny1 if marker == copper else ny2
  alfa = alfa1 if marker == copper else alfa2

  D = cfc.hooke(ptype, E, ny)[np.ix_([0, 1, 3], [0, 1, 3])]

  ep_0 = (1+ny)*alfa*dT*np.array([[1], [1], [0]])
  es, et = cfc.plants(el_x, el_y, ep, D, ed)
  es = np.transpose(es)
  sigma = es - D*ep_0
  sigma = np.transpose(sigma)

  sigma_zz = ny*(sigma[0, 0] + sigma[0, 1]) - alfa*E*dT

  von_mises_el[idx] = np.math.sqrt(
    pow(sigma[0, 0], 2) + pow(sigma[0, 1], 2) + pow(sigma_zz, 2) 
    - sigma[0, 0]*sigma[0, 1] - sigma[0, 0]*sigma_zz - sigma[0, 1]*sigma_zz 
    + 3*pow(sigma[0, 2], 2) 
  )

for i, dof_pair in enumerate(dofs):
  row_matches, col_matches = np.nonzero(edof == dof_pair[0])
  von_mises_node[i] = np.sum(von_mises_el[row_matches])/len(row_matches)

# Unfortunately draw_nodal_values[_shaded] need a different edof based on node number
el_nodes = np.empty((edof.shape[0], 3))

for i, el_topo in enumerate(edof):
  node_inds = []
  for el_dof in el_topo:
    row_matches, col_matches = np.nonzero(dofs == el_dof)
    node_inds.append(row_matches[0])

  node_inds = [*set(node_inds)] #Remove duplicates
  el_nodes[i, :] = node_inds

# +1 because node number is node index +1  
cfv.draw_nodal_values_shaded(von_mises_node, coords, el_nodes+1, 'von Mises stress field [Pa]', mesh.dofs_per_node, mesh.el_type)
cfv.plt.xlabel('X Position [m]')
cfv.plt.ylabel('Y Position [m]')
cfv.colorbar()
cfv.figure()
cfv.draw_displacements(a, coords, edof, mesh.dofs_per_node, mesh.el_type, True, 1)
cfv.show()