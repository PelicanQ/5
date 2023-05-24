import calfem.core as cfc
import calfem.mesh as cfm
import calfem.utils as cfu
import calfem.vis_mpl as cfv
import numpy as np
from constants import (T_inf, ny2, element_size_factor, E1, E2, ny1, ny2, ptype, ep, alfa1, alfa2)
from geom import copper, geom, left_wall_line, fixed_wall, right_sym, top_sym
from util import mirror_coords, mirror_diplacements
import main 

# Now for getting displacements and stresses from stationary temperature distrubution 

# We must make a new mesh with 2 dofs per node
mesh = cfm.GmshMesh(geom)
mesh.el_type = 2
mesh.dofs_per_node = 2
mesh.el_size_factor = 0.01
mesh.return_boundary_elements = True
coords, edof, dofs, bdofs, elementmarkers, belms = mesh.create()

nDofs = np.size(dofs)
ex, ey = cfc.coordxtr(edof, coords, dofs)

# Initialize K matrix
K = np.empty((nDofs, nDofs))

# Assemble K from element stiffness matricies
for el_x, el_y, el_edof, marker in zip(ex, ey, edof, elementmarkers):
  E = E1 if marker == copper else E2
  ny = ny1 if marker == copper else ny2
  D = cfc.hooke(ptype, E, ny)[np.ix_([0, 1, 3], [0, 1, 3])]
  K_e = cfc.plante(el_x, el_y, ep, D)
  cfc.assem(el_edof, K, K_e)

# Initialize BC's
bc = np.array([], 'i')
bcVal = np.array([], 'i')

# Fix one wall in both dimensions 
bc, bcVal = cfu.applybc(bdofs, bc, bcVal, fixed_wall, 0.0)
bc, bcVal = cfu.applybc(bdofs, bc, bcVal, left_wall_line, 0.0)

# Constrain the symmerty lines in the perpendicular dimension
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

displacements, f = cfc.solveq(K, F, bc, bcVal)

# We now have displacements. Lets compute von Mises stresses

element_displacements = cfc.extract_ed(edof, displacements)
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

print(von_mises_el)

#Quote from lab manual: 'The stress of the nodal points can be approximated by taking the mean value of the stresses in the elements connected to a node'
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


print(f'Peak von Mises stress {np.max(von_mises_node)}Pa')

quadrants_coords = mirror_coords(coords)
quadrants_displacements = mirror_diplacements(displacements)

cfv.figure(fig_size=(12,5))
# +1 because node number is node index +1
for quadrant_coords in quadrants_coords:
  cfv.draw_nodal_values_shaded(von_mises_node, quadrant_coords, el_nodes+1, 'von Mises stress field [Pa]', mesh.dofs_per_node, mesh.el_type)

cfv.plt.xlabel('X Position [m]')
cfv.plt.ylabel('Y Position [m]')
cfv.colorbar()
cfv.show()

# Visualize displacements
cfv.figure(fig_size=(12,5))

for q_coords, q_disp in zip(quadrants_coords, quadrants_displacements):
  cfv.draw_displacements(q_disp, q_coords, edof, mesh.dofs_per_node, mesh.el_type, True, 1, title='Displaced mesh (black) overlaid undisplaced mesh (gray)')

cfv.plt.xlabel('X Position [m]')
cfv.plt.ylabel('Y Position [m]')
cfv.show()