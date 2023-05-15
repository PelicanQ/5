import numpy as np
import matplotlib.pylab as plt
import calfem.core as cfc
import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis_mpl as cfv
import calfem.utils as cfu
from geom import geom, fixed_wall, left_wall_line, copper, dirichlet_wall
from lines_along_spline import extract_lines, line_length
from constants import D1, D2, ny2, ny1, h, t, thickness, T_inf, h_flow, a_c

#1 copper
#2 nylon

mesh = cfm.GmshMesh(geom)

mesh.el_type = 2
mesh.dofs_per_node = 1 # three cuz temp and poss'es
mesh.el_size_factor = 0.15
mesh.return_boundary_elements = True
coords, edof, dofs, bdofs, elementmarkers, belms = mesh.create()

# This makes it easier to accees eg. the coordinate of dof 45
dof_to_index = {}
for idx, dof in enumerate(dofs):
  dof_to_index[dof[0]] = idx

# Draw the mesh and geom.

# cfv.draw_geometry(g)

# cfv.drawMesh(
#   coords=coords,
#   edof=edof,
#   dofs_per_node=mesh.dofsPerNode,
#   el_type=mesh.elType,
#   filled=True,
#   title="Yobunga",
#   show_nodes=True
# )

# cfv.show_and_wait()

nDofs = np.size(dofs)
ex, ey = cfc.coordxtr(edof, coords, dofs)

# Create the global stiffness matrix
K = np.zeros([nDofs,nDofs])
for eltopo, elx, ely, marker in zip(edof, ex, ey, elementmarkers):
  D = D1 if marker == copper else D2
  Ke = cfc.flw2te(elx, ely, [thickness], D) # 3x3 för stat
  cfc.assem(eltopo, K, Ke)

# Now we consider load vector. We need coordinates of dofs along a marked spline
lines_coords_flux = []
lines_dofs_flux = []
lines_coords_conv = []
lines_dofs_conv = []

extract_lines(belms[left_wall_line], lines_coords_flux, lines_dofs_flux, coords, dof_to_index)
extract_lines(belms[dirichlet_wall], lines_coords_conv, lines_dofs_conv, coords, dof_to_index)

# Init the boundary convection and flux
f_c = np.zeros((nDofs, 1))
f_h = np.zeros((nDofs, 1))

# First f_h
for line_dofs, line_coords in zip(lines_dofs_flux, lines_coords_flux):
  f_val = -h_flow * thickness * line_length(line_coords)/2 # Was minus sign but resulted in colder than 18
  f_h[dof_to_index[line_dofs[0]]] += f_val
  f_h[dof_to_index[line_dofs[1]]] += f_val
  
# Then f_c and K_c
K_c = np.zeros([nDofs,nDofs])

for line_dofs, line_coords in zip(lines_dofs_conv, lines_coords_conv):
  common_factor = a_c * thickness * line_length(line_coords)
  dof0_idx = dof_to_index[line_dofs[0]]
  dof1_idx = dof_to_index[line_dofs[1]] 

  K_c[dof0_idx][dof0_idx] += common_factor/3
  K_c[dof0_idx][dof1_idx] += common_factor/6
  K_c[dof1_idx][dof0_idx] += common_factor/6
  K_c[dof1_idx][dof1_idx] += common_factor/3
  
  f_val = a_c * thickness * T_inf * line_length(line_coords)/2 # Was minus sign but resulted in colder than 18
  f_c[dof0_idx] += f_val
  f_c[dof1_idx] += f_val

# cfv.plt.spy(K_c)

# We now set the Dirchlet boundary conditions. Create two lists of dofs and corresponding vale
bc_dofs = np.array(bdofs[dirichlet_wall]) 
bc_vals = np.ones_like(bc_dofs) * T_inf

# Solve the matrix equation
a = np.linalg.solve(K + K_c, f_h + f_c)

cfv.draw_nodal_values_shaded(a, coords, edof, 'Yani gås', mesh.dofs_per_node, mesh.el_type, True)
cfv.colorbar()
cfv.show_and_wait()


# Först en frihetsgrad på stationärt
# Kolla Konvektionsproblemet i boken
# Fb integralen, som sedan adderas till Ka

#transient:
#C-matris från plantml
#välj deltat rimligt
#trapetsmetod med theta=1, nudiffa allt:
#a_n+1=(C+delta_t)^-1*(Ca_n+delta_t*f_k+1)
#