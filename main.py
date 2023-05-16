import calfem.core as cfc
import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.utils as cfu
import calfem.vis_mpl as cfv
import numpy as np

from constants import (D1, D2, T_inf, a_c, cp1, cp2, h_flow, ny2, rho1, rho2, thickness, T_0, E1, E2, ny1, ny2, element_size_factor)
from geom import copper, dirichlet_wall, geom, left_wall_line
from lines_along_spline import extract_lines, line_length
from plantml import plantml

#1 copper
#2 nylon

mesh = cfm.GmshMesh(geom)

mesh.el_type = 2
mesh.dofs_per_node = 1 # three cuz temp and poss'es
mesh.el_size_factor = element_size_factor
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
K_h = np.zeros([nDofs,nDofs])

for eltopo, elx, ely, marker in zip(edof, ex, ey, elementmarkers):
  D = D1 if marker == copper else D2
  Ke = cfc.flw2te(elx, ely, [thickness], D) # 3x3 för stat
  cfc.assem(eltopo, K_h, Ke)

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

f = f_c + f_h
K = K_h + K_c

# cfv.plt.spy(K_c)

# We now set the Dirchlet boundary conditions. Create two lists of dofs and corresponding vale
bc_dofs = np.array(bdofs[dirichlet_wall]) 
bc_vals = np.ones_like(bc_dofs) * T_inf

# Solve the matrix equation
temps_stat = np.linalg.solve(K, f)
T_stat_max = np.max(temps_stat)

# Draw stationary solution
# cfv.draw_nodal_values_shaded(a, coords, edof, 'Stationary solution', mesh.dofs_per_node, mesh.el_type, True)
# cfv.colorbar()
# cfv.show_and_wait()

# input('Press any key for task B')

#################### Lets go task B!##################

C = np.zeros([nDofs, nDofs])

for el_edof, el_ex, el_ey, el_marker in zip(edof, ex, ey, elementmarkers):
  rho = rho1 if el_marker == copper else rho2
  cp = cp1 if el_marker == copper else cp2
  C_e = plantml(el_ex, el_ey, thickness * rho * cp)
  
  cfc.assem(el_edof, C, C_e)

# cfv.plt.spy(C)
t0 = 0
t1 = 80
dt = 0.05

tt = np.arange(t0, t1, dt)
temps = np.empty((nDofs, tt.size))
temps[:, 0] = np.ones(nDofs) * T_0

for idx in range(tt.size - 1):
  prev_temps = temps[:, idx] 
  f_1D = np.reshape(f, (nDofs))
  next_temps = np.linalg.solve(C + dt*K, C@prev_temps + dt*f_1D)
  temps[:, idx+1] = next_temps

max_temps = np.max(temps, axis=0)
index_T90 = np.argmax(max_temps >  0.9 * T_stat_max)  

five_temps = [
  temps[:, 0],
  temps[:, int(index_T90 * 0.03 / 4) * 1],
  temps[:, int(index_T90 * 0.03 / 4) * 2],
  temps[:, int(index_T90 * 0.03 / 4) * 3],
  temps[:, int(index_T90 * 0.03)]
]
# for temps in five_temps: 
#   cfv.figure()
#   cfv.draw_nodal_values_shaded(temps, coords, edof, 'asdf', mesh.dofs_per_node, mesh.el_type)
#   cfv.colorbar()
#   cfv.show()

print('Time to 90% max ' + str(tt[index_T90]) + ' seconds')

# Först en frihetsgrad på stationärt
# Kolla Konvektionsproblemet i boken
# Fb integralen, som sedan adderas till Ka

#transient:
#C-matris från plantml
#välj deltat rimligt
#trapetsmetod med theta=1, nudiffa allt:
#a_n+1=(C+delta_t)^-1*(Ca_n+delta_t*f_k+1)