import calfem.core as cfc
import calfem.mesh as cfm
import calfem.vis_mpl as cfv
import numpy as np

from constants import (D1, D2, T_inf, a_c, cp1, cp2, h_flow, rho1, rho2, thickness, T_0, element_size_factor)
from geom import copper, dirichlet_wall, geom, left_wall_line
from util import extract_lines, line_length, mirror_coords
from plantml import plantml

# 1 = copper
# 2 = nylon

mesh = cfm.GmshMesh(geom)

mesh.el_type = 2
mesh.dofs_per_node = 1 # Temperature is the only DOF
mesh.el_size_factor = element_size_factor
mesh.return_boundary_elements = True
coords, edof, dofs, bdofs, elementmarkers, belms = mesh.create()

# This makes it easier to accees eg. the coordinate of dof 45. No assumption of dof number order
dof_to_index = {}
for idx, dof in enumerate(dofs):
  dof_to_index[dof[0]] = idx

# Draw the mesh to verify

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

# Now let's find the stationary temperature distrubution

nDofs = np.size(dofs)
ex, ey = cfc.coordxtr(edof, coords, dofs)

# Create the global stiffness matrix (no convection)
K_h = np.zeros([nDofs,nDofs])

# Assemble global stiffness from element stiffness matricies
for eltopo, elx, ely, marker in zip(edof, ex, ey, elementmarkers):
  D = D1 if marker == copper else D2
  Ke = cfc.flw2te(elx, ely, [thickness], D) # 3x3
  cfc.assem(eltopo, K_h, Ke)

# Now we consider load vector. We need coordinates of dofs along a marked spline. No assumption of equal lengths
lines_coords_flux = []
lines_dofs_flux = []
lines_coords_conv = []
lines_dofs_conv = []

# We want the pairs of coordinates and dofs that make up splines along a marker
extract_lines(belms[left_wall_line], lines_coords_flux, lines_dofs_flux, coords, dof_to_index)
extract_lines(belms[dirichlet_wall], lines_coords_conv, lines_dofs_conv, coords, dof_to_index)

# Init the boundary convection and flux
f_c = np.zeros((nDofs, 1))
f_h = np.zeros((nDofs, 1))

# First heat flux, f_h
for line_dofs, line_coords in zip(lines_dofs_flux, lines_coords_flux):
  f_val = -h_flow * thickness * line_length(line_coords)/2 
  f_h[dof_to_index[line_dofs[0]]] += f_val
  f_h[dof_to_index[line_dofs[1]]] += f_val
  
# Then convection load vector and stiffness matrix f_c and K_c
K_c = np.zeros([nDofs,nDofs])

for line_dofs, line_coords in zip(lines_dofs_conv, lines_coords_conv):
  common_factor = a_c * thickness * line_length(line_coords)
  dof0_idx = dof_to_index[line_dofs[0]]
  dof1_idx = dof_to_index[line_dofs[1]] 

  # Manual assembly
  K_c[dof0_idx][dof0_idx] += common_factor/3
  K_c[dof0_idx][dof1_idx] += common_factor/6
  K_c[dof1_idx][dof0_idx] += common_factor/6
  K_c[dof1_idx][dof1_idx] += common_factor/3
  
  f_val = a_c * thickness * T_inf * line_length(line_coords)/2 
  f_c[dof0_idx] += f_val
  f_c[dof1_idx] += f_val

# Total
f = f_c + f_h
K = K_h + K_c

# cfv.plt.spy(K_c)

# We now set the Dirchlet boundary conditions. Create a list of dofs and on list of corresponding dirichlet values
bc_dofs = np.array(bdofs[dirichlet_wall]) 
bc_vals = np.ones_like(bc_dofs) * T_inf

# Solve the matrix equation
temps_stat = np.linalg.solve(K, f)

# Draw stationary solution. Mirror the plot into all quadrants 
quadrants = mirror_coords(coords)

fig_size=(12, 5)
cfv.matplotlib.rc('font', size=14)
# cfv.figure(fig_size=fig_size)

# for quadrant_coords in quadrants:
#   cfv.draw_nodal_values_shaded(temps_stat, quadrant_coords, edof, 'Stationary temperature distrubution [°C]', mesh.dofs_per_node, mesh.el_type, True)


# cfv.plt.xlabel('X Position [m]')
# cfv.plt.ylabel('Y Position [m]')
# cbar = cfv.colorbar()

# cfv.show_and_wait()

# input('Press Enter key for task B')



#################### Now task B ####################

C = np.zeros([nDofs, nDofs])

for el_edof, el_ex, el_ey, el_marker in zip(edof, ex, ey, elementmarkers):
  rho = rho1 if el_marker == copper else rho2
  cp = cp1 if el_marker == copper else cp2
  C_e = plantml(el_ex, el_ey, thickness * rho * cp)
  
  cfc.assem(el_edof, C, C_e)

# cfv.plt.spy(C)


t0 = 0
t1 = 90
dt = 0.05 # seconds

tt = np.arange(t0, t1, dt) # Time grid
temps = np.empty((nDofs, tt.size)) # Along one row will be temperature in one dof at each time in time grid
temps[:, 0] = np.ones(nDofs) * T_0 # At first time, all dofs are surrounding temp

f_1D = np.reshape(f, (nDofs))

# Trapezoidal time step. This is NumDiff
for idx in range(tt.size - 1):
  prev_temps = temps[:, idx] 
  next_temps = np.linalg.solve(C + dt*K, C@prev_temps + dt*f_1D)
  temps[:, idx+1] = next_temps


max_temps = np.max(temps, axis=0) # Find the hottest point for every time.
T_stat_max = np.max(temps_stat) # Find the hottest point in stationary state
index_T90 = np.argmax(max_temps >  T_0 + 0.9 * (T_stat_max - T_0))  # When do we reach 90% of stationary max?
print(f'Max stationary temperature: {T_stat_max}C')

# Grab five (approximately) equally spaced snapshots.
time_indicies = [ int(index_T90 * 0.03 * i / 4) for i in range(5) ]
five_temps = [ temps[:, time_index] for time_index in time_indicies ]

# for idx, (trans_temps, time_index) in enumerate(zip(five_temps, time_indicies)): 
#   if idx % 2 == 0:
#     cfv.figure(fig_size=(12, 4))
  
#   cfv.subplot(1, 2, (idx % 2) + 1)
  
#   for quadrant_coords in quadrants: 
#     cfv.draw_nodal_values_shaded(trans_temps, quadrant_coords, edof, f'Temperature distrubution at {round(tt[time_index], 3)} seconds [°C]', mesh.dofs_per_node, mesh.el_type)

#   cfv.plt.xlabel('X Position [m]')
#   cfv.plt.ylabel('Y Position [m]')
#   cfv.colorbar(location='top')

# cfv.show()


print('Time to 90% max ' + str(tt[index_T90]) + ' seconds')