import numpy as np
import calfem.core as cfc
import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis_mpl as cfv
import calfem.utils as cfu
from geom import g, left_wall, left_wall_line, copper
from lines_along_spline import lines_along_spline
from constants import D1, D2, ny2, ny1, ep

#1 copper
#2 nylon

mesh = cfm.GmshMesh(g)

mesh.el_type = 2
mesh.dofs_per_node = 1 # three cuz temp and poss'es
mesh.el_size_factor = 0.15

coords, edof, dofs, bdofs, elementmarkers = mesh.create()

# Draw the mesh and geom.
cfv.draw_geometry(g)
cfv.drawMesh(
  coords=coords,
  edof=edof,
  dofs_per_node=mesh.dofsPerNode,
  el_type=mesh.elType,
  filled=True,
  title="Yobunga",
)

cfv.show_and_wait()
nDofs = np.size(dofs)

ex, ey = cfc.coordxtr(edof, coords, dofs)
K = np.zeros([nDofs,nDofs])

lines_dof, line_coords = lines_along_spline(dofs, bdofs, edof, coords, left_wall_line)
print(lines_dof)
print(line_coords)


## Non functional code below


for eltopo, elx, ely, marker in zip(edof, ex, ey, elementmarkers):
  D = D1 if marker == copper else D2
  Ke = cfc.flw2te(elx, ely, ep, D) # 3x3 för stat
  cfc.assem(eltopo, K, Ke)

bc = np.array([], 'i')
bcVal = np.array([], 'f')

bc, bcVal = cfu.applybc(bdofs, bc, bcVal, left_wall, 0.0, 2)
bc, bcVal = cfu.applybc(bdofs, bc, bcVal, left_wall, 0.0, 2)
print(bc)
print(bcVal)


# Först en frihetsgrad på stationärt
# Kolla Konvektionsproblemet i boken
# Fb integralen, som sedan adderas till Ka