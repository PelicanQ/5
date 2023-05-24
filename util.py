import numpy as np

def extract_lines(belms, lines_coords, lines_dofs, coords, dof_to_index):
  for belm in belms:
    dof0 = belm['node-number-list'][0]
    dof1 = belm['node-number-list'][1]
    
    dof0_coord = coords[dof_to_index[dof0]]
    dof1_coord = coords[dof_to_index[dof1]]
    
    # Very readable!
    lines_dofs.append([dof0, dof1])
    lines_coords.append([dof0_coord, dof1_coord])

def line_length(line_coords):
  y0 = line_coords[0][1]
  y1 = line_coords[1][1]
  x0 = line_coords[0][0]
  x1 = line_coords[1][0]
  return np.sqrt((y1-y0)**2 + (x1-x0)**2)

def mirror_coords(coords):
  q1_coords = np.empty_like(coords)
  q2_coords = np.empty_like(coords)
  q3_coords = coords
  q4_coords = np.empty_like(coords)

  sym_x = 0.005
  sym_y = 0.0025
  for idx, coords in enumerate(coords):
    q1_coords[idx, :] = [2*sym_x - coords[0], 2*sym_y - coords[1]]          # 1st quadrant
    q2_coords[idx, :] = [coords[0],         2*sym_y - coords[1]]          # 2nd quadrant
    q4_coords[idx, :] = [2*sym_x - coords[0], coords[1]]                  # 4th quadrant

  return q1_coords, q2_coords, q3_coords, q4_coords

def mirror_diplacements(displacements):
  q1_disp = np.copy(displacements)
  q2_disp = np.copy(displacements)
  q3_disp = np.copy(displacements)
  q4_disp = np.copy(displacements)

  q1_disp = -q1_disp             # 1st quadrant
  q2_disp[1::2] = -q2_disp[1::2] # 2nd quadrant
  q4_disp[0::2] = -q4_disp[0::2] # 4th quadrant

  return q1_disp, q2_disp, q3_disp, q4_disp
