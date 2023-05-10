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

# Get all pairs of dofs in lines along a marked spline.
def lines_along_spline(dofs, bdofs, edof, coords, spline_marker):
  lines_dofs = []
  line_coords = []
  for bdof in bdofs[spline_marker]:
    for EDOF in edof:
      if bdof in EDOF:
        handleEdof(bdofs, EDOF, spline_marker, lines_dofs)

  for line in lines_dofs:
    line_dof_coords = [0, 0] # Placeholder zeroes
    for coord, dof in zip(coords, dofs):
      if dof not in line: 
        continue
      line_dof_coords[line.index(dof)] = coord
      
    line_coords.append(line_dof_coords)

  return lines_dofs, line_coords


def handleEdof(bdofs, EDOF, spline_marker, lines_dofs):
  def line_exists(line):
    for ldofs in lines_dofs:
      if line[0] in ldofs and line[1] in ldofs:
        return True
    return False

  matches = []
  for el_dof in EDOF:
    if el_dof not in bdofs[spline_marker]:
      continue
    matches.append(el_dof)

  if len(matches) == 2 and not line_exists(matches):
    lines_dofs.append(matches)
    return True
  return False
