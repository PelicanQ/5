import numpy as np

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
