from pymol.cgo import *
from pymol import cmd
# 1. Define your parameters
center = [12.385, 48.136, 26.209]
size = [30.0, 30.0, 30.0]

# 2. Calculate min and max bounds for each axis
min_x, max_x = center[0] - size[0] / 2.0, center[0] + size[0] / 2.0
min_y, max_y = center[1] - size[1] / 2.0, center[1] + size[1] / 2.0
min_z, max_z = center[2] - size[2] / 2.0, center[2] + size[2] / 2.0

# 3. Construct 12 wireframe box edges using CGO LINES
box = [
    LINEWIDTH, 2.0,
    BEGIN, LINES,
    COLOR, 1.0, 1.0, 0.0,  # Yellow wireframe (R, G, B: 0.0 to 1.0)

    # Bottom square (z = min_z)
    VERTEX, min_x, min_y, min_z, VERTEX, max_x, min_y, min_z,
    VERTEX, max_x, min_y, min_z, VERTEX, max_x, max_y, min_z,
    VERTEX, max_x, max_y, min_z, VERTEX, min_x, max_y, min_z,
    VERTEX, min_x, max_y, min_z, VERTEX, min_x, min_y, min_z,

    # Top square (z = max_z)
    VERTEX, min_x, min_y, max_z, VERTEX, max_x, min_y, max_z,
    VERTEX, max_x, min_y, max_z, VERTEX, max_x, max_y, max_z,
    VERTEX, max_x, max_y, max_z, VERTEX, min_x, max_y, max_z,
    VERTEX, min_x, max_y, max_z, VERTEX, min_x, min_y, max_z,

    # Vertical connectors
    VERTEX, min_x, min_y, min_z, VERTEX, min_x, min_y, max_z,
    VERTEX, max_x, min_y, min_z, VERTEX, max_x, min_y, max_z,
    VERTEX, max_x, max_y, min_z, VERTEX, max_x, max_y, max_z,
    VERTEX, min_x, max_y, min_z, VERTEX, min_x, max_y, max_z,

    END
]

# 4. Load the CGO object into PyMOL
cmd.load_cgo(box, "docking_box")