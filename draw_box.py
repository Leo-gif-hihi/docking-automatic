from pymol.cgo import *
from pymol import cmd

def draw_box(x=12.385, y=48.136, z=26.209, sx=30.0, sy=30.0, sz=30.0, name="docking_box", r=1.0, g=1.0, b=0.0, linewidth=2.0):
    """
    DESCRIPTION
        Draws a 3D wireframe docking box in PyMOL.

    USAGE
        draw_box [x, y, z, [sx, sy, sz, [name, [r, g, b, [linewidth]]]]]

    EXAMPLES
        draw_box
        draw_box 10.5, 25.0, -5.2, 20, 20, 20
        draw_box 10.5, 25.0, -5.2, 20, 24, 18, name=my_box, r=0, g=1, b=1
    """
    # 1. Convert string arguments from PyMOL command line to floats
    cx, cy, cz = float(x), float(y), float(z)
    size_x, size_y, size_z = float(sx), float(sy), float(sz)
    r, g, b = float(r), float(g), float(b)
    lw = float(linewidth)

    # 2. Calculate min and max bounds for each axis
    min_x, max_x = cx - size_x / 2.0, cx + size_x / 2.0
    min_y, max_y = cy - size_y / 2.0, cy + size_y / 2.0
    min_z, max_z = cz - size_z / 2.0, cz + size_z / 2.0

    # 3. Construct 12 wireframe box edges using CGO LINES
    box = [
        LINEWIDTH, lw,
        BEGIN, LINES,
        COLOR, r, g, b,

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

    # 4. Delete existing object with the same name so it updates cleanly, then load
    if name in cmd.get_names("objects"):
        cmd.delete(name)
    cmd.load_cgo(box, name)
    
    print(f"Loaded '{name}': Center=({cx}, {cy}, {cz}), Size=({size_x}, {size_y}, {size_z})")

# Register the function as a PyMOL command
cmd.extend("draw_box", draw_box)