"""
TODO
"""

import numpy as np
import pyvista as pv

# TODO
pv.global_theme.notebook = False 


def _show_grid(environment, settings):
    """
    TODO
    """
    grid = _get_grid(environment.GRID)

    plotter = pv.Plotter()
    _ = plotter.add_mesh(grid, show_edges=True)
    plotter.show()

    return None


def _show_mask(environment, settings):
    """
    TODO
    """
    vertices = environment.MASK.polygon.vertices
    mask     = environment.MASK.mask

    full_z0 = np.full((1, vertices.shape[0]), environment.GRID.z0)
    points = np.concatenate((vertices, full_z0.T), axis=1)
    
    poly = pv.PolyData(points)
    lines0 = np.ones(len(points)-1, dtype=int) * 2
    lines1 = np.arange(len(points)-1, dtype=int)
    lines2 = np.arange(1,len(points), dtype=int)
    lines = np.column_stack((lines0, lines1, lines2))
    poly.lines = lines

    grid = _get_grid(environment.GRID)
    grid.cell_data["values"] = mask.flatten(order="F")
    grid = grid.cast_to_unstructured_grid()
    ghosts = np.argwhere(grid["values"] == 0)
    grid_ = grid.remove_cells(ghosts)

    plotter = pv.Plotter()
    _ = plotter.add_mesh(grid_, show_edges=True)
    _ = plotter.add_mesh(poly, color="tan", line_width=10, show_edges=True)
    plotter.show()

    return None


def _show_geology(environment, settings):
    """
    TODO
    """
    grid = _get_grid(environment.GRID)
    grid.cell_data["values"] = environment.GEOLOGY.data.flatten(order="F")

    if 'ghost' in settings:
        grid = _get_ghost_grid(grid, settings['ghost'])

    plotter = pv.Plotter()
    _ = plotter.add_mesh(grid, show_edges=True)
    plotter.show(cpos="xy")

    return None

def _show_topography(environment, settings):
#     """
#     TODO
#     """
#     grid = _get_grid(environment.GRID)
#     grid.cell_data["values"] = environment.GEOLOGY.data.flatten(order="F")
#     grid = grid.cast_to_unstructured_grid()
#     ghosts = np.argwhere(grid["values"] == 0)
#     grid_ = grid.remove_cells(ghosts)

#     plotter = pv.Plotter()
#     _ = plotter.add_mesh(grid_, show_edges=True)
#     plotter.show(cpos="xy")

    return None

def _show_numpy_array(environment, settings):
    """
    TODO
    """
    # array = settings['data'] # TODO
    # grid = _get_grid(environment.GRID)
    # grid.cell_data["values"] = array.flatten(order="F")
    
    # grid_ = grid.cast_to_unstructured_grid()
    # ghosts = np.argwhere(grid_["values"] == 0)
    # grid_ = grid_.remove_cells(ghosts)

    # plotter = pv.Plotter()
    # _ = plotter.add_mesh(grid, show_edges=True)
    # plotter.show(cpos="xy")

#     FRACTURES = domain.FRACTURES.location
# print('Nbr fractures:', len(FRACTURES))
# POLYGONS = []
# for fracture in FRACTURES:
#     x, y, z = fracture.get_position()
#     a, b, c = fracture.get_normal()
#     rad     = fracture.radius
#     POLYGONS.append(pv.Polygon(center=(x, y, z), radius=rad, normal=(a, b, c), n_sides=definition))

    return None

######################################################################

def _get_grid(grid):
    """
    TODO
    """
    mesh = pv.UniformGrid()
    mesh.dimensions = np.array((grid.nx, grid.ny, grid.nz)) + 1
    mesh.origin     = (grid.x0 - grid.dx/2, grid.y0 - grid.dy/2, grid.z0 - grid.dz/2)
    mesh.spacing    = (grid.dx, grid.dy, grid.dz)
    return mesh

def _get_ghost_grid(mesh, ghost_list):
    """
    TODO
    """
    mesh   = mesh.cast_to_unstructured_grid()
    ghosts = np.argwhere(np.isin(mesh["values"], ghost_list))
    mesh   = mesh.remove_cells(ghosts)
    return mesh
