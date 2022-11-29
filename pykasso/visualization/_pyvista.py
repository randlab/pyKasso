"""
TODO
"""

import numpy as np
import pyvista as pv

# TODO
pv.global_theme.notebook = False 

def _show_data(environment, feature, settings):
    """
    TODO
    """
    # Initial settings
    kwargs = {
        'scalars'    : None,
        'show_edges' : True,
    }
    #     'cmap' : "viridis",

    # Get the grid
    grid = _get_grid(environment.GRID)

    # Get the data
    if feature.lower() != 'grid':
        grid = _set_data(grid, getattr(environment, feature), 'data', 'data')
        kwargs['scalars'] = 'data'

    # Get the mask
    if getattr(environment, 'MASK') is not None:
        grid = _set_data(grid, getattr(environment, 'MASK'), 'mask', 'mask')

    # Ghost the data
    if 'ghost' in settings:
        if 'mask' in settings:
            ghosts = np.argwhere(np.logical_or(np.isin(grid["data"], settings['ghost']), (grid["mask"] == 0)))
        else:
            print(2)
            ghosts = np.argwhere(np.isin(grid["data"], settings['ghost']))
    else:
        if 'mask' in settings:
            ghosts = np.argwhere(grid["mask"] == 0)
        else:
            ghosts = None
    
    if ghosts is not None:
        grid = grid.remove_cells(ghosts)
    
    # Plot the data
    plotter = pv.Plotter()
    _ = plotter.add_mesh(grid, **kwargs)
    plotter.show(cpos='xy')

    return None


#     grid = _get_grid(environment.GRID)
#     array = settings['data'] # TODO
#     grid.cell_data["values"] = array.flatten(order="F")
    
#     grid_ = grid.cast_to_unstructured_grid()
#     ghosts = np.argwhere(grid_["values"] == 0)
#     grid_ = grid_.remove_cells(ghosts)

#     plotter = pv.Plotter()
#     _ = plotter.add_mesh(grid, show_edges=True)
#     plotter.show(cpos="xy")

#     return None

#     FRACTURES = domain.FRACTURES.location
# print('Nbr fractures:', len(FRACTURES))
# POLYGONS = []
# for fracture in FRACTURES:
#     x, y, z = fracture.get_position()
#     a, b, c = fracture.get_normal()
#     rad     = fracture.radius
#     POLYGONS.append(pv.Polygon(center=(x, y, z), radius=rad, normal=(a, b, c), n_sides=definition))

    # return None

######################################################################

def _get_grid(grid):
    """
    TODO
    """
    mesh = pv.UniformGrid()
    mesh.dimensions = np.array((grid.nx, grid.ny, grid.nz)) + 1
    mesh.origin     = (grid.x0 - grid.dx/2, grid.y0 - grid.dy/2, grid.z0 - grid.dz/2)
    mesh.spacing    = (grid.dx, grid.dy, grid.dz)
    mesh            = mesh.flip_z(inplace=False)
    mesh            = mesh.cast_to_unstructured_grid()
    return mesh

def _set_data(mesh, geologic_feature, attribute, label):
    """
    TODO
    """
    data = getattr(geologic_feature, attribute)
    mesh.cell_data[label] = data.flatten(order="F")
    return mesh

# def _show_mask(environment, settings):
#     """
#     TODO
#     """
#     vertices = environment.MASK.polygon.vertices
#     mask     = environment.MASK.mask

#     full_z0 = np.full((1, vertices.shape[0]), environment.GRID.z0)
#     points = np.concatenate((vertices, full_z0.T), axis=1)
    
#     poly = pv.PolyData(points)
#     lines0 = np.ones(len(points)-1, dtype=int) * 2
#     lines1 = np.arange(len(points)-1, dtype=int)
#     lines2 = np.arange(1,len(points), dtype=int)
#     lines = np.column_stack((lines0, lines1, lines2))
#     poly.lines = lines

#     grid = _get_grid(environment.GRID)
#     grid.cell_data["values"] = mask.flatten(order="F")
#     grid = grid.cast_to_unstructured_grid()
#     ghosts = np.argwhere(grid["values"] == 0)
#     grid_ = grid.remove_cells(ghosts)

#     plotter = pv.Plotter()
#     _ = plotter.add_mesh(grid_, show_edges=True)
#     _ = plotter.add_mesh(poly, color="tan", line_width=10, show_edges=True)
#     plotter.show()

#     return None