"""
TODO
"""

import numpy as np
import pyvista as pv

# TODO
# - incorporer cette option
pv.global_theme.notebook = False
# - 'show_edges' : True
# - 'cmap' : "viridis",

def _show_data(environment, feature, settings):
    """
    TODO
    """
    ########################
    ### Initial settings ###
    ########################

    kwargs = {
        'scalars'    : None,
        'show_edges' : False,
    }

    attributes = [
        'slice',
        'show_grid',
        'inlets',
        'outlets',
        'tracers',
    ]

    for attribute in attributes:
        if attribute not in settings:
            settings[attribute] = False

 
    #########################   
    ### Geologic features ###
    #########################

    # Get the grid
    grid_i = _get_grid(environment.GRID)
    grid = grid_i

    # Get the data
    if feature.lower() != 'grid':
        grid = _get_data(grid, getattr(environment, feature), 'data', 'data')
        kwargs['scalars'] = 'data'

    # Get the mask
    if getattr(environment, 'MASK') is not None:
        grid = _get_data(grid, getattr(environment, 'MASK'), 'mask', 'mask')

    # Ghost the data
    if 'ghost' in settings:
        if 'mask' in settings:
            ghosts = np.argwhere(np.logical_or(np.isin(grid["data"], settings['ghost']), (grid["mask"] == 0)))
        else:
            ghosts = np.argwhere(np.isin(grid["data"], settings['ghost']))
    else:
        if 'mask' in settings:
            ghosts = np.argwhere(grid["mask"] == 0)
        else:
            ghosts = None
    
    if ghosts is not None:
        grid = grid.remove_cells(ghosts)

    ##############
    ### Points ###
    ##############

    if settings['inlets']:
        inlets = _get_points(environment.POINTS, 'inlets')

    if settings['outlets']:
        outlets = _get_points(environment.POINTS, 'outlets')

    # TODO
    if settings['tracers']:
        pass

    #####################
    ### Plot the data ###
    #####################

    plotter = pv.Plotter()

    # Data
    if settings['slice']:
        slices = grid.slice_orthogonal()
        _ = plotter.add_mesh(slices, **kwargs)
    else:
         _ = plotter.add_mesh(grid, **kwargs)

    # Grid
    _ = plotter.add_mesh(grid_i.outline(), color="k")

    if settings['show_grid']:
        plotter.show_grid()

    # Points
    if settings['inlets']:
        plotter.add_points(inlets, render_points_as_spheres=True, point_size=20, color='r')
    if settings['outlets']:
        plotter.add_points(outlets, render_points_as_spheres=True, point_size=20, color='b')

    plotter.show(cpos='xy')

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
    mesh            = mesh.flip_z(inplace=False)
    mesh            = mesh.cast_to_unstructured_grid()
    return mesh

def _get_data(mesh, geologic_feature, attribute, label):
    """
    TODO
    """
    data = getattr(geologic_feature, attribute)
    mesh.cell_data[label] = data.flatten(order="F")
    return mesh

def _get_points(points, kind):
    """
    TODO
    """
    points = points[points['label'] == kind][['x', 'y', 'z']].values
    points = points.astype('float32')
    cloud  = pv.wrap(points)
    return cloud


# FRACTURES = domain.FRACTURES.location
# print('Nbr fractures:', len(FRACTURES))
# POLYGONS = []
# for fracture in FRACTURES:
#     x, y, z = fracture.get_position()
#     a, b, c = fracture.get_normal()
#     rad     = fracture.radius
#     POLYGONS.append(pv.Polygon(center=(x, y, z), radius=rad, normal=(a, b, c), n_sides=definition))

    # return None

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


##################################################3
### DEBUG ###
#############

# def debug_plot_initialize(environment):
#     """
#     TODO
#     """

#     ### Call the plotter
#     p = pv.Plotter(shape=(2, 4), border=True, notebook=False)

#     ### Construct the grid
#     vtk = pv.UniformGrid()
#     vtk.dimensions = np.array((environment.GRID.nx, environment.GRID.ny, environment.GRID.nz)) + 1
#     vtk.origin     = (environment.GRID.x0 - environment.GRID.dx/2, environment.GRID.y0 - environment.GRID.dy/2, environment.GRID.z0 - environment.GRID.dz/2)
#     vtk.spacing    = (environment.GRID.dx, environment.GRID.dy, environment.GRID.dz)

#     features = ['geology', 'faults', 'fractures', 'field']

#     for (i, feature) in enumerate(features):
#         vtk = vtk.copy()
#         p.subplot(0, i)
#         p.add_text(feature, font_size=24)
#         try:
#             vtk['values'] = self.geology.geologic_features[feature].data.flatten(order="F")
#         except:
#             vtk['values'] = np.full((self.GRID.nx, self.GRID.ny, self.GRID.nz), 0).flatten(order="F")
#         p.add_mesh(vtk, show_edges=False)
#         bounds = p.show_bounds(mesh=vtk)
#         p.add_actor(bounds)

#     for (i, feature) in enumerate(features):
#         vtk = vtk.copy()
#         p.subplot(1, i)
#         p.add_text(feature, font_size=24)
#         try:
#             vtk['values'] = self.geology.geologic_features[feature].data.flatten(order="F")
#         except:
#             vtk['values'] = np.full((self.GRID.nx, self.GRID.ny, self.GRID.nz), 0).flatten(order="F")
#         slices = vtk.slice_orthogonal()
#         p.add_mesh(slices, show_edges=False)
#         bounds = p.show_bounds(mesh=vtk)
#         p.add_actor(bounds)

#     p.link_views()
#     p.show(cpos='xy')
#     return None