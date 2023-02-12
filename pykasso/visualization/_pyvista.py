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

def _show_data(environment, feature, settings, show=True):
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
    
    domain_features = ['delimitation', 'topography', 'bedrock', 'water_level']

    if 'iteration' not in settings:
        if hasattr(environment, 'iteration'):
            settings['iteration'] = environment.iteration

    for attribute in attributes:
        if attribute not in settings:
            settings[attribute] = False

    #########################   
    ### Geologic features ###
    #########################

    # Get the grid
    grid_i = _get_grid(environment.grid)
    grid = grid_i

    # Get the data
    # if feature != 'grid':
    if feature in ['cost', 'alpha', 'beta', 'time', 'karst']:
        grid = _get_data_from_dict(grid, getattr(environment, 'maps'), feature, 'data', iteration=settings['iteration'])
    elif feature in domain_features:
        domain = environment.domain
        grid = _get_data_from_attribute(grid, getattr(domain, feature), 'data_volume', 'data')
    else:
        grid = _get_data_from_attribute(grid, getattr(environment, feature), 'data_volume', 'data')
        
    kwargs['scalars'] = 'data'

    # Get the domain
    if getattr(environment, 'domain') is not None:
        grid = _get_data_from_attribute(grid, getattr(environment, 'domain'), 'data_volume', 'domain')

    # Ghost the data
    if 'ghost' in settings:
        if 'domain' in settings:
            ghosts = np.argwhere(np.logical_or(np.isin(grid["data"], settings['ghost']), (grid["domain"] == 0)))
        else:
            ghosts = np.argwhere(np.isin(grid["data"], settings['ghost']))
    else:
        if 'domain' in settings:
            ghosts = np.argwhere(grid["domain"] == 0)
        else:
            ghosts = None
    
    if ghosts is not None:
        grid = grid.remove_cells(ghosts)

    ##############
    ### Points ###
    ##############

    if settings['inlets']:
        inlets = _get_points(environment.inlets)

    if settings['outlets']:
        outlets = _get_points(environment.outlets)

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
    plotter.add_mesh(grid_i.outline(), color="k")

    if settings['show_grid']:
        plotter.show_grid()

    # Points
    if settings['inlets']:
        plotter.add_points(inlets, render_points_as_spheres=False, point_size=20, color='r')
    if settings['outlets']:
        plotter.add_points(outlets, render_points_as_spheres=False, point_size=20, color='b')

    if show:
        plotter.show(cpos='xy')

    return _


##########################################################
### DEBUG ###
#############

def _debug_plot_model(environment, settings):
    """
    TODO
    """
    ### Call the plotter
    plotter = pv.Plotter(shape=(2, 4), border=True)

    features = ['geology', 'faults', 'fractures', 'beddings']

    for (i, feature) in enumerate(features):
        plotter.subplot(0, i)
        plotter.add_text(feature, font_size=24)

        if hasattr(environment, feature) and (getattr(environment, feature) is not None):
            actor = _show_data(environment, feature, {}, show=False)
            plotter.add_actor(actor, reset_camera=True)

            plotter.subplot(1, i)
            actor, misc = _show_data(environment, feature, {'slice':True}, show=False)
            plotter.add_actor(actor, reset_camera=True)

    plotter.link_views()
    plotter.show()
    return None


def _debug_plot_fmm(environment, settings):
    """
    TODO
    """
    ### Initializes ...
    if 'iterations' not in settings:
        settings['iterations'] = [0]

    if environment.fmm['algorithm'] == 'Isotropic3':
        features = ['cost', 'time', 'karst']
    elif environment.fmm['algorithm'] == 'Riemann3':
        features = ['cost', 'alpha', 'beta', 'time', 'karst']

    ### Call the plotter
    row = len(settings['iterations'])
    col = len(features)
    plotter = pv.Plotter(shape=(row, col), border=True)

    if hasattr(environment, 'maps'):
        for (i, feature) in enumerate(features):

            for j, iteration in enumerate(settings['iterations']):
            
                plotter.subplot(j, i)
                
                text = feature + ' - iteration : {}'.format(iteration)
                plotter.add_text(text, font_size=10)

                if feature == 'karst':
                    settings_ = {'iteration' : iteration, 'ghost' : [0]}
                else:
                    settings_ = {'iteration' : iteration}

                actor = _show_data(environment, feature, settings_, show=False)
                plotter.add_actor(actor, reset_camera=True)

    plotter.link_views()
    plotter.show()

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
    # mesh            = mesh.flip_z(inplace=False)
    mesh            = mesh.cast_to_unstructured_grid()
    return mesh

def _get_data_from_attribute(mesh, geologic_feature, attribute, label):
    """
    TODO
    """  
    data = getattr(geologic_feature, attribute)
    mesh.cell_data[label] = data.flatten(order="F")
    return mesh

def _get_data_from_dict(mesh, dict, attribute, label, iteration):
    """
    TODO
    """
    data = dict[attribute][iteration]
    mesh.cell_data[label] = data.flatten(order="F")
    
    return mesh

def _get_points(points):
    """
    TODO
    """
    points = points[['x', 'y', 'z']].values
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


#################################################################################

def _show_array(array):
    """
    TODO
    """
    if len(array.shape) == 2:
        nx, ny = array.shape
        nz = 1
    elif len(array.shape) == 3:
        nx, ny, nz = array.shape
    else:
        msg = "TODO"
        raise ValueError(msg)
    
    mesh = pv.UniformGrid()
    mesh.dimensions = np.array((nx, ny, nz)) + 1
    mesh.cell_data['data'] = array.flatten(order="F")
    mesh = mesh.cast_to_unstructured_grid()
    
    ghosts = np.argwhere(np.isin(mesh["data"], [0]))
    if ghosts is not None:
        mesh = mesh.remove_cells(ghosts)
        
    plotter = pv.Plotter()
    _ = plotter.add_mesh(mesh)
    # plotter.add_mesh(mesh.outline(), color="k")
    plotter.show(cpos='xy')
    
    return None