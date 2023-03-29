"""
TODO
"""

# TODO
# - incorporer cette option
# - 'show_edges' : True
# - 'cmap' : "viridis",

### External dependencies
import numpy as np
import pyvista as pv

### Local dependencies
from ._pyplot import PyplotVisualizer
from .main import DOMAIN_FEATURES, ANISOTROPIC_FEATURES

### Typing
# pyvista grid
# actors


class PyvistaVisualizer(PyplotVisualizer):
    """
    TODO
    """
    DEFAULT_SETTINGS = {
        'ghosts': [],
        'surfaces': [],
        'outline': True,
        'grid': False,
        'inlets': False,
        'outlets': False,
        'cpos': 'xz',
        # 'domain',
        'slice': False,
        # 'tracers',
        # other,
        'orientation': False
    }
    
    def __init__(self, project_directory: str, *args, **kwargs):
        """
        TODO
        """
        super().__init__(project_directory, *args, **kwargs)
        
        # Sets global_theme.notebooks
        # pv.global_theme.notebook = True
        
        # Sets the pyvista grid
        self._set_pyvista_grid()

    def _set_pyvista_grid(self) -> None:
        # Initialization
        x0, y0, z0 = self._get_grid_origin()
        nx, ny, nz = self._get_grid_dimensions()
        dx, dy, dz = self._get_grid_spacing()
        grid = pv.UniformGrid()
        # Constructs grid
        grid.origin = (x0 - dx / 2, y0 - dy / 2, z0 - dz / 2)
        grid.dimensions = np.array((nx, ny, nz)) + 1
        grid.spacing = (dx, dy, dz)
        self.grid = grid.cast_to_unstructured_grid()
        return None
    
    def show(self, simulations: list, features: list,
             settings: list = [{}]) -> None:
        """
        TODO
        """
        # Updates the project state
        self._update_project_state()
        
        # Checks settings type
        if not isinstance(settings, (list)):
            settings = [settings] * len(simulations)
        
        # Adds default settings values
        pyvista_settings = []
        for settings_ in settings:
            for key, value in PyvistaVisualizer.DEFAULT_SETTINGS.items():
                if key not in settings_:
                    settings_[key] = value
            pyvista_settings.append(settings_)
        
        # # Sets cpos
        # nx, ny, nz = self._get_grid_dimensions()
        # if nx == 1:
        #     pyvista_settings['cpos'] = 'yz'
        # elif ny == 1:
        #     pyvista_settings['cpos'] = 'xz'
        # elif nz == 1:
        #     pyvista_settings['cpos'] = 'xy'
        # else:
        #     pyvista_settings['cpos'] = 'xz'
        
        # Creates plotter
        if len(features) == 1:
            if len(simulations) <= 3:
                rows = 1
                columns = len(simulations)
            else:
                side = int(np.ceil(np.sqrt(len(simulations))))
                rows = side
                columns = side
        else:
            rows = len(features)
            columns = len(simulations)
        shape = (rows, columns)
        border = True
        plotter = pv.Plotter(shape=shape, border=border)
        
        # For each simulation, prints the required plots
        
        if len(features) == 1:
            for row in range(rows):
                for column in range(columns):
                    i = row * rows + column - 1
                    n_sim = simulations[i]
                    feature = features[0]
                    settings = pyvista_settings[i]
                    plotter.subplot(row, column)
                    plotter = self._fill_plotter(plotter, n_sim,
                                                 feature, settings)
        else:
            for (row, (n_sim, settings)) in enumerate(zip(simulations, pyvista_settings)):
                for (column, feature) in enumerate(features):
                    plotter.subplot(row, column)
                    plotter = self._fill_plotter(plotter, n_sim,
                                                 feature, settings)
                    
        plotter.link_views()
        plotter.show(cpos='xz')

        return None
    
    def _fill_plotter(self, plotter, n_sim, feature, settings):
        """ """
        text = 'Simulation {} - {}'.format(n_sim, feature)
        plotter.add_text(text, font_size=20)
        
        simulation_data = self._get_simulation_data(n_sim)
        
        actors = self._show_feature(simulation_data, feature, settings)
        for actor in actors:
            plotter.add_actor(actor, reset_camera=True)
            
        return plotter
    
    def create_gif(self, simulation: int, feature: str, location: str,
                   zoom: float = 1, ghosts: list = [], n_points: int = 24,
                   fps: int = 10, window_size=[1024, 768],
                   background_color: str = 'white') -> None:
        """
        TODO
        """
        
        ### Method based on those examples:
        # https://docs.pyvista.org/examples/02-plot/orbit.html#orbiting
        
        ### Gets the simulation data
        simulation_data = self._get_simulation_data(simulation)
        
        ### Gets the mesh
        mesh = self._get_data_from_feature(simulation_data, feature)
        if len(ghosts) > 0:
            mesh = self._get_ghosted_data(mesh, ghosts)
        
        ### Constructs the plotter
        plotter = pv.Plotter(off_screen=True, window_size=window_size)
        plotter.store_image = True
        plotter.add_mesh(mesh, lighting=False)
        plotter.remove_scalar_bar()
        plotter.camera.zoom(zoom)
        plotter.set_background(background_color)
        
        ### Generates the GIF
        path = plotter.generate_orbital_path(n_points=n_points,
                                             shift=mesh.length)
        plotter.open_gif(location, fps=fps)
        plotter.orbit_on_path(path, write_frames=True)
        plotter.close()
        
        return None
    
    def _get_simulation_data(self, n: int) -> dict:
        simulation_directory = self.project_state['simulation_locations'][n]
        simulation_data_path = simulation_directory + 'results.pickle'
        simulation_data = self._read_pickle(simulation_data_path)
        return simulation_data
        
    def _show_feature(self, simulation, feature: str, settings: dict = {}):
        
        # Gets the data
        feature_data = self._get_data_from_feature(simulation, feature)
        feature_data_ = feature_data.copy()
        
        # Gets the domain
        pass
        
        # Ghosts the data
        if len(settings['ghosts']) > 0:
            feature_data = self._get_ghosted_data(feature_data,
                                                  settings['ghosts'])
        
        # Cuts the data with model limits
        pass
        
        ### Creates the plot
        plotter = pv.Plotter()
        actors = []
        
        # Plots the outline of the domain
        if settings['outline']:
            _ = plotter.add_mesh(self.grid.outline(), color="k")
            actors.append(_)
            
        # Plots the grid of the domain
        if settings['grid']:
            _ = plotter.show_grid()
            actors.append(_)
        
        # Plots the points
        kwargs = {
            'render_points_as_spheres': False,
            'point_size': 20,
        }
        if settings['inlets']:
            kwargs['color'] = 'r'
            inlets = self._get_points(simulation['inlets'])
            _ = plotter.add_points(inlets, **kwargs)
            actors.append(_)
        if settings['outlets']:
            kwargs['color'] = 'b'
            outlets = self._get_points(simulation['outlets'])
            _ = plotter.add_points(outlets, **kwargs)
            actors.append(_)
            
        # Plots the surfaces
        kwargs = {
            'smooth_shading': True,
            'opacity': 0.66,
            'scalars': 'data',
        }
        if len(settings['surfaces']) > 0:
            for surface in settings['surfaces']:
                surface = self._get_surface(simulation, surface)
                _ = plotter.add_mesh(surface, **kwargs)
                actors.append(_)
        
        # Plots the data
        kwargs = {
            'scalar_bar_args': {'title': 'Vol1'},
            'scalars': 'data',
        }
        if settings['slice']:
            slices = feature_data.copy().slice_orthogonal()
            _ = plotter.add_mesh(slices, **kwargs)
        else:
            _ = plotter.add_mesh(feature_data.copy(), **kwargs)
        actors.append(_)
        
        # TODO - Plots the orientation arrows
        if settings['orientation']:
            kwargs = {
                'scalar_bar_args': {'title': 'Orientation'},
                'lighting': False,
            }
            nx, ny, nz = self._get_grid_dimensions()
            nx = nx + 1
            ny = ny + 1
            nz = nz + 1
            vx, vy, vz = simulation['orientation']
            vx = np.resize(vx, (nx, ny, nz))
            vy = np.resize(vy, (nx, ny, nz))
            vz = np.resize(vz, (nx, ny, nz))
            vectors = np.column_stack((vx.ravel(), vy.ravel(), vz.ravel()))
            feature_data_['orientation'] = vectors
            feature_data_.set_active_vectors("orientation")
            
            # contours = feature_data_.contour(8, scalars="data")
            # arrows = contours.glyph(orient="orientation")# factor=200.0)
            _ = plotter.add_mesh(feature_data_.arrows) #, **kwargs)
            actors.append(_)
        
        # Plots the scalar bar
        _ = plotter.add_scalar_bar()
        actors.append(_)

        return actors
  
    def _get_data_from_feature(self, simulation, feature: str,
                               iteration: int = -1):
        mesh = self.grid.copy()
        
        if feature in ANISOTROPIC_FEATURES:
            data = simulation['maps'][feature][iteration]
        elif feature in DOMAIN_FEATURES:
            feature_object = getattr(simulation['domain'], feature)
            data = getattr(feature_object, 'data_volume')
        else:
            feature_object = simulation[feature]
            data = getattr(feature_object, 'data_volume')
            
        mesh.cell_data['data'] = data.flatten(order="F")
        return mesh
    
    def _get_ghosted_data(self, data, ghosts: list):
        ghosted_cells = np.argwhere(np.isin(data["data"], ghosts))
        data = data.remove_cells(ghosted_cells)
        return data
    
    def _get_surface(self, simulation, surface: str):
        grid = simulation['grid']
        x = grid.X[:, :, 0]
        y = grid.Y[:, :, 0]
        
        if surface == 'topography':
            z = simulation['domain'].topography.data_surface
        elif surface == 'water_level':
            z = simulation['domain'].water_level.data_surface
        elif surface == 'bedrock':
            z = simulation['domain'].bedrock.data_surface

        surf = pv.StructuredGrid(x, y, z)
        surf['data'] = z.flatten(order="F")

        return surf

    def _get_points(self, points):
        labels = points.index.values.tolist()
        points = points[['x', 'y', 'z']].values
        points = points.astype('float32')
        cloud = pv.wrap(points)
        cloud['labels'] = labels
        return cloud
    

def _show_array(array, ghost=False):
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

    if ghost:
        ghosts = np.argwhere(np.isin(mesh["data"], [0]))
        mesh = mesh.remove_cells(ghosts)

    plotter = pv.Plotter()
    plotter.show_grid()
    kwargs = {}
    kwargs['scalars'] = 'data'
    _ = plotter.add_mesh(mesh, **kwargs)
    # plotter.add_mesh(mesh.outline(), color="k")
    plotter.show(cpos='xy')

    return None


# def show_simulation(self):
#         """
#         TODO
#         """
        # def _show_data(environment, feature, settings):
        
        #     # Gets the domain
        #     if hasattr(environment, 'domain') and (getattr(environment, 'domain') is not None):
        #         grid = _get_data_from_attribute(grid, getattr(environment, 'domain'), 'data_volume', 'domain')

        #     # Ghost the data
        #     if 'ghost' in settings:
        #         if settings['domain'] == True:
        #             ghosts = np.argwhere(np.logical_or(np.isin(grid["data"], settings['ghost']), (grid["domain"] == 0)))
        #         else:
        #             ghosts = np.argwhere(np.isin(grid["data"], settings['ghost']))
        #     else:
        #         if settings['domain'] == True:
        #             ghosts = np.argwhere(grid["domain"] == 0)
        #         else:
        #             ghosts = None
            
        #     if ghosts is not None:
        #         grid = grid.remove_cells(ghosts)

        # return None








# ##########################################################
# ### DEBUG ###
# #############

# def _debug_plot_model(environment, settings):
#     """
#     TODO
#     """
#     ### Call the plotter
#     plotter = pv.Plotter(shape=(2, 4), border=True)

#     features = ['geology', 'faults', 'fractures', 'beddings']

#     for (i, feature) in enumerate(features):
#         plotter.subplot(0, i)
#         plotter.add_text(feature, font_size=24)

#         if hasattr(environment, feature) and (getattr(environment, feature) is not None):
#             actor = _show_data(environment, feature, {}, show=False)
#             plotter.add_actor(actor, reset_camera=True)

#             plotter.subplot(1, i)
#             actor, misc = _show_data(environment, feature, {'slice':True}, show=False)
#             plotter.add_actor(actor, reset_camera=True)

#     plotter.link_views()
#     plotter.show()
#     return None


# def _debug_plot_fmm(environment, settings):
#     """
#     TODO
#     """
#     ### Initializes ...
#     if 'iterations' not in settings:
#         settings['iterations'] = [0]

#     if environment.fmm['algorithm'] == 'Isotropic3':
#         features = ['cost', 'time', 'karst']
#     elif environment.fmm['algorithm'] == 'Riemann3':
#         features = ['cost', 'alpha', 'beta', 'time', 'karst']

#     ### Call the plotter
#     row = len(settings['iterations'])
#     col = len(features)
#     plotter = pv.Plotter(shape=(row, col), border=True)

#     if hasattr(environment, 'maps'):
#         for (i, feature) in enumerate(features):

#             for j, iteration in enumerate(settings['iterations']):
            
#                 plotter.subplot(j, i)
                
#                 text = feature + ' - iteration : {}'.format(iteration)
#                 plotter.add_text(text, font_size=10)

#                 if feature == 'karst':
#                     settings_ = {'iteration' : iteration, 'ghost' : [0]}
#                 else:
#                     settings_ = {'iteration' : iteration}

#                 actor = _show_data(environment, feature, settings_, show=False)
#                 plotter.add_actor(actor, reset_camera=True)

#     plotter.link_views()
#     plotter.show()

#     return None

# ######################################################################

# # FRACTURES = domain.FRACTURES.location
# # print('Nbr fractures:', len(FRACTURES))
# # POLYGONS = []
# # for fracture in FRACTURES:
# #     x, y, z = fracture.get_position()
# #     a, b, c = fracture.get_normal()
# #     rad     = fracture.radius
# #     POLYGONS.append(pv.Polygon(center=(x, y, z), radius=rad, normal=(a, b, c), n_sides=definition))

#     # return None


###############################################################################

# TODO - meilleur ghosting





# def show_average_paths(self):
#     """
#     todo
#     """
#     ### Call the plotter
#     p = pv.Plotter(notebook=False)

#     ### Construct the grid
#     vtk = pv.UniformGrid()
#     vtk.dimensions = np.array((self.GRID.nx, self.GRID.ny, self.GRID.nz)) + 1
#     vtk.origin     = (self.GRID.x0 - self.GRID.dx/2, self.GRID.y0 - self.GRID.dy/2, self.GRID.z0 - self.GRID.dz/2)
#     vtk.spacing    = (self.GRID.dx, self.GRID.dy, self.GRID.dz)

#     vtk['values'] = self.karst_prob.flatten(order="F")

#     mesh = vtk.cast_to_unstructured_grid()
#     ghosts = np.argwhere(vtk['values'] < 1.0)
#     mesh.remove_cells(ghosts)
#     p.add_mesh(mesh, show_edges=False)

#     ### Plotting
#     # p.add_title(feature)
#     p.add_axes()
#     bounds = p.show_bounds(mesh=vtk)
#     p.add_actor(bounds)
#     p.show(cpos='xy')

#     return None