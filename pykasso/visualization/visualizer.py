"""
TODO
"""

# TODO
# - spécifier le type des retours de chaque fonction
# plt.figure() : options settable
# voxels : option extent

### Internal dependencies
import copy

### External dependencies
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from IPython import get_ipython

### Optional dependencies
try:
    import pyvista as pv
except ImportError:
    _has_pyvista = False
else:
    _has_pyvista = True

### Local dependencies
from pykasso.core._namespaces import ANISOTROPIC_FEATURES, DOMAIN_FEATURES

### Typing
from pykasso._typing import Project


def requires_pyvista():
    def _(function):
        def _wrapper(*args, **kwargs):
            if not _has_pyvista:
                raise ImportError("PyVista is required to do this.")
            result = function(*args, **kwargs)
            return result
        return _wrapper
    return _


class Visualizer():
    """
    This class manages pyKasso's project and provides methods to plot
    simulated karstic conduit networks.
    """
    
    DEFAULT_SETTINGS = {
        'surfaces': [],
        'show_inlets': False,
        'show_outlets': False,
        'show_colorbar': False
    }
    
    PV_DEFAULT_SETTINGS = {
        'ghosts': [],
        'ghost_values': [],
        'ghost_subdomains': [],
        'surfaces': {
            'topography': False,
            'water_level': False,
            'bedrock': False
        },
        'outline': False,
        'grid': False,
        # 'inlets': False,  # TODO - rename (see above)
        # 'outlets': False,  # TODO - rename (see above)
        'cpos': 'xz',
        # 'domain',
        'slice': False,
        'orientation': False,  # TODO - to remove ?
        # 'show_colorbar': False  # TODO - rename (see above)
    }
    
    MPL_DEFAULT_SETTINGS = {}
    
    def __init__(self, project: Project, notebook: bool = False,
                 *args, **kwargs) -> None:
        """
        Initialize the class.
        """
        self.project = project
        self._notebook = notebook
        
        # Set notebook value
        self._set_notebook_value(notebook)
        
        if _has_pyvista:
            # Set the pyvista grid
            self._set_pyvista_grid()
        
    @property
    def notebook(self):
        return self._notebook

    @notebook.setter
    def notebook(self, value):
        self._notebook = value
        self._set_notebook_value(value)
        
    def _set_notebook_value(self, boolean):
        """
        Render plots within Notebooks if 'boolean' is true.
        """
        # *** MATPLOTLIB *** #
        if boolean:
            get_ipython().run_line_magic('matplotlib', 'inline')
        else:
            get_ipython().run_line_magic('matplotlib', 'qt')
        
        # *** PYVISTA *** #
        if _has_pyvista:
            # Set global_theme.notebooks
            pv.global_theme.notebook = boolean
        
        return None
    
    # def _set_default_kwargs(self, kwargs) -> dict:
    #     for (key, value) in Visualizer.DEFAULT_SETTINGS.items():
    #         kwargs.setdefault(key, value)
    #     return kwargs
    
    # ************************************************************* #
    # ****          STATIC METHODS      *************************** #
    # ************************************************************* #
    
    # **** MPL **** #
    
    @staticmethod
    def mpl_plot_array_2D(array: np.ndarray, axis: str = 'z',
                          n_slice: int = -1, imshow_options: dict = {},
                          show_colorbar: bool = True, return_ax: bool = False):
        """
        TODO
        """
        # Slice the array
        if array.ndim == 2:
            array = array.T
        elif array.ndim == 3:
            array = Visualizer._slice_array(array, axis, n_slice)
        else:
            msg = "Array dimension is not valid."  # TODO - better error msg ?
            raise ValueError(msg)
        
        # Create the figure
        fig, ax = Visualizer._create_figure_2D(axis)
        
        # Plot the array
        im = ax.imshow(array, origin="lower", **imshow_options)
        
        # Set the colorbar
        if show_colorbar:
            fig.colorbar(im, ax=ax)
        
        if return_ax:
            return (fig, ax)
        else:
            return fig
    
    @staticmethod
    def _slice_array(array: np.ndarray, axis: str, n_slice: int) -> np.ndarray:
        """
        DOC
        """
        if axis == 'x':
            np_slice = (n_slice, slice(None), slice(None))
        elif axis == 'y':
            np_slice = (slice(None), n_slice, slice(None))
        elif axis == 'z':
            np_slice = (slice(None), slice(None), n_slice)
            
        sliced_data = array[np_slice].T
        
        return sliced_data
    
    @staticmethod
    def _create_figure_2D(axis: str) -> tuple:
        """
        Declare the corresponding 2 dimensional matplotlib figure
        according to the axis.
        """
        
        # Declare the figure
        fig = plt.figure()
        ax = plt.subplot()
        
        if axis == 'x':
            ax.set_xlabel('y')
            ax.set_ylabel('z')
        elif axis == 'y':
            ax.set_xlabel('x')
            ax.set_ylabel('z')
        elif axis == 'z':
            ax.set_xlabel('x')
            ax.set_ylabel('y')
        
        out = (fig, ax)
        return out
    
    @staticmethod
    def mpl_plot_array_3D(array: np.ndarray, voxels_options: dict = {},
                          return_ax: bool = False):
        """
        TODO
        """
        # Create the figure
        fig, ax = Visualizer._create_figure_3D()
        
        # Plot the array
        ax.voxels(array, **voxels_options)
        
        # Set the colors
        # TODO
        
        # Set aspect on 'equal'
        ax.set_aspect('equal')
        
        if return_ax:
            return (fig, ax)
        else:
            return fig
    
    @staticmethod
    def _create_figure_3D() -> tuple:
        """
        DOC
        """
        # Declare the figure
        fig = plt.figure()
        ax = plt.subplot(projection='3d')
        
        # Set axis labels
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        
        # Return the figure
        out = (fig, ax)
        return out
    
    # **** PV **** #
    
    @requires_pyvista()
    @staticmethod
    def pv_plot_array(array, ghost=False):
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
    
    # ************************************************************* #
    # ****          MATPLOTLIB          *************************** #
    # ************************************************************* #
    
    def mpl_plot_2D(self, feature: str, n_sim: int = -1, n_iteration: int = -1,
                    axis: str = 'z', n_slice: int = -1,
                    feature_options: dict = {}, show_colorbar: bool = True,
                    inlets_options: dict = None, outlets_options: dict = None):
        """
        TODO
        """
        # Retrieve the data
        sim_data = self.project._get_simulation_data(n_sim)
        
        # Retrieve the data feature
        array = self._get_feature_data(sim_data, feature, n_iteration)
        
        # Plot the feature data
        extent = self.project.grid.extent
        feature_options.setdefault('extent', extent)
        fig, ax = self.mpl_plot_array_2D(array, axis=axis, n_slice=n_slice,
                                         imshow_options=feature_options,
                                         show_colorbar=show_colorbar,
                                         return_ax=True)
        
        # Plot the inlets
        if inlets_options is not None:
            inlets = sim_data['inlets']
            x = inlets.x
            y = inlets.y
            ax.scatter(x, y, **inlets_options)
        
        # Plot the outlets
        if outlets_options is not None:
            outlets = sim_data['outlets']
            x = outlets.x
            y = outlets.y
            ax.scatter(x, y, **outlets_options)
        
        return fig
    
    def mpl_plot_3D(self, feature: str, n_sim: int = -1, n_iteration: int = -1,
                    voxels_options: dict = {}):
        """
        DOC
        """
        # Retrieve the data
        sim_data = self.project._get_simulation_data(n_sim)
        
        # Retrieve the data feature
        array = self._get_feature_data(sim_data, feature, n_iteration)
        
        # Plot the feature data
        fig, ax = self.mpl_plot_array_3D(array, voxels_options, return_ax=True)
        
        return fig
    
    # **** #
    
    def mpl_plot_karstic_network(self, n_sim: int = -1) -> None:
        """
        DOC
        """
        # Retrieve the data
        sim_data = self.project._get_simulation_data(n_sim)
        
        # Create the figure
        fig, ax = self._create_figure_3D()
        
        # Gets nodes
        nodes_ = sim_data['vectors']['nodes']
        nodes = {}
        for node in nodes_:
            nodes[node] = nodes_[node][:3]
            
        # Get edges
        edges = sim_data['vectors']['edges'].values()
        
        # Plot the karstic network
        ax = self._plot_graph(ax, nodes, edges)
        
        return fig
    
    def _plot_graph(self, ax, nodes: dict, edges: list):
        """
        DOC
        """
        
        # Plot edges
        nodes_to_plot = []
        for edge in edges:
            n_node1, n_node2 = edge
            nodes_to_plot.extend([n_node1, n_node2])
            node1 = nodes[n_node1]
            node2 = nodes[n_node2]
            x1, y1, z1 = node1
            x2, y2, z2 = node2
            xs = [x1, x2]
            ys = [y1, y2]
            zs = [z1, z2]
            ax.plot(xs, ys, zs)
        
        # Plot nodes
        nodes_list = []
        for node_to_plot in nodes_to_plot:
            nodes_list.append(nodes[node_to_plot])
        xs, ys, zs = zip(*nodes_list)
        ax.scatter(xs, ys, zs, marker='o')
        return ax
    
    # **** #
    
    def _get_feature_data(self, data: dict, feature: str,
                          iteration: str) -> np.ndarray:
        """
        DOC
        """
        if feature in ['geology', 'faults', 'fractures']:
            featured_data = data[feature].data_volume
        elif feature in ['karst']:
            featured_data = data['maps'][feature][iteration]
        elif feature in ['anisotropic']:  # TODO
            pass
        else:
            msg = 'feature keyword error'
            raise ValueError(msg)
        return featured_data
    
    # def show_karst_system(self, n_sim: int = -1,
    #                       settings: dict = {}) -> None:
        
    #     # Adds default settings values
    #     settings = self._get_mpl_default_settings(settings)
        
    #     # Retrieves the data
    #     simulation = self._get_simulation_data(n_sim)
        
    #     # Creates the plot
    #     # fig = plt.figure(figsize=None, dpi=None)  # TODO
    #     fig, ax = self._create_figure()
        
    #     # Plots the surfaces
    #     # surfs = []
    #     for surface in settings['surfaces']:
    #         x, y, z = self._get_mpl_surface(simulation, surface)
    #         ax.plot_surface(x, y, z, linewidth=0, antialiased=False)
        # surf = ax.plot_surface(x, y, z, cmap=cm.coolwarm,
        #                        linewidth=0, antialiased=False)
        
    #     # Plots the inlets
    #     if settings['inlets']:
    #         inlets = simulation['inlets']
    #         points = inlets[['x', 'y', 'z']].values
    #         xs, ys, zs = zip(*points)
    #         ax.scatter(xs, ys, zs, marker='o', c='r', s=40)
        
    #     # Plots the outlets
    #     if settings['outlets']:
    #         outlets = simulation['outlets']
    #         points = outlets[['x', 'y', 'z']].values
    #         xs, ys, zs = zip(*points)
    #         ax.scatter(xs, ys, zs, marker='o', c='b', s=40)
        
    #     edges = simulation['vectors']['edges']
    #     nodes = simulation['vectors']['nodes']
    #     nodes_list = list(nodes.values())
    #     xs, ys, zs, t = zip(*nodes_list)
    #     ax.scatter(xs, ys, zs, marker='o')

    #     return None
    
    # def _get_mpl_default_settings(self, settings):
    #     for key, value in PyplotVisualizer.DEFAULT_SETTINGS.items():
    #         if key not in settings:
    #             settings[key] = value
    #     return settings
    
    # def _get_mpl_surface(self, simulation, surface: str):
    #     grid = simulation['grid']
    #     x = grid.X[:, :, 0]
    #     y = grid.Y[:, :, 0]

    #     if surface == 'topography':
    #         z = simulation['domain'].topography.data_surface
    #     elif surface == 'water_level':
    #         z = simulation['domain'].water_level.data_surface
    #     elif surface == 'bedrock':
    #         z = simulation['domain'].bedrock.data_surface
        
    #     return (x, y, z)

    # ************************************************************* #
    # ****          PYVISTA             *************************** #
    # ************************************************************* #
    
    def _set_pyvista_grid(self) -> None:
        """
        DOC
        """
        # Initialisation
        p = self.project.grid.get_grid_parameters()
        x0, y0, z0 = p['x0'], p['y0'], p['z0']
        nx, ny, nz = p['nx'], p['ny'], p['nz']
        dx, dy, dz = p['dx'], p['dy'], p['dz']
        pv_grid = pv.UniformGrid()
        # Construct grid
        pv_grid.origin = (x0 - dx / 2, y0 - dy / 2, z0 - dz / 2)
        pv_grid.dimensions = np.array((nx, ny, nz)) + 1
        pv_grid.spacing = (dx, dy, dz)
        self.pv_grid = pv_grid.cast_to_unstructured_grid()
        return None
    
    @requires_pyvista()
    def pv_plot(self, n_sim, feature, settings={}, show=True,
                return_plotter=False):
        """
        DOC
        """
        ######################
        ### Initialisation ###
        ######################
        
        # Create plotter
        plotter = pv.Plotter()
        
        # Check 'settings' parameter
        settings = self._check_pv_settings(settings)
        if 'text' not in settings:
            text = 'Simulation {} - {}'.format(n_sim, feature)
            settings['text'] = text
            
        # Get simulation data
        simulation_data = self.project._get_simulation_data(n_sim)
        
        ###############
        ### Ploting ###
        ###############
        
        # Plot actors
        plotter = self._fill_plotter(plotter, simulation_data,
                                     feature, settings)
        
        # Plot figure
        if show:
            plotter.show()
        
        if return_plotter:
            return plotter
        else:
            return None
    
    def _check_pv_settings(self, settings):
        """
        DOC
        """
        for key, value in Visualizer.PV_DEFAULT_SETTINGS.items():
            if key not in settings:
                settings[key] = value
        for key, value in Visualizer.PV_DEFAULT_SETTINGS['surfaces'].items():
            if key not in settings['surfaces']:
                settings['surfaces'][key] = value
        return settings
    
    def _fill_plotter(self, plotter, simulation_data, feature, settings):
        """
        DOC
        """
        # Print text
        plotter.add_text(settings['text'], font_size=20)
        
        # Generate actors
        actors = self._get_actors(simulation_data, feature, settings)
        
        # Print actors
        for actor in actors:
            plotter.add_actor(actor, reset_camera=True)
            
        return plotter
    
    def _get_actors(self, simulation_data, feature, settings: dict):
        """
        DOC
        """
        
        # Get feature data
        iteration = -1  # TODO
        feature_data = self._get_feature_data_from_dict(simulation_data,
                                                        feature,
                                                        iteration)
        
        # Get mesh data
        data = self._get_mesh_from_feature(feature_data)
            
        # Ghost the data
        test_ghost_values = (len(settings['ghost_values']) > 0)
        test_ghost_subdomains = (len(settings['ghost_subdomains']) > 0)
        if test_ghost_values or test_ghost_subdomains:
            
            domain = simulation_data['domain']
            # If necessery, retrieves the subdomains
            if test_ghost_subdomains > 0:
                subdomains = []
                for subdomain in settings['ghost_subdomains']:
                    if subdomain[-2:] == '_r':
                        subdomain = subdomain[:-2]
                        data_domain = domain.get_subdomain(subdomain)
                        data_domain = np.logical_not(data_domain).astype(int)
                    else:
                        data_domain = domain.get_subdomain(subdomain)
                    subdomains.append(data_domain)
            else:
                subdomains = []
            
            data = self._ghost_data(data, settings['ghost_values'], subdomains)
        
        ### Creates the plot
        plotter = pv.Plotter()
        actors = []
        
        # Plots the outline of the domain
        if settings['outline']:
            _ = plotter.add_mesh(self.pv_grid.outline(), color="k")
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
            inlets = self._get_points(simulation_data['inlets'])
            _ = plotter.add_points(inlets, **kwargs)
            actors.append(_)
        if settings['outlets']:
            kwargs['color'] = 'b'
            outlets = self._get_points(simulation_data['outlets'])
            _ = plotter.add_points(outlets, **kwargs)
            actors.append(_)
            
        # Plots the surfaces
        kwargs = {
            'smooth_shading': True,
            'opacity': 0.66,
            'scalars': 'data',
        }
        for surface_name, boolean in settings['surfaces'].items():
            if boolean:
                surface = self._get_surface(simulation_data, surface_name)
                _ = plotter.add_mesh(surface, **kwargs)
                actors.append(_)
        
        # Plots the data
        kwargs = {
            'scalar_bar_args': {'title': 'Vol1'},
            'scalars': 'data',
        }
        if settings['slice']:
            slices = data.copy().slice_orthogonal()
            _ = plotter.add_mesh(slices, **kwargs)
        else:
            _ = plotter.add_mesh(data.copy(), **kwargs)
        actors.append(_)
        
        # Plots the scalar bar
        if settings['colorbar']:
            _ = plotter.add_scalar_bar()
            actors.append(_)

        return actors
    
    def _get_feature_data_from_dict(self, data_dict, feature, iteration=-1):
        """
        DOC
        """
        if feature in ANISOTROPIC_FEATURES:
            data = data_dict['maps'][feature][iteration]
        elif feature in DOMAIN_FEATURES:
            feature_object = getattr(data_dict['domain'], feature)
            data = getattr(feature_object, 'data_volume')
        else:
            feature_object = data_dict[feature]
            data = getattr(feature_object, 'data_volume')
        return data
        
    def _get_mesh_from_feature(self, data):
        """
        DOC
        """
        mesh = self.pv_grid.copy()
        # nx, ny, nz = data.shape
        # data = data.reshape(nx, nz, ny, order="F")
        data = data.flatten(order="F")
        # data = np.rot90(data, k=1
        # data = np.swapaxes(data, 1, 2)
        mesh.cell_data['data'] = data
        return mesh
    
    def _ghost_data(self, data, ghost_values: list, ghost_subdomains: list):
        """
        DOC
        """
        ghosted_cells = []
        
        # Ghost cells according to data value
        if len(ghost_values) > 0:
            test = np.isin(data["data"], ghost_values)
            ghosted_cells.extend(np.argwhere(test))
        
        # Ghost cells according to subdomain
        if len(ghost_subdomains) > 0:
            tests = []
            for subdomain in ghost_subdomains:
                tests.append(subdomain == 1)
            test_subdomains = np.logical_or.reduce(tests)
            data_domain = test_subdomains.flatten(order="F")
            data.cell_data['domain'] = data_domain
            ghosted_cells.extend(np.argwhere(np.isin(data["domain"], 1)))
            
        data = data.remove_cells(ghosted_cells)
        return data
    
    def _get_points(self, points):
        """
        DOC
        """
        labels = points.index.values.tolist()
        points = points[['x', 'y', 'z']].values
        points = points.astype('float32')
        cloud = pv.wrap(points)
        cloud['labels'] = labels
        return cloud
    
    def _get_surface(self, simulation_data, surface_name: str):
        X, Y, Z = self.project.grid.get_meshgrids()
        x = X[:, :, 0]
        y = Y[:, :, 0]
        
        surface_object = getattr(simulation_data['domain'], surface_name)
        try:
            z = surface_object.data_surface
            surface = pv.StructuredGrid(x, y, z)
            surface['data'] = z.flatten(order="F")
            return surface
        except AttributeError:
            msg = "Simulation has no '{}' surface.".format(surface_name)
            raise AttributeError(msg)
    
    @requires_pyvista()
    def pv_show(self, simulations: list, features: list,
                settings: list = [{}], savefig: str = None) -> None:
        """
        TODO
        """
        ######################
        ### Initialisation ###
        ######################
        
        ### Test parameters validity
        
        # Check type of 'settings' parameter
        if not isinstance(settings, (list)):
            settings = [settings] * len(simulations)
        
        ### Insert default settings values
        pv_settings = []
        for settings_dict in settings:
            settings_dict = self._check_pv_settings(settings_dict)
            pv_settings.append(copy.deepcopy(settings_dict))
        
        ### Create plotter
        if len(features) == 1:
            if len(simulations) <= 3:
                rows = 1
                columns = len(simulations)
            else:
                side = int(np.ceil(np.sqrt(len(simulations))))
                if side % len(simulations) == 0:
                    rows = side - 1
                else:
                    rows = side
                columns = side
        else:
            rows = len(simulations)
            columns = len(features)
            
        shape = (rows, columns)
        border = True
        off_screen = False
        if savefig is not None:
            off_screen = True
        plotter = pv.Plotter(shape=shape, border=border, off_screen=off_screen)
        
        ###############
        ### Ploting ###
        ###############
        
        # For each simulation, print the required plots
        if len(features) == 1:
            for (i, (n_sim, settings)) in enumerate(zip(simulations,
                                                        pv_settings)):
                row = int(np.floor(i / columns))
                column = i % columns
                settings = pv_settings[i]
                feature = features[0]
                
                if 'text' not in settings:
                    text = 'Simulation {} - {}'.format(n_sim, feature)
                    settings['text'] = text
                data = self.project._get_simulation_data(n_sim)
                plotter.subplot(row, column)
                plotter = self._fill_plotter(plotter, data,
                                             feature, settings)
        else:
            for (row, (n_sim, settings)) in enumerate(zip(simulations,
                                                          pv_settings)
                                                      ):
                for (column, feature) in enumerate(features):
                    if 'text' not in settings:
                        text = 'Simulation {} - {}'.format(n_sim, feature)
                        settings['text'] = text
                    data = self.project._get_simulation_data(n_sim)
                    plotter.subplot(row, column)
                    plotter = self._fill_plotter(plotter, data,
                                                 feature, settings)
                    
        plotter.link_views()
        
        if savefig is not None:  # TODO
            path = self.project_directory + '/outputs/' + savefig
            plotter.show(cpos='xz', screenshot=path)
        else:
            plotter.show(cpos='xz')

        return None
    
    
    #     def create_gif(self, simulation: int, feature: str, location: str,
#                    zoom: float = 1, ghosts: list = [], n_points: int = 24,
#                    fps: int = 10, window_size=[1024, 768],
#                    colormap: str = 'viridis',
#                    background_color: str = 'white',
#                    background_color_top: str = None) -> None:
#         """
#         TODO
        
#         colormap: str
#             https://matplotlib.org/stable/tutorials/colors/colormaps.html
#         """
        
#         ### Method based on those examples:
#         # https://docs.pyvista.org/examples/02-plot/orbit.html#orbiting
        
#         ### Gets the simulation data
#         simulation_data = self._get_simulation_data(simulation)
        
#         ### Gets the mesh
#         mesh = self._get_data_from_feature(simulation_data, feature)
#         # mesh_ = mesh.copy()
#         if len(ghosts) > 0:
#             mesh = self._ghost_data(mesh, ghosts)
        
#         ### Constructs the plotter
#         plotter = pv.Plotter(off_screen=True, window_size=window_size)
#         plotter.store_image = True
        
#         # Plots the data
#         kwargs = {
#             'cmap': colormap,
#             'scalar_bar_args': {'title': 'Vol1'},
#             'scalars': 'data',
#             'lighting': True,
#         }
#         plotter.add_mesh(mesh, **kwargs)
#         plotter.remove_scalar_bar()
        
#         # Sets background color
#         if background_color_top is None:
#             background_color_top = background_color
#         plotter.set_background(background_color)
        
#         ### Sets the initial camera position
#         plotter.camera.zoom(zoom)
#         # plotter.camera.roll = 0
#         plotter.camera.elevation = 0
#         plotter.camera.azimuth = 0
        
#         ### Generates the GIF
        
#         # Open a gif
#         plotter.open_gif(location, fps=fps)
        
#         # Creates the camera positions
#         azimuth_step = 360 / n_points
        
#         # Loops
#         for i in range(n_points):
            
#             # Updates camera position
#             # plotter.camera.roll
#             # plotter.camera.elevation
#             plotter.camera.azimuth += azimuth_step
            
#             # Updates the background
            
#             # Writes the frame in the gif
#             plotter.write_frame()
        
#         # Closes and finalizes movie
#         plotter.close()
#         return None
    
    
    
    
    
    
    
    
    



    
        
        












































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

































































# #     #############################
# #     ### Visualization methods ###
# #     #############################
# #
# #     def show_catchment(self, label='geology', title=None, cmap='binary'):
# #         """
# #         Show the entire study domain.
# #
# #         Parameters
# #         ----------
# #         label : str, optional
# #             Data to show : 'geology', 'topography', 'orientationx', 'orientationy' 'faults' or 'fractures'.
# #             By default : 'geology'.
# #         title : str, optional
# #             Title of the plot. If 'None', 'data' becomes the label.
# #         cmap : str, optional
# #             Color map, 'binary' by default.
# #         """
# #         import matplotlib.patches as mp
# #         fig, ax1 = plt.subplots()
# #         #if title is None:
# #         #    title = label
# #         #fig.suptitle(title, fontsize=16)
# #
# #         # Load data to show
# #         try:
# #             data = [data.data[:,:,0] for data in self.geology.data if data.label==label][-1]
# #         except:
# #             print('no data for indicated label parameter')
# #
# #         im1 = ax1.imshow(data.T, origin="lower", extent=self.grid.extent, cmap=cmap)
# #
# #         fig.colorbar(im1, ax=ax1)
# #         if self.settings['data_has_mask']:
# #             import matplotlib.patches as mp
# #             p1 = mp.PathPatch(self.mask.polygon, lw=2, fill=0, edgecolor='red', label='mask')
# #             ax1.add_patch(p1)
# #
# #         # Points
# #         for pts in self.points.points:
# #             x, y = zip(*pts.points)
# #             ax1.plot(x, y, 'o', label=pts.points_key)
# #
# #         ax1.set_aspect('equal', 'box')
# #         plt.legend(loc='upper right')
# #         plt.show()
# #         return fig
# #
# #     def _show_maps(self, sim=-1, iteration=-1, cmap='binary'):
# #         """
# #         Show the simulated karst network as an image.
# #         """
# #         karst_network = self.karst_simulations[sim]
# #
# #         fig, ([ax1,ax2],[ax3,ax4]) = plt.subplots(2, 2, sharex=True, sharey=True)
# #         fig.suptitle('Karst Network', fontsize=16)
# #
# #         ax1.imshow(karst_network.maps['outlets'], extent=self.grid.extent, origin='lower', cmap=cmap)
# #         ax1.set_title('Outlets')
# #
# #         ax2.imshow(karst_network.maps['cost'][iteration], extent=self.grid.extent, origin='lower', cmap=cmap)
# #         ax2.set_title('Cost')
# #
# #         ax3.imshow(karst_network.maps['time'][iteration], extent=self.grid.extent, origin='lower', cmap=cmap)
# #         ax3.set_title('Time')
# #
# #         ax4.imshow(karst_network.maps['karst'][iteration], extent=self.grid.extent, origin='lower', cmap=cmap)
# #         ax4.set_title('Karst')
# #
# #         fig.subplots_adjust(hspace=0.5)
# #         plt.show()
# #         return None
# #
# #
# #     def show(self, data=None, title=None):
# #         """
# #         Show the entire study domain (defaults to showing most recent simulation).
# #         """
# #         if data is None:
# #             data = self.karst_simulations[-1]
# #
# #         fig = plt.figure(figsize=(20,10))
# #
# #         # Cost map
# #         fig.add_subplot(131, aspect='equal')
# #         d = data.maps['cost'][-1]
# #         plt.xlabel('Cost array'+str(d.shape))
# #         d = np.transpose(d, (1,0)) # imshow read MxN and we have NxM
# #         plt.imshow(d, extent=self.grid.extent, origin='lower', cmap='gray') #darker=slower
# #         plt.colorbar(fraction=0.046, pad=0.04)
# #
# #         # Travel time map
# #         fig.add_subplot(132, aspect='equal')
# #         d = data.maps['time'][-1]
# #         plt.xlabel('Travel time array'+str(d.shape))
# #         d = np.transpose(d, (1,0)) # imshow read MxN and we have NxM
# #         plt.imshow(d, extent=self.grid.extent, origin='lower', cmap='cividis') #darker=faster
# #         plt.colorbar(fraction=0.046, pad=0.04)
# #
# #         # Karst map
# #         fig.add_subplot(133, aspect='equal')
# #         d = data.maps['karst'][-1]
# #         plt.xlabel('Karst array'+str(d.shape))
# #         d = np.transpose(d, (1,0)) # imshow read MxN and we have NxM
# #         plt.imshow(d, extent=self.grid.extent, origin='lower', cmap='gray_r') #darker=conduits
# #         plt.colorbar(fraction=0.046, pad=0.04)
# #         i = plt.scatter(data.points['inlets'].x,  data.points['inlets'].y,  c='orange')
# #         o = plt.scatter(data.points['outlets'].x, data.points['outlets'].y, c='steelblue')
# #         p = matplotlib.patches.Rectangle((0,0),0,0, ec='r', fc='none')
# #         if self.settings['data_has_mask']:
# #             closed_polygon = self.mask.vertices[:]
# #             closed_polygon.append(closed_polygon[0])
# #             x,y = zip(*closed_polygon)
# #             plt.plot(x,y, color='red', label='mask')
# #         #plt.legend([i,o,p], ['inlets', 'outlets', 'catchment'], loc='upper right')
# #         plt.legend([i,o], ['inlets', 'outlets'], loc='upper right')
# #
# #         if title is not None:
# #             fig.suptitle(title, fontsize=16)
# #         plt.show()
# #         return fig
# #
# #     def show_network(self, data=None, simplify=False, ax=None, plot_nodes=True, mask=True, labels=['inlets', 'outlets'], title=None, cmap=None, color='k', alpha=1, legend=True):
# #         """
# #         #Chloe: This is a new function that I use to create all the figures for the paper.
# #         Show the karst network as a graph with nodes and edges. Defaults to showing latest iteration.
# #
# #         Parameters
# #         ----------
# #         data:
# #             karst simulation object containing nodes, edges, points, etc. Can be obtained from self.karst_simulations[i]
# #         ax :
# #             axis to plot on
# #         label :
# #             None or list of strings ['nodes','edges','inlets','outlets'], indicating which components to label
# #         title : str
# #             title of plot
# #         cmap : str
# #             colormap to use when plotting
# #         color : str
# #             single color to use when plotting (cannot have both cmap and color)
# #         alpha : float
# #             opacity to plot with (1=opaque, 0=transparent)
# #         legend : bool
# #             whether to display legend
# #         plot_nodes : bool
# #             whether to display nodes
# #         polygon : bool
# #             whether to display the bounding polygon
# #         """
# #
# #         if ax == None:
# #             fig,ax = plt.subplots(figsize=(10,10))
# #             ax.set_aspect('equal')
# #
# #         if data == None:
# #             data = self.karst_simulations[-1]
# #
# #         if mask == True:
# #             if self.settings['data_has_mask']:
# #                 closed_polygon = self.mask.vertices[:]
# #                 closed_polygon.append(closed_polygon[0])
# #                 x,y = zip(*closed_polygon)
# #                 ax.plot(x,y, color='maroon')
# #                 p = matplotlib.lines.Line2D([0],[0], color='k')
# #
# #         if simplify == True:
# #             nodes = data.network['nodes']   #get all nodes
# #             nodes_simple = data.network['karstnet'].graph_simpl.nodes  #get indices of only the nodes in the simplified graph
# #             nodes_simple = {key: nodes[key] for key in nodes_simple}   #make df of only the nodes in the simplified graph, for plotting
# #             edges = data.network['edges']   #get all edges
# #             edges_simple = data.network['karstnet'].graph_simpl.edges  #get only the edges in the simplified graph
# #             edges_simple = {i: edge for i,edge in enumerate(edges_simple)}   #make df of only the edges in the simplified graph, for p
# #             nodes = pd.DataFrame.from_dict(nodes_simple, orient='index', columns=['x','y','type']) #convert to pandas for easier plotting
# #             edges = pd.DataFrame.from_dict(edges_simple, orient='index', columns=['inNode','outNode'])
# #         else:
# #             nodes = pd.DataFrame.from_dict(data.network['nodes'], orient='index', columns=['x','y','type']) #convert to pandas for easier plotting
# #             edges = pd.DataFrame.from_dict(data.network['edges'], orient='index', columns=['inNode','outNode'])
# #
# #         #Set up data for plotting:
# #         fromX = nodes.x.loc[edges.inNode]      #calculate coordinates for link start and end points
# #         fromY = nodes.y.loc[edges.inNode]
# #         toX   = nodes.x.loc[edges.outNode]
# #         toY   = nodes.y.loc[edges.outNode]
# #
# #         #Plot nodes and edges:
# #         if plot_nodes:
# #             n = ax.scatter(nodes.x,              nodes.y,                  c='k',         alpha=alpha, s=5)  #scatterplot nodes
# #         i = ax.scatter(data.points['inlets'].x,  data.points['inlets'].y,  c='orange',    s=30) #scatterplot inlets
# #         o = ax.scatter(data.points['outlets'].x, data.points['outlets'].y, c='steelblue', s=30) #scatterplot outlets
# #         e = matplotlib.lines.Line2D([0],[0])                                                  #line artist for legend
# #         for ind in edges.index:                                                               #loop over edge indices
# #             if cmap is not None:
# #                 ax.plot((fromX.iloc[ind], toX.iloc[ind]), (fromY.iloc[ind], toY.iloc[ind]), c=plt.cm.get_cmap(cmap)(ind/len(edges)), alpha=alpha)  #plot each edge, moving along color gradient to show order
# #             elif color is not None:
# #                 ax.plot((fromX.iloc[ind], toX.iloc[ind]), (fromY.iloc[ind], toY.iloc[ind]), c=color, alpha=alpha)  #plot each edge in same color
# #
# #         #Add labels:
# #         if labels == None:
# #             pass
# #         else:
# #             if 'nodes' in labels:                                         #label node indices
# #                 for ind in nodes.index:                                   #loop over node indices
# #                     ax.annotate(str(ind), xy=(nodes.y[ind]-10, nodes.x[ind]))  #annotate slightly to left of each node
# #             if 'edges' in labels:
# #                 for ind in edges.index:
# #                     ax.annotate(str(ind), xy=(edges.y[ind]-10, edges.x[ind]))  #annotate slightly to left of each edge
# #             if 'inlets' in labels:
# #                 for index,inlet in data.points['inlets'].iterrows():
# #                     ax.annotate(str(int(inlet.outlet))+'-'+str(int(inlet.inlet_iteration)),  xy=(inlet.x-(6*self.grid.dx),  inlet.y))
# #             if 'outlets' in labels:
# #                 for index,outlet in data.points['outlets'].iterrows():
# #                     ax.annotate(str(int(outlet.name)), xy=(outlet.x-(4*self.grid.dx), outlet.y))
# #
# #         #Add legend & title:
# #         if legend:
# #             if plot_nodes:
# #                 if plot_polygon:
# #                     ax.legend([i,o,n,e,p],['inlets','outlets','nodes','edges','mask'])
# #                 else:
# #                     ax.legend([i,o,n,e],['inlets','outlets','nodes','edges'])
# #             else:
# #                 if plot_polygon:
# #                     ax.legend([i,o,e,p],['inlets','outlets','edges','mask'])
# #                 else:
# #                     ax.legend([i,o,e],['inlets','outlets','edges','mask'])
# #         if title is not None:
# #             ax.set_title(title, fontsize=16)
# #
# #         return None



# #     def show_catchment(self, label='geology', title=None, cmap='binary'):
# #         """
# #         Show the entire study domain.
# #
# #         Parameters
# #         ----------
# #         label : str, optional
# #             Data to show : 'geology', 'topography', 'orientationx', 'orientationy' 'faults' or 'fractures'.
# #             By default : 'geology'.
# #         title : str, optional
# #             Title of the plot. If 'None', 'data' becomes the label.
# #         cmap : str, optional
# #             Color map, 'binary' by default.
# #         """
# #         import matplotlib.patches as mp
# #         fig, ax1 = plt.subplots()
# #         #if title is None:
# #         #    title = label
# #         #fig.suptitle(title, fontsize=16)
# #
# #         # Load data to show
# #         try:
# #             data = [data.data[:,:,0] for data in self.geology.data if data.label==label][-1]
# #         except:
# #             print('no data for indicated label parameter')
# #
# #         im1 = ax1.imshow(data.T, origin="lower", extent=self.grid.extent, cmap=cmap)
# #
# #         fig.colorbar(im1, ax=ax1)
# #         if self.settings['data_has_mask']:
# #             import matplotlib.patches as mp
# #             p1 = mp.PathPatch(self.mask.polygon, lw=2, fill=0, edgecolor='red', label='mask')
# #             ax1.add_patch(p1)
# #
# #         # Points
# #         for pts in self.points.points:
# #             x, y = zip(*pts.points)
# #             ax1.plot(x, y, 'o', label=pts.points_key)
# #
# #         ax1.set_aspect('equal', 'box')
# #         plt.legend(loc='upper right')
# #         plt.show()
# #         return fig