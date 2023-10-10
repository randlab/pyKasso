"""
Module defining a tool for pyKasso results visulization.
"""

### Internal dependencies
import copy
import warnings

### External dependencies
import numpy as np
import pandas as pd
import matplotlib.figure
import matplotlib.axes
import matplotlib.colorbar
import matplotlib.pyplot as plt
from IPython import get_ipython

### Optional dependencies
try:
    import pyvista as pv
    import pyvista.core
except ImportError:
    _has_pyvista = False
else:
    _has_pyvista = True

### Local dependencies
from pykasso.core._namespaces import (GEOLOGICAL_FEATURES,
                                      DOMAIN_FEATURES,
                                      ANISOTROPIC_FEATURES,
                                      AUTHORIZED_FEATURES,
                                      SURFACE_FEATURES)

### Typing
from typing import Union
from typing import NamedTuple

from pykasso.core.project import Project


class Figure(NamedTuple):
    """
    `NamedTuple` storing element from a matplotlib figure.
    
    Attributes
    ----------
    fig
        TODO
    ax
        TODO
    cbar
        TODO
    """
    fig: matplotlib.figure.Figure
    ax: matplotlib.axes._axes.Axes
    cbar: matplotlib.colorbar.Colorbar


def requires_pyvista():
    """
    If 'PyVista' package is not installed, return `ImportError` exception
    when a method requiring 'PyVista' is called.
    """
    def _(function):
        def _wrapper(*args, **kwargs):
            if not _has_pyvista:
                msg = ("PyVista is required to do this. On console : 'pip "
                       "install pyvista'.")
                raise ImportError(msg)
            result = function(*args, **kwargs)
            return result
        return _wrapper
    return _


class Visualizer():
    """
    This class manages pyKasso's project and provides methods to plot
    simulated karstic conduit networks.
    """
    
    DEFAULT_PV_SETTINGS = {
        'n_iteration': -1,
        'text_options': {},
        'ghost_values': [],
        'ghost_subdomains': [],
        'show_grid': True,
        'show_outline': False,
        'data_options': {},
        'surfaces_options': {},
        'inlets_options': None,
        'outlets_options': None,
        'show_slice': False,
        'show_colorbar': True,
        'fractures_options': None,
        # 'cpos': 'xz',
    }
    
    DEFAULT_MPL_IMSHOW_ORIGIN = 'lower'
    DEFAULT_MPL_IMSHOW_EXTENT = None
    PV_SCALAR_KEYWORD = 'data'
    
    def __init__(self,
                 project: Project,
                 notebook: bool = False,
                 ) -> None:
        """
        Initialize the class.

        Parameters
        ----------
        project : Project
            pyKasso project to read.
        notebook : bool, optional
            If `True`, figures are rendered within python notebooks,
            by default `False`
        """
        self.project = project
        self._notebook = notebook
        
        # Set notebook value
        self._set_notebook_value(notebook)
        
        # Set the pyvista grid
        if _has_pyvista:
            self._set_pyvista_grid()
        
    @property
    def notebook(self) -> bool:
        """
        Get the state of the rendering aspect. If `True`, figures are
        rendered within python notebooks.
        """
        return self._notebook

    @notebook.setter
    def notebook(self,
                 boolean: bool,
                 ) -> None:
        """
        When `notebook` is set, modify the way the figures are rendered.
        """
        self._notebook = boolean
        self._set_notebook_value(boolean)
        
    def _set_notebook_value(self,
                            boolean: bool,
                            ) -> None:
        """
        Render plots within python notebooks if `boolean` is `True`.
        """
        # *** MATPLOTLIB *** #
        if boolean:
            get_ipython().run_line_magic('matplotlib', 'inline')
        else:
            get_ipython().run_line_magic('matplotlib', 'qt')
        
        # *** PYVISTA *** #
        if _has_pyvista:
            pv.global_theme.notebook = boolean
        
        return None
    
    # ************************************************************* #
    # ****       Plot Array - Static Methods       **************** #
    # ************************************************************* #
    
    # **** MATPLOTLIB **** #
    
    @staticmethod
    def mpl_plot_array_2D(array: np.ndarray,
                          axis: str = 'z',
                          n_slice: int = -1,
                          imshow_options: dict = {},
                          contour_options: dict = None,
                          show_colorbar: bool = True,
                          ) -> Figure:
        """
        Plot a slice of a 3D array with the matplotlib library.

        Parameters
        ----------
        array : np.ndarray
            3D array to slice and plot.
        axis : str, optional
            Axis used to slice the array, by default 'z'.
        n_slice : int, optional
            Indice for slicing, by default -1.
        imshow_options : dict, optional
            Dictionary containing the arguments for the imshow() matplotlib
            function, by default {}.
            More details here:
            https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.imshow.html
        contour_options : dict, optional
            Dictionary containing the arguments for the contour() matplotlib
            function, by default None.
            More details here:
            https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.contour.html
        show_colorbar : bool, optional
            If `False`, remove the colorbar, by default `True`.
            More details here:
            https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.colorbar.html

        Returns
        -------
        Figure
            `NamedTuple` describing a matplotlib figure.
            - fig: TODO;
            - ax: TODO;
            - cbar: TODO;
            
        Examples
        --------
        >>> TODO
        """
        # Control array dimension
        Visualizer._test_array_dimension(array)
        
        # Slice the array
        if array.ndim == 2:
            array = array.T
        elif array.ndim == 3:
            array = Visualizer._slice_array(array, axis, n_slice)
            
        # Create the figure
        fig, ax = Visualizer._create_figure_2D(axis)
        
        # Plot the array
        options = imshow_options.copy()
        origin = options.pop('origin', Visualizer.DEFAULT_MPL_IMSHOW_ORIGIN)
        extent = options.pop('extent', Visualizer.DEFAULT_MPL_IMSHOW_EXTENT)
        im = ax.imshow(array, origin=origin, extent=extent, **options)
        
        # Plot the contours
        if contour_options is not None:
            cs = ax.contour(array, origin=origin, extent=extent,
                            **contour_options)
        
        # Set the colorbar
        if show_colorbar:
            cbar = fig.colorbar(im, ax=ax)
            if contour_options is not None:
                cbar.add_lines(cs)
        else:
            cbar = None
        
        return Figure(fig, ax, cbar)
        
    @staticmethod
    def _test_array_dimension(array: np.ndarray,
                              ) -> Union[ValueError, None]:
        """
        Test if the array is well a (2,) or (3,) dimensional array.
        """
        if array.ndim not in [2, 3]:
            msg = "Array dimension is not valid."
            raise ValueError(msg)
        else:
            return None
    
    @staticmethod
    def _slice_array(array: np.ndarray,
                     axis: str,
                     n_slice: int,
                     ) -> np.ndarray:
        """
        Slice a 3 dimensional ``array`` in the ``axis`` direction at the
        ``n_slice`` rank.
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
    def _create_figure_2D(axis: str,
                          ) -> tuple[matplotlib.figure.Figure, plt.Axes]:
        """
        Create a 2 dimensional matplotlib figure and set the label names
        according to the `axis`.
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
        return (fig, ax)
    
    @staticmethod
    def mpl_plot_array_3D(array: np.ndarray,
                          voxels_options: dict = {},
                          ) -> Figure:
        """
        Plot a 3D `array` with the matplotlib library.

        Parameters
        ----------
        array : np.ndarray
            3D array to plot.
        voxels_options : dict, optional
            Dictionary containing the arguments for the voxels() matplotlib
            function, by default {}.
            More details here:
            https://matplotlib.org/stable/api/_as_gen/mpl_toolkits.mplot3d.axes3d.Axes3D.voxels.html
            
        Returns
        -------
        Figure
            `NamedTuple` describing a matplotlib figure.
            - fig: TODO;
            - ax: TODO;
            - cbar: TODO;
            
        Examples
        --------
        >>> TODO
        """
        # Create the figure
        fig, ax = Visualizer._create_figure_3D()
        
        # Plot the array
        ax.voxels(array, **voxels_options)
        
        # Set aspect on 'equal'
        aspect = 'equal'
        ax.set_aspect(aspect)
        
        # Set the colorbar
        cbar = None
        
        return Figure(fig, ax, cbar)
    
    @staticmethod
    def _create_figure_3D() -> tuple[matplotlib.figure.Figure, plt.Axes]:
        """
        Create a 3 dimensional matplotlib figure and set the label names.
        """
        # Declare the figure
        fig = plt.figure()
        ax = plt.subplot(projection='3d')
        
        # Set axis labels
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        return (fig, ax)
    
    # **** PYVISTA **** #
    
    @staticmethod
    @requires_pyvista()
    def pv_plot_array(array: np.ndarray,
                      mask_values: list = [],
                      mask_subdomain: list = [],
                      show_grid: bool = True,
                      show_outline: bool = False,
                      threshold_options: dict = {},
                      cpos: Union[str, list[float, float, float]] = 'xy',
                      **kwargs
                      ) -> None:
        """
        Plot a 3D ``array`` with the pyvista library.

        Parameters
        ----------
        array : np.ndarray
            3D array to plot.
        mask_values : list, optional
            List of values to mask, by default ``[]``.
        mask_subdomain : list, optional
            List of pykasso subdomains to mask, by default ``[]``.
        show_grid : bool, optional
            If ``True``, plot the axis ans its values, by default ``True``.
        show_outline : bool, optional
            If ``True``, plot the full delimitation of the array regardless
            the values masked or filtered, by default ``False``.
        threshold_options : dict, optional
            Dictionary containing the arguments for the threshold() matplotlib
            function, by default ``{}``.
            More details here:
            https://docs.pyvista.org/version/stable/api/core/_autosummary/pyvista.DataSetFilters.threshold.html
        cpos : Union[str, list[float, float, float]]
            Initial camera position.
            More details here:
            https://docs.pyvista.org/version/stable/api/plotting/_autosummary/pyvista.Plotter.camera_position.html
        Returns
        -------
        None
        """

        # Control array dimension
        Visualizer._test_array_dimension(array)
        
        # Get array shape
        if array.ndim == 2:
            nx, ny = array.shape
            nz = 1
        elif array.ndim == 3:
            nx, ny, nz = array.shape

        # Create the mesh
        mesh = Visualizer._create_empty_mesh((nx, ny, nz))
        
        # Fill the mesh
        mesh.cell_data[Visualizer.PV_SCALAR_KEYWORD] = array.flatten(order="F")
        
        # Ghost values
        mesh = Visualizer._ghost_values(mesh, mask_values, mask_subdomain)
        
        # Create the plotter
        plotter = pv.Plotter()
        
        # Plot the grid
        if show_grid:
            plotter.show_grid()
            
        # Plot the outline
        if show_outline:
            plotter.add_mesh(mesh.outline(), color="k")
        
        # Filter then plot the data
        kwargs['scalars'] = Visualizer.PV_SCALAR_KEYWORD
        threshold_options.setdefault('scalars', Visualizer.PV_SCALAR_KEYWORD)
        mesh = mesh.threshold(**threshold_options)
        plotter.add_mesh(mesh, **kwargs)
        
        # Show the plotter
        plotter.show(cpos=cpos)

        return None
    
    @staticmethod
    def _create_empty_mesh(shape: tuple[int, int, int]
                           ) -> pyvista.core.pointset.UnstructuredGrid:
        """
        Create an empty pyvista mesh according to the `shape` tuple.

        Parameters
        ----------
        shape : tuple[int, int, int]
            Tuple of the shape of the mesh: (nx, ny, nz).

        Returns
        -------
        UnstructuredGrid
        """
        mesh = pv.UniformGrid()
        mesh.dimensions = np.array(shape) + 1
        mesh = mesh.cast_to_unstructured_grid()
        return mesh
    
    # TODO - should be rewritten
    @staticmethod
    def _ghost_values(mesh: pyvista.core.pointset.UnstructuredGrid,
                      ghost_values: list,
                      ghost_subdomains: list
                      ) -> pyvista.core.pointset.UnstructuredGrid:
        """
        DOC
        """
        ghosted_cells = []
        
        # Ghost cells according to data value
        if len(ghost_values) > 0:
            test = np.isin(mesh["data"], ghost_values)
            ghosted_cells.extend(np.argwhere(test))
        
        # TODO ??? improve the readability of the code
        # Ghost cells according to subdomain
        if len(ghost_subdomains) > 0:
            tests = []
            for subdomain in ghost_subdomains:
                tests.append(subdomain == 1)
            test_subdomains = np.logical_or.reduce(tests)
            data_domain = test_subdomains.flatten(order="F")
            mesh.cell_data['domain'] = data_domain
            ghosted_cells.extend(np.argwhere(np.isin(mesh["domain"], 1)))
            
        mesh = mesh.remove_cells(ghosted_cells)
        return mesh
    
    # ************************************************************* #
    # ****          MATPLOTLIB          *************************** #
    # ************************************************************* #
    
    def mpl_plot_2D(self,
                    n_sim: int = -1,
                    feature: str = 'karst',
                    n_iteration: int = -1,
                    axis: str = 'z',
                    n_slice: int = -1,
                    imshow_options: dict = {},
                    contour_options: dict = None,
                    show_colorbar: bool = True,
                    scatter_inlets_options: dict = None,
                    scatter_outlets_options: dict = None,
                    ) -> Figure:
        """
        Plot a slice of a 3D array of a feature from a specific simulation
        with the matplotlib library.

        Parameters
        ----------
        n_sim : int, optional
            Indice of the simulation to consider, by default -1.
        feature : str, optional
            _description_, by default 'karst'
        n_iteration : int, optional
            _description_, by default -1
        axis : str, optional
            Axis used to slice the array, by default 'z'.
        n_slice : int, optional
            Indice for slicing, by default -1.
        imshow_options : dict, optional
            Dictionary containing the arguments for the imshow() matplotlib
            function, by default {}.
            More details here:
            https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.imshow.html
        contour_options : dict, optional
            Dictionary containing the arguments for the contour() matplotlib
            function, by default None.
            More details here:
            https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.contour.html
        show_colorbar : bool, optional
            If `False`, remove the colorbar, by default `True`.
            More details here:
            https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.colorbar.html
        scatter_inlets_options : dict, optional
            _description_, by default None
        scatter_outlets_options : dict, optional
            _description_, by default None

        Returns
        -------
        Figure
            `NamedTuple` describing a matplotlib figure.
            - fig: TODO;
            - ax: TODO;
            - cbar: TODO;
            
        Examples
        --------
        >>> TODO
        """
        # Retrieve the data
        sim_data = self.project._get_simulation_data(n_sim)
        
        # Retrieve the data feature
        array = Visualizer._get_data_from_dict(simulation_data=sim_data,
                                               feature=feature,
                                               n_iteration=n_iteration)
        
        # Plot the feature data
        extent = self.project.grid.extent
        imshow_options.setdefault('extent', extent)
        f = Visualizer.mpl_plot_array_2D(array=array,
                                         axis=axis,
                                         n_slice=n_slice,
                                         imshow_options=imshow_options,
                                         contour_options=contour_options,
                                         show_colorbar=show_colorbar)
        
        # Plot the inlets
        if scatter_inlets_options is not None:
            inlets = sim_data['inlets']
            x = inlets.x
            y = inlets.y
            f.ax.scatter(x, y, **scatter_inlets_options)
        
        # Plot the outlets
        if scatter_outlets_options is not None:
            outlets = sim_data['outlets']
            x = outlets.x
            y = outlets.y
            f.ax.scatter(x, y, **scatter_outlets_options)
        
        return f
    
    def mpl_plot_3D(self,
                    n_sim: int = -1,
                    feature: str = 'karst',
                    n_iteration: int = -1,
                    voxels_options: dict = {},
                    ) -> Figure:
        """
        Plot a slice of a 3D array of a feature from a specific simulation
        with the matplotlib library.

        Parameters
        ----------
        n_sim : int, optional
            Indice of the simulation to consider, by default -1.
        feature : str, optional
            _description_, by default 'karst'
        n_iteration : int, optional
            _description_, by default -1
        voxels_options : dict, optional
            _description_, by default {}

        Returns
        -------
        Figure
            `NamedTuple` describing a matplotlib figure.
            - fig: TODO;
            - ax: TODO;
            - cbar: TODO;
            
        Examples
        --------
        >>> TODO
        """
        # Retrieve the data
        sim_data = self.project._get_simulation_data(n_sim)
        
        # Retrieve the data feature
        array = Visualizer._get_data_from_dict(sim_data, feature, n_iteration)
        
        # Plot the feature data
        f = self.mpl_plot_array_3D(array, voxels_options)
        
        return f
    
    # **** #
    
    def mpl_plot_karstic_network(self,
                                 n_sim: int = -1
                                 ) -> Figure:
        """
        Plot a simple view of the karstic conduit network's appearance.

        Parameters
        ----------
        n_sim : int, optional
            Indice of the simulation to consider, by default -1.

        Returns
        -------
        Figure
            `NamedTuple` describing a matplotlib figure.
            - fig: TODO;
            - ax: TODO;
            - cbar: TODO;
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
        
        # Set empty colorbar
        cbar = None
        
        return Figure(fig, ax, cbar)
    
    def _plot_graph(self,
                    ax: plt.Axes,
                    nodes: dict,
                    edges: list,
                    ) -> plt.Axes:
        """
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

    # ************************************************************* #
    # ****          PYVISTA             *************************** #
    # ************************************************************* #
    
    def _set_pyvista_grid(self) -> None:
        """
        Create an empty mesh representing the grid of the simulations for
        future pyvista plot.
        """
        # Get grid parameters
        grid = self.project.grid.get_grid_parameters()
        x0, y0, z0 = grid['x0'], grid['y0'], grid['z0']
        nx, ny, nz = grid['nx'], grid['ny'], grid['nz']
        dx, dy, dz = grid['dx'], grid['dy'], grid['dz']
        pv_grid = pv.UniformGrid()
        
        # Construct grid
        pv_grid.origin = (x0 - dx / 2, y0 - dy / 2, z0 - dz / 2)
        pv_grid.dimensions = np.array((nx, ny, nz)) + 1
        pv_grid.spacing = (dx, dy, dz)
        self.pv_grid = pv_grid.cast_to_unstructured_grid()
        return None
    
    @requires_pyvista()
    def pv_show(self,
                simulations: list,
                features: list,
                settings: list = [],
                cpos: str = 'xz',
                savefig: bool = False,
                ) -> None:
        """
        TODO
        """
        
        #############################
        ### Parameters validation ###
        #############################
        
        # If 'simulation' is an integer, transform it in a list
        if isinstance(simulations, int):
            simulations = [simulations]
        
        # Test if 'simulations' is a list
        if not isinstance(simulations, list):
            msg = "Parameter 'simulations' must be of type: list."
            raise TypeError(msg)
        
        # If 'features' is a string, transform it in a list
        if isinstance(features, str):
            features = [features]
        
        # Test if 'features' is a list
        if not isinstance(features, list):
            msg = "Parameter 'features' must be of type: list."
            raise TypeError(msg)
        
        # If 'settings' has not been set
        if settings == []:
            settings = [{}] * len(simulations)
            
        # If 'settings' is a dictionary, transform it in a list
        if isinstance(settings, dict):
            settings = [settings] * len(simulations)
        
        # TODO - to remove
        # # Check type of 'settings' parameter
        # if not isinstance(settings, (list)):
        #     settings = [settings] * len(simulations)
        
        ######################
        ### Initialization ###
        ######################
        
        # Insert default settings values
        pv_settings = []
        for settings_dict in settings:
            settings_dict = Visualizer._set_default_pv_settings(settings_dict)
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
        if savefig:
            off_screen = True
        else:
            off_screen = False
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
                
                text = 'Simulation {} - {}'.format(n_sim, feature)
                settings.setdefault('text_options', {})
                settings['text_options'].setdefault('text', text)
                
                data = self.project._get_simulation_data(n_sim)
                plotter.subplot(row, column)
                plotter = self._fill_plotter(plotter,
                                             data,
                                             feature,
                                             settings)
        else:
            for (row, (n_sim, settings)) in enumerate(zip(simulations,
                                                          pv_settings)
                                                      ):
                for (column, feature) in enumerate(features):
                    
                    text = 'Simulation {} - {}'.format(n_sim, feature)
                    settings.setdefault('text_options', {})
                    settings['text_options'].setdefault('text', text)
                    
                    data = self.project._get_simulation_data(n_sim)
                    plotter.subplot(row, column)
                    plotter = self._fill_plotter(plotter,
                                                 data,
                                                 feature,
                                                 settings)
                    
        plotter.link_views()
        
        if savefig:
            outputs_dir = self.project.core['paths']['outputs_dir']
            filename = outputs_dir + 'pv_plot'
        #     path = self.project_directory + '/outputs/' + savefig
            # plotter.show(cpos=cpos, screenshot=filename)
            plotter.show(cpos=cpos, auto_close=False)
            plotter.show(cpos=cpos, screenshot=filename)
        else:
            plotter.show(cpos=cpos)
            
        return None
    
    @staticmethod
    def _set_default_pv_settings(settings: dict) -> dict:
        """
        To a dictionary, set the default settings for pyvista.
        """
        for (parameter, value) in Visualizer.DEFAULT_PV_SETTINGS.items():
            settings.setdefault(parameter, value)
        return settings
    
    def _fill_plotter(self,
                      plotter: pv.Plotter,
                      simulation_data: dict,
                      feature: str,
                      settings: dict
                      ) -> pv.Plotter:
        """_summary_

        Parameters
        ----------
        plotter : pv.Plotter
            _description_
        simulation_data : dict
            _description_
        feature : str
            _description_
        settings : dict
            _description_

        Returns
        -------
        pv.Plotter
            _description_

        Raises
        ------
        e
            _description_
        """
        text_options = settings.pop('text_options')
        
        # Test if asked feature is valid
        if feature not in AUTHORIZED_FEATURES:
            msg = "Asked feature '{}' is invalid. Authorized features: {}."
            msg = msg.format(feature, AUTHORIZED_FEATURES)
            warnings.warn(msg)
            text_options['text'] = '{} : Invalid feature'.format(feature)
        
        # Otherwise retrieve data
        else:
            try:
                actors = self._get_actors(simulation_data=simulation_data,
                                          feature=feature,
                                          **settings)
            
            # If  ?  # TODO
            except IndexError:
                msg = "Asked feature '{}' has no data.".format(feature)
                warnings.warn(msg)
                text_options['text'] = '{} : No data'.format(feature)
            
            # If ?  # TODO
            except Exception as e:
                raise e
            
            # Generate and plot actors
            else:
                for actor in actors:
                    plotter.add_actor(actor, reset_camera=True)
           
        # Print figure caption
        plotter.add_text(**text_options)
            
        return plotter
    
    def _get_actors(self,
                    simulation_data: dict,
                    feature: str,
                    n_iteration: int = -1,
                    ghost_values: list = [],
                    ghost_subdomains: list = [],
                    show_grid: bool = True,
                    show_outline: bool = False,
                    data_options: dict = {},
                    surfaces_options: dict = {},
                    inlets_options: dict = None,
                    outlets_options: dict = None,
                    show_slice: bool = False,
                    fractures_options: dict = None,
                    show_colorbar: bool = True,
                    ) -> list[pv.Actor]:
        """
        DOC
        """
        # Get data from feature
        feature_data = Visualizer._get_data_from_dict(simulation_data,
                                                      feature,
                                                      n_iteration)
        
        # Create the mesh
        mesh = self.pv_grid.copy()

        # Fill the mesh
        mesh.cell_data['data'] = feature_data.flatten(order="F")
            
        # Ghost the data
        test_ghost_values = (len(ghost_values) > 0)
        test_ghost_subdomains = (len(ghost_subdomains) > 0)
        if test_ghost_values or test_ghost_subdomains:
            
            domain = simulation_data['domain']
            # If necessery, retrieves the subdomains
            if test_ghost_subdomains > 0:
                subdomains = []
                for subdomain in ghost_subdomains:
                    if subdomain[-2:] == '_r':
                        subdomain = subdomain[:-2]
                        data_domain = domain.get_subdomain(subdomain)
                        data_domain = np.logical_not(data_domain).astype(int)
                    else:
                        data_domain = domain.get_subdomain(subdomain)
                    subdomains.append(data_domain)
            else:
                subdomains = []
            
            mesh = Visualizer._ghost_values(mesh, ghost_values, subdomains)
        
        ### Create the plot
        plotter = pv.Plotter()
        actors = []
        
        # Plot the grid of the domain
        if show_grid:
            _ = plotter.show_grid()
            actors.append(_)
        
        # Plot the outline of the domain
        if show_outline:
            _ = plotter.add_mesh(self.pv_grid.outline(), color="k")
            actors.append(_)
            
        # Plot the inlets
        if inlets_options is not None:
            inlets_options.setdefault('render_points_as_spheres', False)
            inlets_options.setdefault('point_size', 20)
            inlets_options.setdefault('color', 'r')

            inlets = Visualizer._df_to_cloud(simulation_data['inlets'])
            _ = plotter.add_points(inlets, **inlets_options)
            actors.append(_)
        
        # Plot the outlets
        if outlets_options is not None:
            outlets_options.setdefault('render_points_as_spheres', False)
            outlets_options.setdefault('point_size', 20)
            outlets_options.setdefault('color', 'b')
            
            outlets = Visualizer._df_to_cloud(simulation_data['outlets'])
            _ = plotter.add_points(outlets, **outlets_options)
            actors.append(_)
            
        # Plot the surfaces
        if surfaces_options is not None:
            meshgrids = self.project.grid.get_meshgrids()
            
            for surface_name, options in surfaces_options.items():
                options.setdefault('smooth_shading', True)
                options.setdefault('opacity', 0.66)
                # options.setdefault('scalars', 'data')
                surface = Visualizer._get_surface(meshgrids,
                                                  simulation_data,
                                                  surface_name)
                if surface is not None:
                    _ = plotter.add_mesh(surface, **options)
                    actors.append(_)
            
        # Plot the fractures
        if fractures_options is not None:
            # Set default parameters
            f = simulation_data['fractures'].fractures
            fractures_options.setdefault('fractures', f)
            fid = f['family_id'].unique().tolist()
            fractures_options.setdefault('family_id', fid)
            fractures_options.setdefault('sort', ['radius'])
            fractures_options.setdefault('n', 250)
            
            if len(fractures_options['fractures']) > 0:
                # Sort fractures
                fractures = Visualizer._sort_fractures(**fractures_options)
                
                # Get polygons
                polygons = Visualizer._fractures_to_polygons(fractures)
                
                # Plot polygons
                options_polygons = {}
                options_polygons.setdefault('smooth_shading', True)
                options_polygons.setdefault('opacity', 0.66)
                for polygon in polygons:
                    _ = plotter.add_mesh(polygon, **options_polygons)
                    actors.append(_)
        
        # Plot the data
        data_options.setdefault('scalar_bar_args', {'title': 'Vol1'})
        data_options.setdefault('scalars', 'data')
        if show_slice:
            _ = plotter.add_mesh_slice_orthogonal(mesh=mesh,
                                                  generate_triangles=False,
                                                  widget_color=None,
                                                  tubing=False,
                                                  interaction_event=45,
                                                  **data_options)
            actors.extend(_)
            
        else:
            _ = plotter.add_mesh(mesh.copy(), **data_options)
            actors.append(_)
        
        # Plot the scalar bar
        if show_colorbar:
            _ = plotter.add_scalar_bar()
            actors.append(_)

        return actors
    
    ### TODO - to rewrite
    @staticmethod
    def _get_data_from_dict(simulation_data: dict,
                            feature: str,
                            n_iteration: int,
                            ) -> np.ndarray:
        """
        Return the data from the ``feature`` key contained in the
        ``simulation_data`` dictionary.

        Parameters
        ----------
        simulation_data : dict
            Dictionary containing data from the simulation.
        feature : str
            Name of the geologic feature to retrieve.
        n_iteration : int
            Iteration to consider to for retrieving data.

        Returns
        -------
        np.ndarray
            Numpy array of the required geologic feature.

        Raises
        ------
        ValueError
            If the geologic feature name is invalid.
        AttributeError
            If the geologic feature required is empty.
        """
        # Select the adequate way to retrieve data
        if feature in GEOLOGICAL_FEATURES:
            featured_data = simulation_data[feature].data_volume
        elif feature in DOMAIN_FEATURES:
            featured_data = getattr(simulation_data['domain'], feature)
            featured_data = featured_data.data_volume
        elif feature in ANISOTROPIC_FEATURES:
            featured_data = simulation_data['maps'][feature][n_iteration]
            # else:  # TODO
                # msg = 'feature keyword error'
                # raise ValueError(msg)
        # except AttributeError:
            
        # except IndexError:
            # msg = "Warning. A simulation computed with "
            
        return featured_data
    
    @staticmethod
    @requires_pyvista()
    def _df_to_cloud(points: pd.DataFrame) -> pv.PolyData:
        """
        Transform a pandas ``DataFrame`` storing points into a pyvista
        ``PolyData``
        
        Parameters
        ----------
        points : pd.DataFrame
            ``DataFrame`` of points to plot.

        Returns
        -------
        pv.PolyData
            pyvista object representing the points to plot.
        """
        labels = points.index.values.tolist()
        points = points[['x', 'y', 'z']].values
        points = points.astype('float32')
        cloud = pv.wrap(points)
        cloud['labels'] = labels
        return cloud
    
    @staticmethod
    @requires_pyvista()
    def _get_surface(meshgrids: tuple[np.ndarray, np.ndarray, np.ndarray],
                     simulation_data: dict,
                     surface_name: str
                     ) -> Union[pv.StructuredGrid, None]:
        """
        Plot a surface feature. If the surface name is invalid, or if the
        surface feature is absent from the current simulation, print a warning
        and return ``None``.

        Parameters
        ----------
        meshgrids : tuple[np.ndarray, np.ndarray, np.ndarray]
            Tuple of numpy arrays representing the meshgrids of the simulation.
        simulation_data : dict
            Data of the simulation.
        surface_name : str
            Name of the surface to plot. Valid names : 'topography',
            'water_level', 'bedrock'.

        Returns
        -------
        Union[pv.StructuredGrid, None]
            The pyvista surface object. ``None`` if the surface does not exist.
        """
        X, Y, Z = meshgrids
        x = X[:, :, 0]
        y = Y[:, :, 0]
        # Control validity of the surface name
        if surface_name not in SURFACE_FEATURES:
            msg = ("'{}' is not a valid name surface. Valid names: {}"
                   .format(surface_name, SURFACE_FEATURES))
            warnings.warn(msg)
            return None
        try:
            surface_object = getattr(simulation_data['domain'], surface_name)
            z = surface_object.data_surface
            surface = pv.StructuredGrid(x, y, z)
            surface['data'] = z.flatten(order="F")
            return surface
        except AttributeError:
            msg = "Simulation has no '{}' surface.".format(surface_name)
            warnings.warn(msg)
            return None
    
    @staticmethod
    def _sort_fractures(fractures: pd.DataFrame,
                        family_id: list[int],
                        sort: list[str],
                        max_number: int,
                        ) -> list:
        """
        Filter the pandas ``DataFrame`` storing the fractures.

        Parameters
        ----------
        fractures : pd.DataFrame
            ``DataFrame`` storing the fractures.
        family_id : list[int]
            List of family ids to plot.
        sort : list[str]
            List of ``DataFrame`` attributes considered to sort the values.
        max_number : int
            Maximum number of fractures to plot.

        Returns
        -------
        list
            List of fractures passing the filter.
        """
        # Keep fractures according to the 'family_id' parameter
        f = fractures[fractures['family_id'].isin(family_id)]
        
        # Sort and keep fractures to the 'sort' and 'n' parameters
        f = f.sort_values(sort, ascending=False).iloc[:max_number]
        
        # Retrieve data from fractures and convert it to python list
        f = f[['x', 'y', 'z', 'radius', 'normal']].to_numpy().tolist()
        return f
        
    @staticmethod
    @requires_pyvista()
    def _fractures_to_polygons(fractures: list) -> list[pv.PolyData]:
        """
        Transform a list of fractures in a list of pyvista polygons.

        Parameters
        ----------
        fractures : list
            List of fractures to plot.

        Returns
        -------
        list[pv.PolyData]
            List of pyvista polygons.
        """
        # Generate the polygons
        polygons = []
        for (x, y, z, r, n) in fractures:
            poly = pv.Polygon(center=(x, y, z), radius=r, normal=n, n_sides=4)
            polygons.append(poly)
    
        return polygons




class SetVisibilityCallback:
    """Helper callback to keep a reference to the actor being modified."""

    def __init__(self, actor):
        self.actor = actor

    def __call__(self, state):
        self.actor.SetVisibility(state)




# ????
#     @requires_pyvista()
#     def pv_plot(self, n_sim, feature, settings={}, show=True,
#                 return_plotter=False):
#         """
#         DOC
#         """
#         ######################
#         ### Initialisation ###
#         ######################
        
#         # Create plotter
#         plotter = pv.Plotter()
        
#         # Check 'settings' parameter
#         settings = self._check_pv_settings(settings)
#         if 'text' not in settings:
#             text = 'Simulation {} - {}'.format(n_sim, feature)
#             settings['text'] = text
            
#         # Get simulation data
#         simulation_data = self.project._get_simulation_data(n_sim)
        
#         ###############
#         ### Ploting ###
#         ###############
        
#         # Plot actors
#         plotter = self._fill_plotter(plotter, simulation_data,
#                                      feature, settings)
        
#         # Plot figure
#         if show:
#             plotter.show()
        
#         if return_plotter:
#             return plotter
#         else:
#             return None
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
#     #     def create_gif(self, simulation: int, feature: str, location: str,
# #                    zoom: float = 1, ghosts: list = [], n_points: int = 24,
# #                    fps: int = 10, window_size=[1024, 768],
# #                    colormap: str = 'viridis',
# #                    background_color: str = 'white',
# #                    background_color_top: str = None) -> None:
# #         """
# #         TODO
        
# #         colormap: str
# #             https://matplotlib.org/stable/tutorials/colors/colormaps.html
# #         """
        
# #         ### Method based on those examples:
# #         # https://docs.pyvista.org/examples/02-plot/orbit.html#orbiting
        
# #         ### Gets the simulation data
# #         simulation_data = self._get_simulation_data(simulation)
        
# #         ### Gets the mesh
# #         mesh = self._get_data_from_feature(simulation_data, feature)
# #         # mesh_ = mesh.copy()
# #         if len(ghosts) > 0:
# #             mesh = self._ghost_data(mesh, ghosts)
        
# #         ### Constructs the plotter
# #         plotter = pv.Plotter(off_screen=True, window_size=window_size)
# #         plotter.store_image = True
        
# #         # Plots the data
# #         kwargs = {
# #             'cmap': colormap,
# #             'scalar_bar_args': {'title': 'Vol1'},
# #             'scalars': 'data',
# #             'lighting': True,
# #         }
# #         plotter.add_mesh(mesh, **kwargs)
# #         plotter.remove_scalar_bar()
        
# #         # Sets background color
# #         if background_color_top is None:
# #             background_color_top = background_color
# #         plotter.set_background(background_color)
        
# #         ### Sets the initial camera position
# #         plotter.camera.zoom(zoom)
# #         # plotter.camera.roll = 0
# #         plotter.camera.elevation = 0
# #         plotter.camera.azimuth = 0
        
# #         ### Generates the GIF
        
# #         # Open a gif
# #         plotter.open_gif(location, fps=fps)
        
# #         # Creates the camera positions
# #         azimuth_step = 360 / n_points
        
# #         # Loops
# #         for i in range(n_points):
            
# #             # Updates camera position
# #             # plotter.camera.roll
# #             # plotter.camera.elevation
# #             plotter.camera.azimuth += azimuth_step
            
# #             # Updates the background
            
# #             # Writes the frame in the gif
# #             plotter.write_frame()
        
# #         # Closes and finalizes movie
# #         plotter.close()
# #         return None




# ??????????????????????????????????
# def show_simulation(self):
#         """
#         TODO
#         """
#         def _show_data(environment, feature, settings):
        
#             # Gets the domain
#             if hasattr(environment, 'domain') and (getattr(environment, 'domain') is not None):
#                 grid = _get_data_from_attribute(grid, getattr(environment, 'domain'), 'data_volume', 'domain')

#             # Ghost the data
#             if 'ghost' in settings:
#                 if settings['domain'] == True:
#                     ghosts = np.argwhere(np.logical_or(np.isin(grid["data"], settings['ghost']), (grid["domain"] == 0)))
#                 else:
#                     ghosts = np.argwhere(np.isin(grid["data"], settings['ghost']))
#             else:
#                 if settings['domain'] == True:
#                     ghosts = np.argwhere(grid["domain"] == 0)
#                 else:
#                     ghosts = None
            
#             if ghosts is not None:
#                 grid = grid.remove_cells(ghosts)

#         return None