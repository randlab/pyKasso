"""
Module defining a tool for pyKasso results visulization.
"""

### Internal dependencies
import datetime
import copy
import warnings
import platform

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
        'text_options': {},
        'n_iteration': -1,
        'mask_values': [],
        'mask_subdomains': [],
        'show_grid': True,
        'show_outline': False,
        'threshold_options': {},
        'inlets_options': None,
        'outlets_options': None,
        'surfaces_options': {},
        'fractures_options': None,
        'data_options': {},
        'show_slice': False,
        'show_scalar_bar': True,
        'scalar_bar_args': {},
        'discrete_scalar_bar': False,
        'background_color': 'white',
        'cmap': 'viridis',
    }
    
    DEFAULT_MPL_IMSHOW_ORIGIN = 'lower'
    DEFAULT_MPL_IMSHOW_EXTENT = None
    PV_SCALAR_KEYWORD = 'data'
    
    def __init__(self,
                 project: Project,
                 notebook: bool = False,
                 jupyter_backend: str = 'static',
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
        jupyter_backend : str, optional
            TODO
        """
        self.project = project
        self._notebook = notebook
        self._jupyter_backend = jupyter_backend
        
        # Set notebook value
        self._set_notebook_value(notebook)
        
        # Set the pyvista grid
        if _has_pyvista:
            self.pv_grid = self._get_pyvista_grid_from_project()
            
        # Set the jupyter backend
        if _has_pyvista:
            self._set_pyvista_jupyter_backend(jupyter_backend)
        
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
        if not isinstance(boolean, bool):
            msg = "Input must be a boolean."
            raise TypeError(msg)
        
        self._notebook = boolean
        self._set_notebook_value(boolean)
    
    @staticmethod
    def _set_notebook_value(boolean: bool) -> None:
        """
        Render plots within python notebooks if `boolean` is `True`.
        """
        
        # *** MATPLOTLIB *** #
        
        # First test system
        if platform.system() == 'Linux':
            get_ipython().run_line_magic('matplotlib', 'inline')
        else:
            if boolean:
                get_ipython().run_line_magic('matplotlib', 'inline')
            else:
                get_ipython().run_line_magic('matplotlib', 'qt')
            
        # *** PYVISTA *** #
        if _has_pyvista:
            pv.global_theme.notebook = boolean
        
        return None
    
    @property
    def jupyter_backend(self) -> bool:
        return self._jupyter_backend

    @jupyter_backend.setter
    def jupyter_backend(self,
                        mode: str,
                        ) -> None:
        self._set_pyvista_jupyter_backend(mode)
        self._jupyter_backend = mode
        
    @staticmethod
    def _set_pyvista_jupyter_backend(mode: str) -> None:
        pv.set_jupyter_backend(mode)
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
                          show_colorbar: bool = False,
                          show_grid: bool = False,
                          show_axis: bool = True,
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
        show_grid : bool, optional
            If `False`, remove the grid, by default `False`.
        show_axis : bool, optional
            If `False`, remove the axis, by default `True`.
            
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
        ax.grid(show_grid)
        
        # Plot the array
        options = imshow_options.copy()
        origin = options.pop('origin', Visualizer.DEFAULT_MPL_IMSHOW_ORIGIN)
        extent = options.pop('extent', Visualizer.DEFAULT_MPL_IMSHOW_EXTENT)
        im = ax.imshow(array, origin=origin, extent=extent, **options)
        
        # Plot the contours
        if contour_options is not None:
            cs = ax.contour(array, origin=origin, extent=extent,
                            **contour_options)
            
        # Set the axis
        if show_axis is False:
            # plt.axis('off')
            ax.set_xlabel(None)
            ax.set_ylabel(None)
            ax.set_xticklabels([])
            ax.set_yticklabels([])
                
            ticks_ = [ax.xaxis.get_major_ticks(), ax.yaxis.get_major_ticks()]
            for ticks in ticks_:
                for tick in ticks:
                    tick.tick1line.set_visible(False)
                    tick.tick2line.set_visible(False)
                    tick.label1.set_visible(False)
                    tick.label2.set_visible(False)
        
        # Set the colorbar
        if show_colorbar:
            cbar = fig.colorbar(im, ax=ax)
            if contour_options is not None:
                cbar.add_lines(cs)
        else:
            cbar = None
        
        return Figure(fig, ax, cbar)
        
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
                      settings: dict = {},
                      grid_origin: tuple = (0, 0, 0),
                      grid_spacing: tuple = (1, 1, 1),
                      cpos: Union[str, list[float, float, float]] = 'xy',
                      zoom: float = None,
                      window_size: tuple = None,
                      savefig: bool = False,
                      filename: str = None,
                      return_plotter: bool = False,
                      ) -> None:
        """ DOC """
        
        # Control array dimension
        Visualizer._test_array_dimension(array)
        
        # Construct the pyvista grid and the 'simulation_data' dict
        simulation_data = {}
        mesh = Visualizer._get_pyvista_grid_from_array(array,
                                                       grid_origin,
                                                       grid_spacing)
        simulation_data['mesh'] = mesh
        simulation_data['data'] = array
        
        # Construct the plotter
        plotter = Visualizer._construct_plotter(shape=(1, 1),
                                                border=True)
        
        # Fill the plotter
        plotter = Visualizer._fill_plotter(plotter,
                                           simulation_data,
                                           settings)
        plotter.link_views()
        
        # Show the plotter
        result = Visualizer._show_plotter(plotter, cpos, window_size, zoom,
                                          savefig, filename, return_plotter)
        return result
    
    @staticmethod
    def _test_array_dimension(array: np.ndarray,
                              ) -> Union[ValueError, None]:
        """
        Test if the array is well a (2,) or (3,) dimensional array.
        """
        if array.ndim not in [2, 3]:
            msg = "Array dimensions are not valid."
            raise ValueError(msg)
        else:
            return None
        
    @staticmethod
    @requires_pyvista()
    def _get_pyvista_grid_from_array(array: np.array,
                                     origin: tuple,
                                     spacing: tuple,
                                     ):  # pv.ImageData
        """ DOC """
        # Get array shape
        if array.ndim == 2:
            nx, ny = array.shape
            nz = 1
        elif array.ndim == 3:
            nx, ny, nz = array.shape
        x0, y0, z0 = origin
        dx, dy, dz = spacing
        grid = [x0, y0, z0, nx, ny, nz, dx, dy, dz]
        pv_grid = Visualizer._construct_pyvista_grid(*grid)
        return pv_grid
    
    @staticmethod
    @requires_pyvista()
    def _construct_pyvista_grid(x0, y0, z0, nx, ny, nz, dx, dy, dz):
        """ DOC """
        pv_grid = pv.ImageData()
        pv_grid.origin = (x0 - dx / 2, y0 - dy / 2, z0 - dz / 2)
        pv_grid.dimensions = np.array((nx, ny, nz)) + 1
        pv_grid.spacing = (dx, dy, dz)
        pv_grid = pv_grid.cast_to_unstructured_grid()
        return pv_grid
                      
    @staticmethod
    @requires_pyvista()
    def _construct_plotter(shape: tuple,
                           border: bool = True,
                           ):  # -> pv.Plotter
        """ DOC """
        plotter = pv.Plotter(shape=shape, border=border)
        return plotter
    
    @staticmethod
    @requires_pyvista()
    def _fill_plotter(plotter,  # pv.Plotter
                      simulation_data: dict,
                      settings: dict,
                      ):  # -> pv.Plotter
        """ DOC """
        
        # Retrieve text options
        text_options = settings.pop('text_options', None)
        
        # Retrieve data
        actors = Visualizer._get_actors(simulation_data,
                                        **settings)
        
        # Generate and plot actors
        for actor in actors:
            plotter.add_actor(actor, reset_camera=True)
           
        # Print figure caption
        if text_options is not None:
            plotter.add_text(**text_options)
            
        return plotter
    
    @staticmethod
    def _get_actors(simulation_data: dict,
                    mask_values: list = [],
                    mask_subdomains: list = [],
                    show_grid: bool = True,
                    show_outline: bool = False,
                    threshold_options: dict = {},
                    inlets_options: dict = None,
                    outlets_options: dict = None,
                    surfaces_options: dict = {},
                    fractures_options: dict = None,
                    data_options: dict = {},
                    show_slice: bool = False,
                    show_scalar_bar: bool = True,
                    scalar_bar_args: dict = {},
                    discrete_scalar_bar: bool = False,
                    background_color: str = 'white',
                    cmap: str = 'viridis',  # TODO Union
                    ):  # -> list[pv.Actor]
        """
        DOC
        """
    #     """
    #     Plot a 3D ``array`` with the pyvista library.

    #     Parameters
    #     ----------
    #     array : np.ndarray
    #         3D array to plot.
    #     mask_values : list, optional
    #         List of values to mask, by default ``[]``.
    #     mask_subdomain : list, optional
    #         List of pykasso subdomains to mask, by default ``[]``.
    #     show_grid : bool, optional
    #         If ``True``, plot the axis ans its values, by default ``True``.
    #     show_outline : bool, optional
    #         If ``True``, plot the full delimitation of the array regardless
    #         the values masked or filtered, by default ``False``.
    #     threshold_options : dict, optional
    #         Dictionary containing the arguments for the threshold() matplotlib
    #         function, by default ``{}``.
    #         More details here:
    #         https://docs.pyvista.org/version/stable/api/core/_autosummary/pyvista.DataSetFilters.threshold.html
    #     text_options: dict, optional
    #         __doc__ , by default ``None``.
    #         More details here:
    #         https://docs.pyvista.org/version/stable/api/plotting/_autosummary/pyvista.Plotter.add_text.html#pyvista.Plotter.add_text.html
    #     show_scalar_bar: bool, optional
    #         __doc__ , by default ``True``.
    #     scalar_bar_args: dict, optional
    #         __doc__ , by default ``{}``.
    #     discrete_scalar_bar: bool, optional
    #         __doc__ , by default ``False``.
    #     cmap: str, optional
    #         __doc__ , by default ``viridis``.
    #     cpos : Union[str, list[float, float, float]], optional
    #         Initial camera position.
    #         More details here:
    #         https://docs.pyvista.org/version/stable/api/plotting/_autosummary/pyvista.Plotter.camera_position.html
    #     background_color : str, optional
    #         __doc__, by default ``'white'``.
    #     savefig : bool, optional
    #         __doc__, by default ``False``.
    #     filename : str, optional
    #         __doc__, by default ``None``.
            
    #     """
    
        # Retrieve the mesh
        mesh = simulation_data['mesh']
        mesh0 = mesh.copy()
        
        # Retrieve the data and fill the mesh
        data = simulation_data[Visualizer.PV_SCALAR_KEYWORD]
        mesh.cell_data[Visualizer.PV_SCALAR_KEYWORD] = data.flatten(order="F")
            
        # Ghost the data
        test_ghost_values = (len(mask_values) > 0)
        test_ghost_subdomains = (len(mask_subdomains) > 0)
        if test_ghost_values or test_ghost_subdomains:
            
            if 'domain' in simulation_data:
                domain = simulation_data['domain']
            else:
                domain = data.copy()
            # If necessery, retrieves the subdomains
            if test_ghost_subdomains > 0:
                subdomains = []
                for subdomain in mask_subdomains:
                    if subdomain[-2:] == '_r':
                        subdomain = subdomain[:-2]
                        data_domain = domain.get_subdomain(subdomain)
                        data_domain = np.logical_not(data_domain).astype(int)
                    else:
                        data_domain = domain.get_subdomain(subdomain)
                    subdomains.append(data_domain)
            else:
                subdomains = []
            
            mesh = Visualizer._mask_values(mesh, mask_values, subdomains)
        
        ### Create the plot and the list of actors
        plotter = pv.Plotter()
        plotter.background_color = background_color
        actors = []
        
        # Plot the grid of the domain
        if show_grid:
            _ = plotter.show_grid()
            actors.append(_)
        
        # Plot the outline of the domain
        if show_outline:
            _ = plotter.add_mesh(mesh0.outline(), color="k")
            actors.append(_)
           
        # Filter the data
        threshold_options.setdefault('scalars', Visualizer.PV_SCALAR_KEYWORD)
        mesh = mesh.threshold(**threshold_options)
        
        if discrete_scalar_bar:
            n = len(np.unique(mesh.cell_data[Visualizer.PV_SCALAR_KEYWORD]))
            cmap = plt.cm.get_cmap(cmap, n)
        
        # Plot the inlets
        if (inlets_options is not None) and ('inlets' in simulation_data):
            inlets_options.setdefault('render_points_as_spheres', False)
            inlets_options.setdefault('point_size', 20)
            inlets_options.setdefault('color', 'r')

            inlets = Visualizer._df_to_cloud(simulation_data['inlets'])
            _ = plotter.add_points(inlets, **inlets_options)
            actors.append(_)
        
        # Plot the outlets
        if (outlets_options is not None) and ('outlets' in simulation_data):
            outlets_options.setdefault('render_points_as_spheres', False)
            outlets_options.setdefault('point_size', 20)
            outlets_options.setdefault('color', 'b')
            
            outlets = Visualizer._df_to_cloud(simulation_data['outlets'])
            _ = plotter.add_points(outlets, **outlets_options)
            actors.append(_)
            
        # Plot the surfaces
        if surfaces_options is not None:
            for surface_name, options in surfaces_options.items():
                options.setdefault('smooth_shading', True)
                options.setdefault('opacity', 0.66)
                # options.setdefault('scalars', 'data')
                surface = Visualizer._get_surface(simulation_data,  # TODO
                                                  surface_name)
                if surface is not None:
                    _ = plotter.add_mesh(surface, **options)
                    actors.append(_)
            
        # Plot the fractures
        if (fractures_options is not None) and ('fractures' in simulation_data):
            # Set default parameters
            f = simulation_data['fractures'].fractures
            fractures_options.setdefault('fractures', f)
            fid = f['family_id'].unique().tolist()
            fractures_options.setdefault('family_id', fid)
            fractures_options.setdefault('sort', ['radius'])
            fractures_options.setdefault('max_number', 250)
            
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
       
        ### Plot the data
        scalar_keyword = Visualizer.PV_SCALAR_KEYWORD
        data_options.setdefault('scalars', scalar_keyword)
        
        # Set the scalar bar
        if discrete_scalar_bar:
            n = len(np.unique(mesh.cell_data[scalar_keyword]))
            cmap = plt.cm.get_cmap(cmap, n)
        
        # Plot the sliced data
        if show_slice:
            _ = plotter.add_mesh_slice_orthogonal(
                mesh=mesh,
                generate_triangles=False,
                widget_color=None,
                tubing=False,
                interaction_event=45,
                cmap=cmap,
                **data_options
            )
            actors.extend(_)
            
        else:
            _ = plotter.add_mesh(mesh.copy(),
                                 show_scalar_bar=show_scalar_bar,
                                 scalar_bar_args=scalar_bar_args,
                                 cmap=cmap,
                                 **data_options)
            actors.append(_)
           
        # Set the colorbar
        if show_scalar_bar:
            # scalar_bar_args.setdefault('title', 'Vol1')
            _ = plotter.add_scalar_bar(**scalar_bar_args)
            actors.append(_)

        return actors
    
    @staticmethod
    @requires_pyvista()
    def _show_plotter(plotter,
                      cpos: Union[str, list[float, float, float]] = 'xy',
                      window_size: tuple = None,
                      zoom: float = None,
                      savefig: bool = False,
                      filename: str = None,
                      return_plotter: bool = False,
                      ):
        """ DOC """
        
        # Set camera position
        plotter.camera_position = cpos
        if zoom is not None:
            plotter.camera.zoom(zoom)
            
        if savefig:
            # outputs_dir = self.project.core['paths']['outputs_dir']
            date = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
            
            if filename is None:
                filename = 'pv_plot_' + date
                # filename = outputs_dir + 'pv_plot_' + date
        
            plotter.show(
                window_size=window_size,
                screenshot=filename,
            )
        else:
            if return_plotter:
                return plotter
            else:
                plotter.show(window_size=window_size)
        return None

    # TODO - should be rewritten
    @staticmethod
    def _mask_values(mesh,  # pyvista.core.pointset.UnstructuredGrid
                     mask_values: list,
                     mask_subdomains: list,
                     ):  # pyvista.core.pointset.UnstructuredGrid
        """
        DOC
        """
        ghosted_cells = []
        
        # Ghost cells according to data value
        if len(mask_values) > 0:
            test = np.isin(mesh["data"], mask_values)
            ghosted_cells.extend(np.argwhere(test))
        
        # TODO ??? improve the readability of the code
        # Ghost cells according to subdomain
        if len(mask_subdomains) > 0:
            tests = []
            for subdomain in mask_subdomains:
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
                    show_grid: bool = False,
                    show_axis: bool = True,
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
        show_grid : bool, optional
            If `True`, print the grid, by default `False`.
        show_axis : bool, optional
            If `False`, remove the axis, by default `True`.
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
                                         show_colorbar=show_colorbar,
                                         show_grid=show_grid,
                                         show_axis=show_axis)
        
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
    
    @requires_pyvista()
    def _get_pyvista_grid_from_project(self):
        """ DOC """
        grid = self.project.grid.get_grid_parameters()
        x0, y0, z0 = grid['x0'], grid['y0'], grid['z0']
        nx, ny, nz = grid['nx'], grid['ny'], grid['nz']
        dx, dy, dz = grid['dx'], grid['dy'], grid['dz']
        grid = [x0, y0, z0, nx, ny, nz, dx, dy, dz]
        pv_grid = Visualizer._construct_pyvista_grid(*grid)
        return pv_grid
    
    @requires_pyvista()
    def pv_show(self,
                simulations: list,
                features: list,
                settings: list = [],
                cpos: Union[str, list] = 'xz',
                window_size: tuple = None,
                zoom: float = None,
                savefig: bool = False,
                filename: str = None,
                return_plotter: bool = False,
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
        
        ######################
        ### Initialization ###
        ######################
        
        # Insert default settings values
        pv_settings = []
        for settings_dict in settings:
            settings_dict = Visualizer._set_default_pv_settings(settings_dict)
            pv_settings.append(copy.deepcopy(settings_dict))
        
        ### Construct the plotter
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
        
        plotter = Visualizer._construct_plotter(shape=(rows, columns),
                                                border=True)
        
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
                
                plotter = self._prepare_plotter(plotter,
                                                feature,
                                                n_sim,
                                                settings,
                                                (row, column))
        else:
            for (row, (n_sim, settings)) in enumerate(zip(simulations,
                                                          pv_settings)
                                                      ):
                for (column, feature) in enumerate(features):
                    
                    plotter = self._prepare_plotter(plotter,
                                                    feature,
                                                    n_sim,
                                                    settings,
                                                    (row, column))

        plotter.link_views()
        
        # Show the plotter
        result = Visualizer._show_plotter(plotter,
                                          cpos,
                                          window_size,
                                          zoom,
                                          savefig,
                                          filename,
                                          return_plotter)
        return result
    
    def _prepare_plotter(self,
                         plotter,
                         feature: str,
                         n_sim: int,
                         settings: dict,
                         position: tuple,
                         ):
        """ DOC """
        text_options = settings.pop('text_options', {})
        
        # Test feature validity
        if Visualizer._test_feature_validity(feature):
            text = 'Simulation {} - {}'.format(n_sim, feature)
        else:
            text = 'Simulation {} - {} : Invalid feature'.format(n_sim,
                                                                 feature)
        
        # Retrieve simulation data
        n_iteration = settings.pop('n_iteration', -1)
        simulation_data = self._retrieve_simulation_data(feature,
                                                         n_sim,
                                                         n_iteration)
        
        # Fill the plotter
        row, column = position
        plotter.subplot(row, column)
        plotter = Visualizer._fill_plotter(plotter,
                                           simulation_data,
                                           settings)
        
        # Print text
        text_options.setdefault('text', text)
        if text_options['text'] is not None:
            plotter.add_text(**text_options)
            
        return plotter
    
    @staticmethod
    def _test_feature_validity(feature: str) -> bool:
        """ DOC """
        # Test if asked feature is valid
        if feature not in AUTHORIZED_FEATURES:
            msg = "Asked feature '{}' is invalid. Authorized features: {}."
            msg = msg.format(feature, AUTHORIZED_FEATURES)
            warnings.warn(msg)
            return False
        else:
            return True
            
    def _retrieve_simulation_data(self,
                                  feature: str,
                                  n_sim: int,
                                  n_iteration: int,
                                  ) -> dict:
        """ DOC """
        # Retrieve the simulation data
        simulation_data = self.project._get_simulation_data(n_sim)
        
        # Get the pyvista grid
        simulation_data['grid'] = self.project.grid.data_volume
        simulation_data['mesh'] = self.pv_grid.copy()
        
        # Get data from feature
        feature_data = Visualizer._get_data_from_dict(simulation_data,
                                                      feature,
                                                      n_iteration)
        simulation_data['data'] = feature_data
        simulation_data['meshgrids'] = self.project.grid.get_meshgrids()
        
        return simulation_data
    
    @staticmethod
    def _set_default_pv_settings(settings: dict) -> dict:
        """
        To a dictionary, set the default settings for pyvista.
        """
        for (parameter, value) in Visualizer.DEFAULT_PV_SETTINGS.items():
            settings.setdefault(parameter, value)
        return settings

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
        try:
            # Select the adequate way to retrieve data
            if feature in GEOLOGICAL_FEATURES:
                featured_data = simulation_data[feature].data_volume
            elif feature in DOMAIN_FEATURES:
                featured_data = getattr(simulation_data['domain'], feature)
                featured_data = featured_data.data_volume
            elif feature in ANISOTROPIC_FEATURES:
                featured_data = simulation_data['maps'][feature][n_iteration]

        # except AttributeError:
        #     pass
                
        # except IndexError:
        #     pass
            
        except Exception:
            print('error')
            featured_data = simulation_data['grid'].copy()
            
        return featured_data
    
    @staticmethod
    @requires_pyvista()
    def _df_to_cloud(points: pd.DataFrame):  # pv.PolyData
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
    def _get_surface(simulation_data: dict,
                     surface_name: str
                     ):  # -> Union[pv.StructuredGrid, None]
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
            'water_table', 'bedrock'.

        Returns
        -------
        Union[pv.StructuredGrid, None]
            The pyvista surface object. ``None`` if the surface does not exist.
        """
        X, Y, Z = simulation_data['meshgrids']
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
    def _fractures_to_polygons(fractures: list):  # -> list[pv.PolyData]
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
    
    @requires_pyvista()
    def create_gif(self,
                   filename: str,
                   n_sim: int = -1,
                   feature: str = 'karst',
                   settings: dict = {},
                   window_size: list[int] = [1024, 768],
                   background_color: str = 'white',
                   background_color_top: str = None,
                   zoom: float = 1,
                   fps: int = 10,
                   n_points: int = 24,
                   ) -> None:
        """
        TODO
        
        colormap: str
            https://matplotlib.org/stable/tutorials/colors/colormaps.html
        """
        
        ### Method based on those examples:
        # https://docs.pyvista.org/examples/02-plot/orbit.html#orbiting
        
        ### Get the simulation data
        simulation_data = self.project._get_simulation_data(n_sim)
        
        ### Create the plotter
        plotter = pv.Plotter(off_screen=True, window_size=window_size)
        
        ### Fill the plotter
        plotter = self._fill_plotter(plotter=plotter,
                                     simulation_data=simulation_data,
                                     feature=feature,
                                     settings=settings)
        
        # Set background color
        if background_color_top is None:
            background_color_top = background_color
        plotter.set_background(background_color)
        
        ### Set the initial camera position
        plotter.camera.zoom(zoom)
        # plotter.camera.roll = 0
        plotter.camera.elevation = 0
        plotter.camera.azimuth = 0
        
        ### Generate the GIF
        
        # Open a gif
        plotter.open_gif(filename, fps=fps)
        
        # Create the camera positions
        azimuth_step = 360 / n_points
        
        # Loops
        for _ in range(n_points):
            
            # Update camera position
            # plotter.camera.roll
            # plotter.camera.elevation
            plotter.camera.azimuth += azimuth_step
            
            # Update the background
            
            # Write the frame in the gif
            plotter.write_frame()
        
        # Close and finalize movie
        plotter.close()
        return None
