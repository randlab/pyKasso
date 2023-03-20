"""
TODO
"""
# todo - automatic update when visualization class is already created

####################
### Dependencies ###
####################
import pickle

### External dependencies
import yaml
import numpy as np

import pyvista as pv # ???
from pykasso.visualization import _pyvista

################################
### Visualization sub-module ###
################################

AUTHORIZED_FEATURES = [
    'grid',
    'domain',
    'delimitation',
    'topography',
    'bedrock',
    'water_level',
    'geology',
    # 'beddings' - TODO
    'faults',
    'fractures',
    'cost',
    'alpha',
    'beta',
    'time',
    'karst',
]

class Visualization():
    """"""
    
    ######################
    ### Initialization ###
    ######################
    
    def __init__(self, project_state:str) -> None:
        """"""
        with open(project_state, 'r') as f:
            self.project_state = yaml.safe_load(f)
            self.grid = self._get_grid()
            
    def _get_grid(self):
        grid_settings = self.project_state['grid']
        grid = pv.UniformGrid()
        grid.dimensions = np.array((grid_settings['nx'], grid_settings['ny'], grid_settings['nz'])) + 1
        grid.origin     = (grid_settings['x0'] - grid_settings['dx']/2, grid_settings['y0'] - grid_settings['dy']/2, grid_settings['z0'] - grid_settings['dz']/2)
        grid.spacing    = (grid_settings['dx'], grid_settings['dy'], grid_settings['dz'])
        grid = grid.cast_to_unstructured_grid()
        return grid
    
    #################
    ### Show data ###
    #################
    
    def show(self, simulations:list, features:list, settings={}):
        """"""
        # Sets default settings values
        default_settings = {
            'ghosts'   : [],
            'surfaces' : [],
            'outline'  : True,
            'grid'     : False,
            'inlets'   : False,
            'outlets'  : False,
            'cpos'     : None,
        }
        for key, value in default_settings.items():
            if key not in settings:
                settings[key] = value
        
        plotter = pv.Plotter(shape=(len(features), len(simulations)), border=True)
        ACTORS = []
        
        for i, n_simulation in enumerate(simulations):
            
            for j, feature in enumerate(features):
                
                simulation_data = self._get_simulation_data(n_simulation)
                
                plotter.subplot(j, i)
                text = 'Simulation {} - {}'.format(n_simulation, feature)
                plotter.add_text(text, font_size=20)
                
                actors = self._show_feature(simulation_data, feature, settings)
                for actor in actors:
                    plotter.add_actor(actor, reset_camera=True)
        
        plotter.link_views()      
        plotter.show(cpos=settings['cpos'])

        return None
    
    ###########
    ### GIF ###
    ###########

    def create_gif(self, simulation:int, feature:str, location:str, zoom:float=1, ghosts:list=[], n_points:int=24, fps:int=10, window_size=[1024, 768], background_color:str='white') -> None:
        
        ### Method based on those examples: https://docs.pyvista.org/examples/02-plot/orbit.html#orbiting
        
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
        path = plotter.generate_orbital_path(n_points=n_points, shift=mesh.length)
        plotter.open_gif(location, fps=fps)
        plotter.orbit_on_path(path, write_frames=True)
        plotter.close()
        
        return None
    
    #######################
    ### Private methods ###
    #######################
    
    def _get_simulation_data(self, n:int) -> dict:
        simulation_data_path = self.project_state['simulation_locations'][n] + 'results.pickle'
        simulation_data      = self._read_pickle(simulation_data_path)
        return simulation_data
        
    def _show_feature(self, simulation, feature:str, settings:dict={}):
        
        # Gets the data
        feature_data = self._get_data_from_feature(simulation, feature)
        
        # Gets the domain
        pass
        
        # Ghosts the data
        if len(settings['ghosts']) > 0:
            feature_data = self._get_ghosted_data(feature_data, settings['ghosts'])
           
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
        if settings['inlets']:
            inlets = self._get_points(simulation['inlets'])
            _ = plotter.add_points(inlets, render_points_as_spheres=False, point_size=20, color='r')
            actors.append(_)
        if settings['outlets']:
            outlets = self._get_points(simulation['outlets'])
            _ = plotter.add_points(outlets, render_points_as_spheres=False, point_size=20, color='b')
            actors.append(_)
            
        # Plots the surfaces
        if len(settings['surfaces']) > 0:
            for surface in settings['surfaces']:
                surface = self._get_surface(simulation, surface)
                _ = plotter.add_mesh(surface, smooth_shading=True, opacity=0.66, scalars='data')
                actors.append(_)
        
        # Plots the data
        _ = plotter.add_mesh(feature_data.copy(), scalars='data', scalar_bar_args={'title': 'Vol1'})
        actors.append(_)
        _ = plotter.add_scalar_bar()
        actors.append(_)

        return actors
    
    def _get_data_from_feature(self, simulation, feature:str, iteration:int=-1):
        """"""
        mesh = self.grid
        
        if feature in ['cost', 'alpha', 'beta', 'time', 'karst']:
            data = simulation['maps'][feature][iteration]
        elif feature in ['delimitation', 'topography', 'bedrock', 'water_level']:
            feature_object = getattr(simulation['domain'], feature)
            data = getattr(feature_object, 'data_volume')
        else:
            feature_object = simulation[feature]
            data = getattr(feature_object, 'data_volume')
            
        mesh.cell_data['data'] = data.flatten(order="F")
        return mesh
    
    def _get_ghosted_data(self, data, ghosts:list):
        """"""        
        ghosted_cells = np.argwhere(np.isin(data["data"], ghosts))
        data = data.remove_cells(ghosted_cells)
        return data
    
    def _get_surface(self, simulation, surface:str):
        """"""
        grid = simulation['grid']
        x = grid.X[:,:,0]
        y = grid.Y[:,:,0]
        
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
        cloud  = pv.wrap(points)
        cloud['labels'] = labels
        return cloud
    
    def _read_pickle(self, path:str):
        """"""
        with open(path, 'rb') as handle:
            return pickle.load(handle)
        
    ####### DELETE ???
    # def compare(self, simulations:list, settings:dict={}):
        
    #     plotter = pv.Plotter(shape=(1, len(simulations)), border=True)
        
    #     for i in simulations:
            
    #         simulation_data_path = self.project_state['simulation_locations'][i] + 'results.pickle'
    #         simulation_data      = self._read_pickle(simulation_data_path)
    #         plotter.subplot(0, i)
    #         plotter.add_text(str(i), font_size=24)
    #         actor = _pyvista._show_data(simulation_data, 'karst', {'ghost':[0]}, show=False)
    #         plotter.add_actor(actor, reset_camera=True)
        
    #     return None
        
    
        
############################################

def show(environment, feature, engine='pyvista', settings={}):
    """
    TODO
    """
    ### Selects appropriates functions according to selected engine
    if engine == 'matplotlib':
        from pykasso.visualization import _matplotlib
        _show_data = _matplotlib._show_data
    elif engine == 'pyvista':
        from pykasso.visualization import _pyvista
        _show_data = _pyvista._show_data
    else:
        raise ValueError("ERROR : selected engine not recognised")

    ### Controls validity of settings
    attributes = ['domain', 'inlets', 'outlets', 'tracers']
    for attribute in attributes:
        if attribute in settings:
            if not hasattr(environment, attribute):
                print('WARNING : No {} available'.format(attribute))
                del settings[attribute]
            else:
                if getattr(environment, attribute) is None:
                    print('WARNING : No {} available'.format(attribute))
                    del settings[attribute]

    ### Plot feature
    if _is_feature_valid(environment, feature):
        _show_data(environment, feature, settings)

    return None

def _is_feature_valid(environment, feature):
    """
    TODO
    """
    if feature not in AUTHORIZED_FEATURES:
        raise ValueError("ERROR : selected feature has been not recognized (authorized features : {})".format(AUTHORIZED_FEATURES))

    if feature in ['grid','domain','geology','faults','fractures']:

        if not hasattr(environment, feature):
            raise ValueError("ERROR : environment has no '{}' attribute.".format(feature))
    
        if getattr(environment, feature) is None:
            raise ValueError("ERROR : '{}' attribute is type of None.".format(feature))
        
    elif feature in ['delimitation', 'topography', 'bedrock', 'water_level']:
        
        if not hasattr(environment, 'domain'):
            raise ValueError("ERROR : environment has no '{}' attribute.".format(feature))
        
        if getattr(environment, 'domain') is None:
            raise ValueError("ERROR : '{}' attribute is type of None.".format(feature))
        
        if not hasattr(environment.domain, feature):
            raise ValueError("ERROR : domain's environment has no '{}' attribute.".format(feature))
        
        if getattr(environment.domain, feature) is None:
            raise ValueError("ERROR : '{}' domain's attribute is type of None.".format(feature))
    
    elif feature in ['cost','alpha','beta','time','karst']:

        if not hasattr(environment, 'maps'):
            raise ValueError("ERROR : environment has no 'maps' attribute.")
        
        if feature in ['alpha','beta']:
            if environment.fmm['algorithm'] != 'Riemann3':
                raise ValueError("ERROR : environment has no '{}' data.".format(feature))
            
    return True



def show_array(array, engine='pyvista', ghost=False):
    """
    TODO
    ARRAY VISUALIZATION
    """
    if engine == 'pyvista':
        from pykasso.visualization import _pyvista
        _pyvista._show_array(array, ghost)
    
    
    return None
    
###################
### DEBUG PLOTS ###
###################

def debug(environment, step, engine='matplotlib', settings={}):
    """
    TODO
    """
    if engine == 'pyvista':
        from pykasso.visualization import _pyvista

        if step == 'model':
            _pyvista._debug_plot_model(environment, settings)
        elif step == 'fmm':
            _pyvista._debug_plot_fmm(environment, settings)

    return None







#     #############################
#     ### Visualization methods ###
#     #############################
#
#     def show_catchment(self, label='geology', title=None, cmap='binary'):
#         """
#         Show the entire study domain.
#
#         Parameters
#         ----------
#         label : str, optional
#             Data to show : 'geology', 'topography', 'orientationx', 'orientationy' 'faults' or 'fractures'.
#             By default : 'geology'.
#         title : str, optional
#             Title of the plot. If 'None', 'data' becomes the label.
#         cmap : str, optional
#             Color map, 'binary' by default.
#         """
#         import matplotlib.patches as mp
#         fig, ax1 = plt.subplots()
#         #if title is None:
#         #    title = label
#         #fig.suptitle(title, fontsize=16)
#
#         # Load data to show
#         try:
#             data = [data.data[:,:,0] for data in self.geology.data if data.label==label][-1]
#         except:
#             print('no data for indicated label parameter')
#
#         im1 = ax1.imshow(data.T, origin="lower", extent=self.grid.extent, cmap=cmap)
#
#         fig.colorbar(im1, ax=ax1)
#         if self.settings['data_has_mask']:
#             import matplotlib.patches as mp
#             p1 = mp.PathPatch(self.mask.polygon, lw=2, fill=0, edgecolor='red', label='mask')
#             ax1.add_patch(p1)
#
#         # Points
#         for pts in self.points.points:
#             x, y = zip(*pts.points)
#             ax1.plot(x, y, 'o', label=pts.points_key)
#
#         ax1.set_aspect('equal', 'box')
#         plt.legend(loc='upper right')
#         plt.show()
#         return fig
#
#     def _show_maps(self, sim=-1, iteration=-1, cmap='binary'):
#         """
#         Show the simulated karst network as an image.
#         """
#         karst_network = self.karst_simulations[sim]
#
#         fig, ([ax1,ax2],[ax3,ax4]) = plt.subplots(2, 2, sharex=True, sharey=True)
#         fig.suptitle('Karst Network', fontsize=16)
#
#         ax1.imshow(karst_network.maps['outlets'], extent=self.grid.extent, origin='lower', cmap=cmap)
#         ax1.set_title('Outlets')
#
#         ax2.imshow(karst_network.maps['cost'][iteration], extent=self.grid.extent, origin='lower', cmap=cmap)
#         ax2.set_title('Cost')
#
#         ax3.imshow(karst_network.maps['time'][iteration], extent=self.grid.extent, origin='lower', cmap=cmap)
#         ax3.set_title('Time')
#
#         ax4.imshow(karst_network.maps['karst'][iteration], extent=self.grid.extent, origin='lower', cmap=cmap)
#         ax4.set_title('Karst')
#
#         fig.subplots_adjust(hspace=0.5)
#         plt.show()
#         return None
#
#
#     def show(self, data=None, title=None):
#         """
#         Show the entire study domain (defaults to showing most recent simulation).
#         """
#         if data is None:
#             data = self.karst_simulations[-1]
#
#         fig = plt.figure(figsize=(20,10))
#
#         # Cost map
#         fig.add_subplot(131, aspect='equal')
#         d = data.maps['cost'][-1]
#         plt.xlabel('Cost array'+str(d.shape))
#         d = np.transpose(d, (1,0)) # imshow read MxN and we have NxM
#         plt.imshow(d, extent=self.grid.extent, origin='lower', cmap='gray') #darker=slower
#         plt.colorbar(fraction=0.046, pad=0.04)
#
#         # Travel time map
#         fig.add_subplot(132, aspect='equal')
#         d = data.maps['time'][-1]
#         plt.xlabel('Travel time array'+str(d.shape))
#         d = np.transpose(d, (1,0)) # imshow read MxN and we have NxM
#         plt.imshow(d, extent=self.grid.extent, origin='lower', cmap='cividis') #darker=faster
#         plt.colorbar(fraction=0.046, pad=0.04)
#
#         # Karst map
#         fig.add_subplot(133, aspect='equal')
#         d = data.maps['karst'][-1]
#         plt.xlabel('Karst array'+str(d.shape))
#         d = np.transpose(d, (1,0)) # imshow read MxN and we have NxM
#         plt.imshow(d, extent=self.grid.extent, origin='lower', cmap='gray_r') #darker=conduits
#         plt.colorbar(fraction=0.046, pad=0.04)
#         i = plt.scatter(data.points['inlets'].x,  data.points['inlets'].y,  c='orange')
#         o = plt.scatter(data.points['outlets'].x, data.points['outlets'].y, c='steelblue')
#         p = matplotlib.patches.Rectangle((0,0),0,0, ec='r', fc='none')
#         if self.settings['data_has_mask']:
#             closed_polygon = self.mask.vertices[:]
#             closed_polygon.append(closed_polygon[0])
#             x,y = zip(*closed_polygon)
#             plt.plot(x,y, color='red', label='mask')
#         #plt.legend([i,o,p], ['inlets', 'outlets', 'catchment'], loc='upper right')
#         plt.legend([i,o], ['inlets', 'outlets'], loc='upper right')
#
#         if title is not None:
#             fig.suptitle(title, fontsize=16)
#         plt.show()
#         return fig
#
#     def show_network(self, data=None, simplify=False, ax=None, plot_nodes=True, mask=True, labels=['inlets', 'outlets'], title=None, cmap=None, color='k', alpha=1, legend=True):
#         """
#         #Chloe: This is a new function that I use to create all the figures for the paper.
#         Show the karst network as a graph with nodes and edges. Defaults to showing latest iteration.
#
#         Parameters
#         ----------
#         data:
#             karst simulation object containing nodes, edges, points, etc. Can be obtained from self.karst_simulations[i]
#         ax :
#             axis to plot on
#         label :
#             None or list of strings ['nodes','edges','inlets','outlets'], indicating which components to label
#         title : str
#             title of plot
#         cmap : str
#             colormap to use when plotting
#         color : str
#             single color to use when plotting (cannot have both cmap and color)
#         alpha : float
#             opacity to plot with (1=opaque, 0=transparent)
#         legend : bool
#             whether to display legend
#         plot_nodes : bool
#             whether to display nodes
#         polygon : bool
#             whether to display the bounding polygon
#         """
#
#         if ax == None:
#             fig,ax = plt.subplots(figsize=(10,10))
#             ax.set_aspect('equal')
#
#         if data == None:
#             data = self.karst_simulations[-1]
#
#         if mask == True:
#             if self.settings['data_has_mask']:
#                 closed_polygon = self.mask.vertices[:]
#                 closed_polygon.append(closed_polygon[0])
#                 x,y = zip(*closed_polygon)
#                 ax.plot(x,y, color='maroon')
#                 p = matplotlib.lines.Line2D([0],[0], color='k')
#
#         if simplify == True:
#             nodes = data.network['nodes']   #get all nodes
#             nodes_simple = data.network['karstnet'].graph_simpl.nodes  #get indices of only the nodes in the simplified graph
#             nodes_simple = {key: nodes[key] for key in nodes_simple}   #make df of only the nodes in the simplified graph, for plotting
#             edges = data.network['edges']   #get all edges
#             edges_simple = data.network['karstnet'].graph_simpl.edges  #get only the edges in the simplified graph
#             edges_simple = {i: edge for i,edge in enumerate(edges_simple)}   #make df of only the edges in the simplified graph, for p
#             nodes = pd.DataFrame.from_dict(nodes_simple, orient='index', columns=['x','y','type']) #convert to pandas for easier plotting
#             edges = pd.DataFrame.from_dict(edges_simple, orient='index', columns=['inNode','outNode'])
#         else:
#             nodes = pd.DataFrame.from_dict(data.network['nodes'], orient='index', columns=['x','y','type']) #convert to pandas for easier plotting
#             edges = pd.DataFrame.from_dict(data.network['edges'], orient='index', columns=['inNode','outNode'])
#
#         #Set up data for plotting:
#         fromX = nodes.x.loc[edges.inNode]      #calculate coordinates for link start and end points
#         fromY = nodes.y.loc[edges.inNode]
#         toX   = nodes.x.loc[edges.outNode]
#         toY   = nodes.y.loc[edges.outNode]
#
#         #Plot nodes and edges:
#         if plot_nodes:
#             n = ax.scatter(nodes.x,              nodes.y,                  c='k',         alpha=alpha, s=5)  #scatterplot nodes
#         i = ax.scatter(data.points['inlets'].x,  data.points['inlets'].y,  c='orange',    s=30) #scatterplot inlets
#         o = ax.scatter(data.points['outlets'].x, data.points['outlets'].y, c='steelblue', s=30) #scatterplot outlets
#         e = matplotlib.lines.Line2D([0],[0])                                                  #line artist for legend
#         for ind in edges.index:                                                               #loop over edge indices
#             if cmap is not None:
#                 ax.plot((fromX.iloc[ind], toX.iloc[ind]), (fromY.iloc[ind], toY.iloc[ind]), c=plt.cm.get_cmap(cmap)(ind/len(edges)), alpha=alpha)  #plot each edge, moving along color gradient to show order
#             elif color is not None:
#                 ax.plot((fromX.iloc[ind], toX.iloc[ind]), (fromY.iloc[ind], toY.iloc[ind]), c=color, alpha=alpha)  #plot each edge in same color
#
#         #Add labels:
#         if labels == None:
#             pass
#         else:
#             if 'nodes' in labels:                                         #label node indices
#                 for ind in nodes.index:                                   #loop over node indices
#                     ax.annotate(str(ind), xy=(nodes.y[ind]-10, nodes.x[ind]))  #annotate slightly to left of each node
#             if 'edges' in labels:
#                 for ind in edges.index:
#                     ax.annotate(str(ind), xy=(edges.y[ind]-10, edges.x[ind]))  #annotate slightly to left of each edge
#             if 'inlets' in labels:
#                 for index,inlet in data.points['inlets'].iterrows():
#                     ax.annotate(str(int(inlet.outlet))+'-'+str(int(inlet.inlet_iteration)),  xy=(inlet.x-(6*self.grid.dx),  inlet.y))
#             if 'outlets' in labels:
#                 for index,outlet in data.points['outlets'].iterrows():
#                     ax.annotate(str(int(outlet.name)), xy=(outlet.x-(4*self.grid.dx), outlet.y))
#
#         #Add legend & title:
#         if legend:
#             if plot_nodes:
#                 if plot_polygon:
#                     ax.legend([i,o,n,e,p],['inlets','outlets','nodes','edges','mask'])
#                 else:
#                     ax.legend([i,o,n,e],['inlets','outlets','nodes','edges'])
#             else:
#                 if plot_polygon:
#                     ax.legend([i,o,e,p],['inlets','outlets','edges','mask'])
#                 else:
#                     ax.legend([i,o,e],['inlets','outlets','edges','mask'])
#         if title is not None:
#             ax.set_title(title, fontsize=16)
#
#         return None