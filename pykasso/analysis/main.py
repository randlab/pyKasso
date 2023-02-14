####################
### Dependencies ###
####################
import os
import copy
import pickle

### External dependencies
import yaml
import numpy as np
import pandas as pd
import karstnet as kn

###########################
### Analysis sub-module ###
###########################

class Analysis():
    """"""
    def __init__(self, project_state:str) -> None:
        """"""
        with open(project_state, 'r') as f:
            self.project_state = yaml.safe_load(f)
            self.statistics    = self._load_statistical_reference()
            
    def _load_statistical_reference(self):
        """Reference metrics for statistical karstic network analysis."""
        usecols = None
        # usecols = "B:I,K"
        return pd.read_excel(os.path.dirname(os.path.abspath(__file__)) + '/../../misc/' + 'statistics.xlsx', usecols=usecols).describe()

    def compute_karstnet_metrics(self, verbosity=0):
        """"""
        df_metrics = pd.DataFrame()
        
        for i, path in enumerate(self.project_state['simulation_locations']):
            
            # Extracts data
            results = self._read_pickle(path + 'results.pickle')
            karstnet_edges = list(results['vectors']['edges'].values())     # Converts to format needed by karstnet (list)
            karstnet_nodes = copy.deepcopy(results['vectors']['nodes'])     # Converts to format needed by karstnet (dic with only coordinates) - make sure to only modify a copy!
            for key, value in karstnet_nodes.items():                       # Drops last item in list (the node type) for each dictionary entry
                value.pop()

            # Computes karstnet metrics
            k = kn.KGraph(karstnet_edges, karstnet_nodes)  # Makes graph - edges must be a list, and nodes must be a dic of format {nodeindex: [x,y]}
            metrics = k.characterize_graph(verbosity)
            
            # Concatenates dataframes
            df_ = pd.DataFrame(metrics, index=[i])
            df_metrics = pd.concat([df_metrics, df_])
            
        return df_metrics
    
    def compare_karstnet_metrics(self, dataframe):
        """"""
        if isinstance(dataframe, pd.core.series.Series):
            dataframe = dataframe.to_frame().T
        
        # Defines background/text coloring function             (# return 'background-color: red')
        def _bg_color(x, min_val, max_val, mean_val, std_val):
            if (x < min_val) or (x > max_val):
                return 'color: red'
            else:
                inf = mean_val - std_val
                sup = mean_val + std_val
                if (x < inf) or (x > sup): 
                    return 'color: #FFFF00'
                else:
                    return 'color: #00FF00'
        
        # Iterates in the dataframe columns
        styler = dataframe.style
        for column_name in dataframe:
            styler = styler.applymap(
                _bg_color, 
                min_val  = self.statistics[column_name]['min'], 
                max_val  = self.statistics[column_name]['max'],
                mean_val = self.statistics[column_name]['mean'],
                std_val  = self.statistics[column_name]['std'], 
                subset = [column_name])

        return styler
    
    def compute_average_karstic_network(self):
        """Compute the mean of all the simulations."""
        nx, ny, nz = self.project_state['grid']['nx'], self.project_state['grid']['ny'], self.project_state['grid']['nz']
        karst_map = np.zeros((nx, ny, nz))
        for path in self.project_state['simulation_locations']:
            results = self._read_pickle(path + 'results.pickle')
            karst_map = np.add(karst_map, results['maps']['karst'][-1]).copy()
        return karst_map / self.project_state['n_simulation']
    
    
    def _read_pickle(self, path:str):
        """"""
        with open(path, 'rb') as handle:
            return pickle.load(handle)
           
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