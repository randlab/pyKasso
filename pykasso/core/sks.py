# ############
# ### TODO ###
# ############
# # 4 - Set rand_seed
# # 6 - Better 'debug'/'verbosity' management
# # LT - Methods documentation
# # Anisotropic fast-marching
# #   - How to compute alpha map ?
# #   - How to compute beta map ?

# export des points : inlets et outlets
# import copy
# import karstnet as kn
# import numpy.ma as ma

# ### Fast-Marching package
# import agd
# from agd import Eikonal
# from agd.Metrics import Riemann
import os
import shutil
import logging
# from functools import wraps

### Local dependencies
from ._validations import validate_core_settings, validate_sks_settings, validate_sim_settings
from .grid import Grid
from .mask import Mask
from .geologic_feature import Geology, Topography, Karst, Field, Faults, Fractures
from .points import generate_random_points, inspect_points_validity

### External dependencies
import yaml
import numpy as np
import pandas as pd




# TODO
# Visualiser topography
# Visualisateur : couches visibles en argument
# Tester la coupe des fractures en fonction de la topography


### SKS class
class SKS():
    """
    Class storing all the parameters, the data and the resulting simulated stochastic karst networks.

    TODO
    """
    def __init__(self, project_directory:str):
        
        # Initializes SKS settings
        self.CORE_SETTINGS = {}
        self.SKS_SETTINGS  = {}
        self.SIM_SETTINGS  = {}
        self.SIMULATIONS   = pd.DataFrame()

        # Creates the project directories
        list_directories = [
            project_directory,
            project_directory + '/inputs/',
            project_directory + '/outputs/',
            project_directory + '/settings/',
        ]
        [os.makedirs(directory, exist_ok=True) for directory in list_directories]

        # Populates the CORE_SETTINGS dictionnary
        self.CORE_SETTINGS['project_directory']  = list_directories[0]
        self.CORE_SETTINGS['inputs_directory']   = list_directories[1]
        self.CORE_SETTINGS['outputs_directory']  = list_directories[2]
        self.CORE_SETTINGS['settings_directory'] = list_directories[3]
        self.CORE_SETTINGS['core_settings_filename'] = 'CORE_SETTINGS.yaml'
        self.CORE_SETTINGS['sks_settings_filename']  = 'SKS_SETTINGS.yaml'
        self.CORE_SETTINGS['sim_settings_filename']  = 'SIM_SETTINGS.yaml'

        ### Populates the settings directory with initial yaml settings file
        misc_directory = os.path.dirname(os.path.abspath(__file__)) + '/../../misc/'

        # CORE_SETTINGS
        core_settings_yaml = self.CORE_SETTINGS['settings_directory'] + self.CORE_SETTINGS['core_settings_filename']
        with open(core_settings_yaml, 'w') as f:
            text = "# pyKasso CORE SETTINGS \n---\n"
            yaml_text = yaml.safe_dump(self.CORE_SETTINGS)
            text = text + yaml_text + "..."
            f.write(text)

        # SKS_SETTINGS
        shutil.copy2(misc_directory + self.CORE_SETTINGS['sks_settings_filename'], project_directory + '/settings/')

        # SIM_SETTINGS
        shutil.copy2(misc_directory + self.CORE_SETTINGS['sim_settings_filename'], project_directory + '/settings/')

        ### Creates a log file
        log_filename = 'pyKasso.log'
        filename = self.CORE_SETTINGS['outputs_directory'] + log_filename
        if os.path.exists(filename):
            try:
                os.remove(filename)
            except:
                pass
        logging.basicConfig(filename=filename, encoding='utf-8', level=logging.INFO, filemode="w")
        logger = logging.getLogger("sks")
        logger.info('pyKasso project initialized')


##########################################################################################################################################################################
### LOADING ###
###############

    def load_settings(self):
        """
        TODO
        """
        self.CORE_SETTINGS = self._load_settings_file('core_settings_filename', validate_core_settings)
        self.SKS_SETTINGS  = self._load_settings_file('sks_settings_filename' , validate_sks_settings)
        self.SIM_SETTINGS  = self._load_settings_file('sim_settings_filename' , validate_sim_settings)
        return None

    def _load_settings_file(self, file, validation_function):
        """
        TODO
        """
        file_destination = self.CORE_SETTINGS['settings_directory'] + self.CORE_SETTINGS[file]
        with open(file_destination, 'r') as f:
            unvalidated_settings = yaml.safe_load(f)
            validated_settings = validation_function(unvalidated_settings)
        return validated_settings

##########################################################################################################################################################################
### WRAPPERS ###
################

    def _decorator_logging(test):
        def _(function):
            # print(function.__name__)
            # @wraps(function)
            def _wrapper_logging(*args, **kwargs):
                logger = logging.getLogger("sks")
                try:
                    result = function(*args, **kwargs)
                except:
                    logger.error(function.__name__)
                    raise
                else:
                    logger.info(function.__name__)
                    return result
                finally:
                    pass
            return _wrapper_logging
        return _

##########################################################################################################################################################################
### CONSTRUCTING STATIC AND DYNAMIC FEATURES ###
################################################

    def construct_features(self):
        """
        TODO
        """
        self._construct_static_features()
        # self._construct_dynamic_features()
        return None

    def _construct_static_features(self):
        """
        TODO
        """
        # 1 - Constructs the grid
        self._construct_feature_grid()

        # 2 - Constructs the mask
        self._construct_feature_mask()
        
        # 3 - Constructs the geology
        self._construct_feature_geology()

        # 4 - Constructs the initial karst network
        # self._construct_feature_karst()

        # 5 - Constructs the field
        # self._construct_feature_field()

        return None

    def _construct_dynamic_features(self):
        """
        TODO
        """
        # n - Constructs the points
        # self._construct_feature_points()

        # n - Constructs the faults
        # self._construct_feature_faults()

        # n - Constructs the fractures
        # self._construct_feature_fractures()
        return None

    @_decorator_logging('grid')
    def _construct_feature_grid(self):
        """
        Constructs the grid.
        """
        self.GRID = Grid(**self.SKS_SETTINGS['grid'])
        return None

    @_decorator_logging('mask')
    def _construct_feature_mask(self):
        """
        Constructs the mask.
        """
        if self.SKS_SETTINGS['mask']['data'] != '':
            self.MASK = Mask(**self.SKS_SETTINGS['mask'])
            self.MASK.validate_vertices(self.GRID)
        else:
            self.MASK = None
        return None

    @_decorator_logging('geology')
    def _construct_feature_geology(self):
        """
        Constructs the geology.
        """
        self.GEOLOGY    = Geology(**self.SKS_SETTINGS['geology'], grid=self.GRID)
        topography      = np.where(self.GEOLOGY.data > 0, 1, 0)
        self.TOPOGRAPHY = Topography(name='Topography', data=topography, grid=self.GRID)
        self.GEOLOGY._compute_surface(self.TOPOGRAPHY.data)
        return None

    @_decorator_logging('karst')
    def _construct_feature_karst(self):
        """
        Constructs the initial karst conduits network.
        """
        if self.SKS_SETTINGS['karst']['data'] != '':
            self.KARST = Karst(**self.SKS_SETTINGS['karst'], grid=self.GRID)
            self.KARST.data = np.where(self.TOPOGRAPHY == 1, self.GEOLOGY.data, 0)
            self.KARST._compute_surface(self.TOPOGRAPHY)
        else:
            self.KARST = None
        return None

    @_decorator_logging('field')
    def _construct_feature_field(self):
        """
        Constructs the field.
        """
        if self.SKS_SETTINGS['field']['data'] != '':
            self.FIELD = Field(**self.SKS_SETTINGS['field'], grid=self.GRID)
        else:
            self.FIELD = None
        return None

    ### SIM_SETTINGS

    @_decorator_logging('points')
    def _construct_feature_points(self):
        """
        TODO
        """
        self.POINTS = pd.DataFrame(columns=['label', 'x', 'y', 'z'])
        dtypes_dict = {
            'label' : 'object',
            'x'     : 'float64',
            'y'     : 'float64',
            'z'     : 'float64',
        }
        self.POINTS = self.POINTS.astype(dtypes_dict)

        # Constructs the outlets
        self._construct_feature_io('outlets')

        # Constructs the inlets
        self._construct_feature_io('inlets')

        # n - Constructs the tracers
        # self._construct_feature_tracers()

        return None

    @_decorator_logging('io')
    def _construct_feature_io(self, kind):
        """
        TODO
        Constructs the inlets / outlets.

        Four cases

        1. No points declared
        2. More points required than provided : Generates additional random points
        3. Less points required than provided : Pick random points among provided ones
        4. Points required equals points declared
        """

        # Controls type data input
        # uncontrolled_points = ''
        # if isinstance(self.SIM_SETTINGS[kind]['data'], str):
        #     location = self.SIM_SETTINGS[kind]['data']
        #     uncontrolled_points = np.genfromtxt(location)
        # elif isinstance(self.SIM_SETTINGS[kind]['data'], list):
        #     uncontrolled_points = self.SIM_SETTINGS[kind]['data']

        # # If points have been declared : inspects validity of points
        # if not uncontrolled_points == '':
        #     points = inspect_points_validity(uncontrolled_points)

        # Case 1 - No points declared
        if (self.SIM_SETTINGS[kind]['data'] == '') or (self.SIM_SETTINGS[kind]['data'] == []):
            points = generate_random_points(self.SIM_SETTINGS[kind]['number'], self.GRID, self.MASK, self.GEOLOGY)

        # Case 2 - More points required than provided
        elif (self.SIM_SETTINGS[kind]['number'] > len(self.SIM_SETTINGS[kind]['data'])):
            n_points = self.SIM_SETTINGS[kind]['number'] - len(self.SIM_SETTINGS[kind]['data'])
            points = self.SIM_SETTINGS[kind]['data'] + generate_random_points(n_points, self.GRID, self.MASK, self.GEOLOGY)

        # Case 3 - Less points required than provide
        elif (self.SIM_SETTINGS[kind]['number'] < len(self.SIM_SETTINGS[kind]['data'])):
            rng = np.random.default_rng()
            points = rng.choice(self.SIM_SETTINGS[kind]['data'], self.SIM_SETTINGS[kind]['number'], replace=False)

        # Case 4 - Points required equals points declared
        else:
            points = self.SIM_SETTINGS[kind]['data']

        ### Populates the DataFrame
        label = [kind] * self.SIM_SETTINGS[kind]['number']
        x, y = zip(*points)
        data = {
            'label' : label,
            'x'     : x,
            'y'     : y,
        }
        POINTS = pd.DataFrame(data=data)
        self.POINTS = pd.concat([self.POINTS, POINTS])

        return None

    @_decorator_logging('tracers')
    def _construct_feature_tracers(self):
        """
        Constructs the tracers.
        """
        return None

    @_decorator_logging('faults')
    def _construct_feature_faults(self):
        """
        TODO
        """
        if self.SIM_SETTINGS['faults']['data'] != '':
            self.FAULTS = Faults(**self.SIM_SETTINGS['faults'], grid=self.GRID)
            self.FAULTS.data = np.where(self.TOPOGRAPHY == 1, self.FAULTS.data, 0)
            self.FAULTS._compute_surface(self.TOPOGRAPHY)
        else:
            self.FAULTS = None
        return None

    @_decorator_logging('fractures')
    def _construct_feature_fractures(self):
        """
        TODO
        """
        if self.SIM_SETTINGS['fractures']['data'] != '':
            self.FRACTURES = Fractures(**self.SIM_SETTINGS['fractures'], grid=self.GRID)
            # self.FRACTURES.data = np.where(self.TOPOGRAPHY == 1, self.FAULTS.data, 0)
            # self.FRACTURES._compute_surface(self.TOPOGRAPHY)
        else:
            self.FRACTURES = None
        return None


##########################################################################################################################################################################
### CONSTRUCTING DYNAMIC FEATURES ###
#####################################




##########################################################################################################################################################################

#     

##########################################################################################################################################################################

#     def compute_karst_networks(self):
#         """
#         TODO
#         """
#         return None
















#     def __init__(self, yaml_settings_file=None):
#         """
#         Constructs a SKS class according to the specified settings datafile.

#         Parameters
#         ----------
#         yaml_settings_file : str
#             YAML settings file location.

#         Examples
#         --------
#         >>> catchment = pk.SKS()
#         >>> catchment = pk.SKS('example.yaml')
#         """

#         print('CAUTION: You are using the development version of this package.') # uncomment this to tell apart the development version and the main version

#         # If no yaml file provided, select the default configuration
#         if yaml_settings_file is None:
#             yaml_settings_file = os.path.dirname(os.path.abspath(__file__)) + '/' + 'default_files/settings.yaml'

#         try:
#             with open(yaml_settings_file, 'r') as stream:
#                 try:
#                     settings = yaml.safe_load(stream)
#                 except yaml.YAMLError as exc:
#                     print(exc)
#                     raise
#         except:
#             print("/!\\ Error : unable to read the YAML settings file!")
#             raise

#         self.SETTINGS    = settings
#         self.SIMULATIONS = []

#         # Force some parameters from the yaml configuration file to be float in order to avoid further issues
#         parameters = self.SETTINGS['geologic_features']['fractures']['fractures_settings']
#         for parameter in parameters:
#             parameters[parameter] = [float(elem) for elem in parameters[parameter]]

#         # Load the reference indicators for statistical karstic network analysis
#         self._get_reference_statistics()

#         # TODO 4
#         # Set random seed according to parameters
#         # if self.SETTINGS['rand_seed'] > 0:
#         #     np.random.seed(self.SETTINGS['rand_seed'])
#         ##########################

#         self._initialize_data()

#         return None


#     ######################
#     ### initialization ###
#     ######################

#     # __init__()
#     def _initialize_data(self):
#         """
#         TODO
#         """
#         ### iii. The Geology
#         self.geology = GeologicFeatureManager(self.GRID)

#         for feature in self.SETTINGS['geologic_features']:
#             if self.SETTINGS['geologic_features'][feature]['is_active']: # TODO
#                 settings = self.SETTINGS['geologic_features'][feature]
#                 settings['label'] = feature
#                 self.geology.create_geologic_feature(**settings)
#                 if self.SETTINGS['debug'][feature]:
#                     self.debug_plot(feature)
#             else:
#                 pass

#         ### iv. The Points
#         self.points = PointFeatureManager(self.GRID, self.mask, self.geology)

#         for feature in self.SETTINGS['points']:
#             if feature != 'options':
#                 settings = self.SETTINGS['points'][feature]
#                 settings['label'] = feature
#                 self.points.create_point_feature(**settings)
#         if self.SETTINGS['debug']['points']:
#             self.debug_plot('points')

#         # Shuffle inlets/outlets
#         if self.SETTINGS['points']['options']['inlets_shuffle']:
#             np.random.shuffle(self.points.point_features['inlets'].points)
#         if self.SETTINGS['points']['options']['outlets_shuffle']:
#             np.random.shuffle(self.points.point_features['outlets'].points)

#         return None

#     # __init__()
#     def _get_reference_statistics(self):
#         """
#         Gets the reference statistics for comparing it with karstnet outputs.
#         """
#         path = os.path.dirname(os.path.abspath(__file__)) + '/' + 'default_files/statistics.xlsx'
#         df = pd.read_excel(path, usecols = "B:I,K")
#         self.reference_statistics = df.describe()
#         return None


#     ########################
#     ### Karst simulation ###
#     ########################

#     ### I - Initializes parameters
#     def initialize(self, verbose=False):
#         """
#         TODO
#         """
#         # i. Initializes parameters from gelogic model
#         self._initialize_geologic_model_parameters()
#         # ii. Initializes parameters from farst-marching method
#         self._initialize_karst_network_parameters()

#         if self.SETTINGS['debug']['test']:
#             self.debug_plot_initialize()

#         return None

#     # i. Initializes parameters from gelogic model
#     def _initialize_geologic_model_parameters(self):
#         """
#         TODO
#         """
#         # Creates a dictionnary with the geologic data
#         self.DATA = {
#             'mask'    : None,
#             'geology' : {
#                 'geology'   : None,
#                 'faults'    : None,
#                 'fractures' : None,
#                 'karst'     : None,
#                 'field'     : None,
#             },
#             'points' : {
#                 'inlets'  : None,
#                 'outlets' : None,
#             },
#         }

#         # Populates the dictionnary
#         self.DATA['mask'] = self.mask
#         for key in self.DATA['geology']:
#             if self.geology.geologic_features[key] is not None:
#                 self.DATA['geology'][key] = self.geology.geologic_features[key]
#         for key in self.DATA['points']:
#             self.DATA['points'][key] = self.points.point_features[key].points

#         return None

#     # ii. Initializes parameters from farst-marching method
#     def _initialize_karst_network_parameters(self):
#         """
#         Initialize the karst network parameters.
#         """
#         ### Inlets - Outlets - Iterations

#         # TODO ###################################################
#         # Proceed to assert on some variables
#                 # features = [self.settings['outlets_importance'], self.settings['inlets_importance'], self.settings['inlets_per_outlet']]
#                 # for feature in features:
#                 #     assert isinstance(feature, list)
#                 # assert len(self.settings['outlets_importance']) == len(self.outlets)
#         ###################################################

#         # Defining some variables
#         self.outlets = self.DATA['points']['outlets']
#         self.inlets  = self.DATA['points']['inlets']
#         outlets_nbr = len(self.outlets)
#         inlets_nbr  = len(self.inlets)
#         self.outlets_importance = self.SETTINGS['points']['options']['outlets_importance']
#         self.inlets_importance  = self.SETTINGS['points']['options']['inlets_importance']
#         self.inlets_per_outlet  = self.SETTINGS['points']['options']['inlets_per_outlet']

#         # Calculating inlets and outlets repartitions
#         self.nbr_iteration  = len(self.outlets_importance) * len(self.inlets_importance)         # total number of iterations that will occur
#         outlets_repartition = self._repartition_points(outlets_nbr, self.outlets_importance)     # correct for outlets_importance not summing to correct number of actual outlets
#         inlets_repartition  = self._repartition_points(inlets_nbr , self.inlets_per_outlet)      # correct for inlets_per_outlet not summing to correct number of actual inlets

#         # Distributing inlets and outlets iterations
#         outlets_distribution = pd.Series([k for (k, n) in enumerate(outlets_repartition) for j in range(n)], name='outlet_iteration')
#         inlets_distribution  = pd.Series([k for (k, n) in enumerate(inlets_repartition)  for j in range(n)], name='outlet_key')
#         self._outlets = pd.concat([self.outlets, outlets_distribution], axis=1) # store as a semi-private variable for internal use only
#         self._inlets  = pd.concat([self.inlets , inlets_distribution] , axis=1) # store as a semi-private variable for internal use only

#         # Distributing iterations for each inlet
#         for (outlet_key, row) in self._outlets.iterrows():
#             inlets_test    = self._inlets['outlet_key']==outlet_key
#             inlets_current = self._inlets[inlets_test]
#             inlets_nbr     = len(inlets_current)
#             repartition  = self._repartition_points(inlets_nbr, self.inlets_importance)
#             distribution = pd.Series([k for (k, n) in enumerate(repartition) for j in range(n)], name='inlet_iteration', index=inlets_current.index)
#             self._inlets.loc[inlets_test, 'inlet_iteration'] = distribution


#         ### Raster maps
#         self.maps            = {}
#         self.maps['outlets'] = np.full((self.GRID.nx, self.GRID.ny, self.GRID.nz), np.nan) # map of null values where each cell with an outlet will have the index of that outlet
#         self.maps['nodes']   = np.full((self.GRID.nx, self.GRID.ny, self.GRID.nz), np.nan) # map of null values where each cell that has a node will be updated with that node index
#         self.maps['cost']    = np.zeros((self.nbr_iteration, self.GRID.nx, self.GRID.ny, self.GRID.nz)) # cost of travel through each cell
#         self.maps['alpha']   = np.zeros((self.nbr_iteration, self.GRID.nx, self.GRID.ny, self.GRID.nz)) # cost of travel along gradient through each cell
#         self.maps['beta']    = np.zeros((self.nbr_iteration, self.GRID.nx, self.GRID.ny, self.GRID.nz)) # cost of travel perpendicular to gradient through each cell
#         self.maps['time']    = np.zeros((self.nbr_iteration, self.GRID.nx, self.GRID.ny, self.GRID.nz)) # travel time to outlet from each cell
#         self.maps['karst']   = np.zeros((self.nbr_iteration, self.GRID.nx, self.GRID.ny, self.GRID.nz)) # presence/absence of karst conduit in each cell


#         ### Vector maps
#         self.nodes     = {} #empty dic to store nodes (key: nodeID, val: [x, y, z, type])
#         self.edges     = {} #empty dic to store edges (key: edgeID, val: [inNode, outNode])
#         self.n         = 0  #start node counter at zero
#         self.e         = 0  #start edge counter at zero
#         self.geodesics = [] #empty list to store raw fast-marching path output

#         ### Set up fast-marching:
#         # Note: AGD-HFM library has different indexing, so model dimensions must be [nz, ny, nx],
#         # and model extent must be [zmin,zmax, ymin,ymax, xmin,xmax] (NOT x0,y0,z0)
#         self.riemannMetric = []                    #this changes at every iteration, but cannot be stored?
#         self.fastMarching = agd.Eikonal.dictIn({
#             'model'             : self.SETTINGS['algorithm'],   # set algorithm from settings file ('Isotropic2', 'Isotropic3', 'Riemann2', 'Riemann3')
#             'order'             : 2,                            # recommended setting: 2 (replace by variable)
#             'exportValues'      : 1,                            # export the travel time field
#             'exportGeodesicFlow': 1                             # export the walker paths (i.e. the conduits)
#         })
#         self.fastMarching.SetRect(                              # give the fast-marching algorithm the model grid
#             sides=[[self.GRID.zmin, self.GRID.zmax],
#                    [self.GRID.ymin, self.GRID.ymax],            # leftmost edge, rightmost edge (NOT centerpoint)
#                    [self.GRID.xmin, self.GRID.xmax]],           # bottom edge,   top edge (NOT centerpoint)
#             dims=[self.GRID.nz, self.GRID.ny, self.GRID.nx])    # number of cells, number of cells, number of cells

#         return None

#     # ii
#     def _repartition_points(self, nbr_points, importance):
#         '''Correct for integers in importance factors list not summing correctly to total number of points'''
#         total_importance = float(sum(importance))                                               # total number of points as assigned (this may not be the actual total)
#         proportion       = [float(i)/total_importance       for i in importance]                # percent of points to use this iteration
#         repartition      = [round(proportion[i]*nbr_points) for i in range(len(proportion))]    # number of points to use this iteration (need to round bc percentage will not result in whole number)
#         repartition[-1]  = nbr_points-sum(repartition[0:-1])                                    # leftover points are assignd to last iteration
#         return repartition




#     def compute_karst_network(self, verbose=False):
#         """
#         Compute the karst network according to the parameters.

#         Append the results to the `karst_simulations` attribute as an element of an array.

#         Parameters
#         ----------
#         TODO
#         """
#         ### i. Computes conduits for each generation & store nodes and edges for network
#         self._compute_iterations_karst_network()
#         if self.SETTINGS['debug']['costmap']:
#             self.debug_plot_fmm_feature('costmap')
#         if self.SETTINGS['debug']['timemap']:
#             self.debug_plot_fmm_feature('timemap')
#         if self.SETTINGS['debug']['karstmap']:
#             self.debug_plot_fmm_feature('karstmap')

#         ### ii. Calculates the karst network statistics indicators with karstnet and save karst network
#         karstnet_edges = list(self.edges.values()) #convert to format needed by karstnet (list)
#         karstnet_nodes = copy.deepcopy(self.nodes) #convert to format needed by karstnet (dic with only coordinates) - make sure to only modify a copy!
#         for key, value in karstnet_nodes.items():  #drop last item in list (the node type) for each dictionary entry
#             value.pop()

#         # Computes karstnet indicators
#         k = kn.KGraph(karstnet_edges, karstnet_nodes)  #make graph - edges must be a list, and nodes must be a dic of format {nodeindex: [x,y]}
#         stats = k.characterize_graph(verbose)

#         ### iii. Store all the relevant data for this network in dictionaries:
#         maps = copy.deepcopy(self.maps)                  #store copy of maps
#         points = {}
#         points['inlets']  = copy.deepcopy(self._inlets)  #store copy of inlets in pandas df format with info about iteration and inlet/outlet assignment (this will be used in plotting later)
#         points['outlets'] = copy.deepcopy(self._outlets) #store copy of outlets in pandas df format with info about iteration and inlet/outlet assignment (this will be used in plotting later)
#         network = {}
#         network['edges'] = copy.deepcopy(self.edges) #store copy of edges list
#         network['nodes'] = copy.deepcopy(self.nodes) #store copy of nodes list
#         network['karstnet'] = copy.deepcopy(k)       #store copy of karstnet network object (including graph)
#         config = copy.deepcopy(self.SETTINGS)        #store copy of settings for the run being stored
#         try:
#             del config['debug']
#         except:
#             pass
#         self.SIMULATIONS.append(KarstNetwork(maps, points, network, stats, config))

#         # TODO ??????????
#         ### iv. - Return inlets and outlets to their original format:
#         # #Chloe: this is to convert inlets and outlets back from pandas dataframes
#         # #if the private self._inlets and self._outlets variables are working correctly, this step may be removed
#         # self.inlets  = np.asarray(self._inlets)[:,0:2]    #return inlets to array format without info about iteration and inlet/outlet assignment (important for later iterations)
#         # self.outlets = np.asarray(self._outlets)[:,0:2]   #return outlets to array format without info about iteration and inlet/outlet assignment (important for later iterations)
#         # return None




#     ### ii. Compute conduits
#     def _compute_iterations_karst_network(self):
#         """
#         Compute each generation of karst conduits.

#         TODO : à terminer
#         """
#         if self.SETTINGS['verbosity'] > 0:
#             print('-START-')

#         # Define outlets map according to outlets emplacements
#         for (i, outlet) in self._outlets.iterrows(): #assign outlet indices. Compute the outlets map (array indicating location of outlets as their index and everywhere else as nan).
#             X = self.GRID.get_i(outlet.x)
#             Y = self.GRID.get_j(outlet.y)
#             Z = self.GRID.get_k(outlet.z)
#             self.maps['outlets'][X][Y][Z] = i

#         # Plot for debugging:
#         if self.SETTINGS['verbosity'] > 1:
#             f, ax = plt.subplots(1, 1, figsize=(10, 10))
#             ax.imshow(np.transpose(self.DATA['geology']['geology'].data[:,:,0], (1,0)), extent=self.GRID.extent, origin='lower', cmap='gray_r', alpha=0.5)

#         # Set up iteration structure:
#         iteration = 0

#         # Iterate in outlets groups
#         for outlet_iteration in range(len(self.outlets_importance)):

#             if self.SETTINGS['verbosity'] > 2:
#                 print('Total Iteration:', iteration, 'Outlet iteration:', outlet_iteration)

#             outlets_current = self._outlets[self._outlets['outlet_iteration'] == outlet_iteration]
#             if self.SETTINGS['verbosity'] > 1:
#                 ax.scatter(outlets_current.x, outlets_current.y, c='c') # Debugging

#             lst = outlets_current.index.values.tolist()
#             inlets_current = self._inlets[self._inlets['outlet_key'].isin(lst)]

#             # Iterate in inlets groups
#             for inlet_iteration in range(len(self.inlets_importance)):
#                 inlets_current_iteration = inlets_current[inlets_current['inlet_iteration'] == inlet_iteration]
#                 self._outlets.loc[outlets_current.index, 'iteration'] = iteration          # Store total iteration counter
#                 self._inlets.loc [inlets_current_iteration.index, 'iteration'] = iteration # Store total iteration counter

#                 #Compute travel time maps and conduit network
#                 if self.SETTINGS['algorithm'] == 'Isotropic3':
#                     self._compute_cost_map(iteration)
#                     self._compute_time_map_isotropic(iteration)
#                     self._compute_karst_map(iteration)
#                 elif self.SETTINGS['algorithm'] == 'Riemann3':
#                     self._compute_cost_map(iteration)
#                     self._compute_alpha_map(iteration)
#                     self._compute_beta_map(iteration)
#                     self._compute_riemann_metric(iteration)
#                     self._compute_time_map_riemann(iteration)
#                     self._compute_karst_map(iteration)
#                 else:
#                     print('Unrecognized algorithm:', self.SETTINGS['algorithm'])
#                 iteration = iteration + 1   #increment total iteration number by 1

#                 if self.SETTINGS['verbosity'] > 0:
#                     print('iteration:{}/{}'.format(iteration, self.nbr_iteration))
#         if self.SETTINGS['verbosity'] > 0:
#             print('- END -')

# #             for (o, outlet) in outlets_current.iterrows():                     #loop over outlets in current outlet iteration
# #
# #                 if self.SETTINGS['verbosity'] > 2:
# #                     print('\t Current outlet index:', outlet.name)             #debugging
# #
# #                 if self.SETTINGS['verbosity'] > 1:
# #                     ax.annotate(str(outlet_iteration), (outlet.x, outlet.y))
# #
# #                 inlets_outlet = self._inlets[self._inlets['outlet_key']==outlet.name]         #get the inlets assigned to current outlet
# #
# #                 if self.SETTINGS['verbosity'] > 1:
# #                     print('\t Inlets assigned to this outlet:\n', inlets_outlet)
# #
# #                 for inlet_iteration in range(len(self.inlets_importance)): #loop over inlet iterations
# #
# #                     if self.SETTINGS['verbosity'] > 2:
# #                         print('\t\t Inlet iteration:', inlet_iteration)
# #
# #                     inlets_current = inlets_outlet[inlets_outlet['inlet_iteration']==inlet_iteration] #get the inlets assigned to the current inlet iteration
# #
# #                     if self.SETTINGS['verbosity'] > 2:
# #                         print(inlets_current)                                                         #debugging
# #
# #                     if self.SETTINGS['verbosity'] > 1:
# #                         ax.scatter(inlets_current.x, inlets_current.y)
# #
# #                     for (i, inlet)in inlets_current.iterrows():                                 #loop over inlet in current inlet iteration
# #
# #                         if self.SETTINGS['verbosity'] > 1:
# #                             ax.annotate(str(outlet_iteration)+'-'+str(inlet_iteration), (inlet.x,inlet.y))  #debugging
# #
# #                         # self._outlets.loc[self._outlets.index == outlet.name, 'iteration'] = iteration   #store total iteration counter
# #                         # self._inlets.loc [self._inlets.index  == inlet.name , 'iteration'] = iteration   #store total iteration counter
# #                         self._outlets.loc[i, 'iteration'] = iteration   #store total iteration counter
# #                         self._inlets.loc [i, 'iteration'] = iteration   #store total iteration counter
# #
#         return None


#     # ii. Iso- and anisotropic case
#     def _compute_cost_map(self, iteration):
#         """
#         Compute the cost map (how difficult it is to traverse each cell).

#         TODO : à terminer
#         """
#         # If it's the first iteration, iniatialize the cost map according to the geological settings.
#         if iteration == 0:
#             # Geology
#             geology = self.DATA['geology']['geology']
#             if geology.mode == 'null':
#                 self.maps['cost'][0] = np.full((self.GRID.nx, self.GRID.ny, self.GRID.nz), self.SETTINGS['cost_aquifer']) #every cell has the same travel cost and is part of the aquifer
#             elif geology.mode == 'image':
#                 self.maps['cost'][0] = np.where(self.DATA['geology']['geology'].data == 1, self.SETTINGS['cost_aquiclude'], self.SETTINGS['cost_aquifer'])
#             elif geology.mode == 'gslib':
#                 for key, cost_name in geology.cost.items():
#                     try:
#                         self.maps['cost'][0] = np.where(self.DATA['geology']['geology'].data == key, self.SETTINGS[cost_name], self.maps['cost'][0])
#                     except:
#                         # TODO
#                         print('oups')

#             # Faults
#             if self.SETTINGS['geologic_features']['faults']['is_active']:
#                 self.maps['cost'][0] = np.where(self.DATA['geology']['faults'].data > 0, self.SETTINGS['cost_faults'], self.maps['cost'][0])

#             # Fractures
#             if self.SETTINGS['geologic_features']['fractures']['is_active']:
#                 self.maps['cost'][0] = np.where(self.DATA['geology']['fractures'].data > 0, self.SETTINGS['cost_fractures'], self.maps['cost'][0])

#             # If out of mask
#             if self.SETTINGS['mask']['is_active']:
#                 if self.DATA['mask'] is not None:
#                     self.maps['cost'][0] = np.where(self.DATA['mask'] == 1, self.SETTINGS['cost_out'], self.maps['cost'][0])

#         else: # If it's not the first iteration
#             self.maps['cost'][iteration] = self.maps['cost'][iteration-1] # set cost map to be same as previous iteration
#             self.maps['cost'][iteration] = np.where(self.maps['karst'][iteration-1] > 0, self.SETTINGS['cost_conduits'], self.maps['cost'][iteration]) # where karst conduits are present from previous iteration, set cost to conduit cost, elsewhere, leave unchanged
#         return None

#     # ii. Isotropic case
#     def _compute_time_map_isotropic(self, iteration):
#         """
#         Compute the travel time map (how long it takes to get to the outlet from each cell),
#         using the isotropic agd-hfm fast-marching algorithm, and store travel time map.
#         Note: the AGD-HFM library uses different indexing, so x and y indices are reversed for inlets and outlets.
#         """
#         # Set the outlets for this iteration
#         seeds_z = self._outlets[self._outlets['iteration']==iteration].z
#         seeds_y = self._outlets[self._outlets['iteration']==iteration].y
#         seeds_x = self._outlets[self._outlets['iteration']==iteration].x
#         seeds   = list(zip(seeds_z, seeds_y, seeds_x))
#         self.fastMarching['seeds'] = seeds

#         # Select inlets for current iteration
#         tips_z = self._inlets[self._inlets['iteration']==iteration].z
#         tips_y = self._inlets[self._inlets['iteration']==iteration].y
#         tips_x = self._inlets[self._inlets['iteration']==iteration].x
#         tips   = list(zip(tips_z, tips_y, tips_x))
#         self.fastMarching['tips'] = tips

#         # Set the isotropic travel cost through each cell
#         self.fastMarching['cost'] = np.transpose(self.maps['cost'][iteration], (2,1,0))

#         # Set verbosity of hfm run
#         self.fastMarching['verbosity'] = self.SETTINGS['verbosity']

#         # Run the fast marching algorithm and store the outputs
#         self.fastMarchingOutput = self.fastMarching.Run()

#         # Store travel time maps
#         self.maps['time'][iteration] = np.transpose(self.fastMarchingOutput['values'], (2,1,0))
#         return None


#     # ii. Anisotropic case
#     def _compute_alpha_map(self, iteration):
#         """
#         Compute the alpha map: travel cost in the same direction as the gradient.
#         Cost map * topography map, so that the cost is higher at higher elevations, encouraging conduits to go downgradient.

#         TODO : à terminer
#         """
#         self.maps['alpha'][iteration] = self.maps['cost'][iteration] * self.GRID.Z

#         # if self.settings['topography_mode'] != 'null':
#         #     self.maps['alpha'][iteration] = self.maps['cost'][iteration] * self.data['topography'].data
#         # elif self.settings['orientation_mode'] == 'surface':
#         #     self.maps['alpha'][iteration] = self.maps['cost'][iteration] * self.data['surface'].data
#         # else:
#         #     self.maps['alpha'][iteration] = self.maps['cost'][iteration]

#         return None

#     # ii. Anisotropic case
#     def _compute_beta_map(self, iteration):
#         """
#         Compute the beta map: travel cost perpendicular to the gradient.
#         If beta is higher than alpha, conduits will follow the steepest gradient.
#         If beta is lower than alpha, conduits will follow contours.
#         """
#         self.maps['beta'][iteration] = self.maps['alpha'][iteration] / self.SETTINGS['cost_ratio']
#         return None

#     # ii. Anisotropic case
#     def _compute_riemann_metric(self, iteration):
#         """
#         Compute the riemann metric: Define the Riemannian metric needed as input for the anisotropic fast marching.

#         TODO : à terminer
#         """
#         x, y, z = self.GRID.Z.shape
#         orientationx = np.transpose(np.gradient(self.GRID.Z, axis=0), (2,1,0))
#         orientationy = np.transpose(np.gradient(self.GRID.Z, axis=1), (2,1,0))
#         if (z == 1):
#             orientationz = np.transpose(self.GRID.Z, (2,1,0)) # TODO ??!!
#         else:
#             orientationz = np.transpose(np.gradient(self.GRID.Z, axis=2), (2,1,0))

#         alpha = np.transpose(self.maps['alpha'][iteration], (2,1,0))
#         beta  = np.transpose(self.maps['beta'][iteration],  (2,1,0))
#         self.riemannMetric = agd.Metrics.Riemann.needle([orientationz, orientationy, orientationx], alpha, beta)
#         return None


#     # ii. Anisotropic case
#     def _compute_time_map_riemann(self, iteration):
#         """
#         Compute the travel time map (how long it takes to get to the outlet from each cell),
#         using the anisotropic agd-hfm fast-marching algorithm, and store travel time map.
#         Note: the AGD-HFM library uses different indexing, so x and y indices are reversed for inlets and outlets.
#         """
#         # Set the outlets for this iteration
#         seeds_z = self._outlets[self._outlets['iteration']==iteration].z
#         seeds_y = self._outlets[self._outlets['iteration']==iteration].y
#         seeds_x = self._outlets[self._outlets['iteration']==iteration].x
#         seeds   = list(zip(seeds_z, seeds_y, seeds_x))
#         self.fastMarching['seeds'] = seeds

#         # Select inlets for current iteration
#         tips_z = self._inlets[self._inlets['iteration']==iteration].z
#         tips_y = self._inlets[self._inlets['iteration']==iteration].y
#         tips_x = self._inlets[self._inlets['iteration']==iteration].x
#         tips   = list(zip(tips_z, tips_y, tips_x))
#         self.fastMarching['tips'] = tips

#         # Set the travel cost through each cell
#         self.fastMarching['metric'] = self.riemannMetric

#         # Set verbosity of hfm run
#         self.fastMarching['verbosity'] = self.SETTINGS['verbosity']

#         # Run the fast marching algorithm and store the outputs
#         self.fastMarchingOutput = self.fastMarching.Run()

#         # Store travel time maps
#         self.maps['time'][iteration] = np.transpose(self.fastMarchingOutput['values'], (2,1,0))

#         # Store fastest travel paths
#         self.geodesics.append(self.fastMarchingOutput['geodesics'])
#         return None


#     # ii.
#     def _compute_karst_map(self, iteration):
#         """
#         Compute the karst map based on the paths from agd-hfm.
#         Array of all zeros, with ones in cells containing a karst conduit.
#         """
#         # Get karst map from previous iteration (except for the very first iteration)
#         if iteration > 0:
#             self.maps['karst'][iteration] = self.maps['karst'][iteration-1]

#         # Debugging plot:
#         # Chloe: this should stay in since it is very useful if there are problems
#         # f1, ax1 = plt.subplots(1, 1, figsize=(10, 10))

#         ### Loop over conduit paths generated by fast marching:
#         for path in self.fastMarchingOutput['geodesics']:   #loop over conduit paths in this iteration (there is one path from each inlet)
#             merge = False                                   #reset indicator for whether this conduit has merged with an existing conduit
#             for p in range(path.shape[1]):                  #loop over points making up this conduit path
#                 point = path[:,p]                           #get coordinates of current point
#                 # print('p', point)
#                 [[iz, iy, ix], error] = self.fastMarching.IndexFromPoint(point) #convert to coordinates to indices, /!\ returning iy first then ix
#                 # ax1.scatter(point[1], point[0], c='g',s=5)  #debugging

#                 ############## WHAT IS HAPPENING HERE?
#                 if ix < 0 or iy < 0 or iz < 0:
#                     print(ix,iy,iz)
#                     continue

#                 #Place nodes and links:
#                 if np.isnan(self.maps['nodes'][ix, iy, iz]):                                    #if there is no existing conduit node here
#                     if ~np.isnan(self.maps['outlets'][ix, iy, iz]):                              #if there is an outlet here (cell value is not nan)
#                         outlet = self._outlets.iloc[int(self.maps['outlets'][ix, iy, iz])]         #get the outlet coordinates using the ID in the outlets map
#                         self.nodes[self.n]             = [outlet.x, outlet.y, outlet.z, 'outfall']           #add a node at the outlet coordinates (with the node type for SWMM)
#                         self.maps['nodes'][ix, iy, iz] = self.n                                   #update node map with node index
#                         # ax1.scatter(outlet.x, outlet.y, marker='o', c='b')                   #debugging
#                         if p > 0:                                                           #if this is not the first point (i.e. the inlet) in the current path
#                             if merge == False:                                               #if this conduit has not merged with an existing conduit
#                                 self.edges[self.e] = [self.n-1, self.n]                       #add an edge connecting the previous node to the current node
#                                 self.e = self.e + 1                                             #increment edge counter up by one
#                                 # ax1.plot((self.nodes[self.n][0], self.nodes[self.n-1][0]), (self.nodes[self.n][1], self.nodes[self.n-1][1]))
#                             else:                                                          #if this conduit HAS merged with an existing conduit
#                                 [[fromiz, fromiy, fromix], error] = self.fastMarching.IndexFromPoint(path[:, p-1]) #get xyz indices of previous point in current conduit path
#                                 n_from = self.maps['nodes'][fromix, fromiy, fromiz]           #get node index of the node already in the cell where the previous point was
#                                 self.edges[self.e] = [n_from, self.n]                         #add an edge connecting existing conduit node to current node
#                                 self.e = self.e + 1                                             #increment edge counter up by one
#                                 # ax1.plot((self.nodes[self.n].x, self.nodes[n_from].x), (self.nodes[self.n].y, self.nodes[n_from].y))
#                         self.n = self.n + 1                                                   #increment node counter up by one
#                     else:                                                                  #if there is NOT an outlet here
#                         if p > 0:                                                           #if this is not the first point in the current path
#                             #possible improvement: if the next point on the path is on an existing point, skip the current point.
#                             self.nodes[self.n] = [point[2], point[1], point[0], 'junction']            #add a junction node here (with the node type for SWMM)
#                             self.maps['nodes'][ix, iy, iz] = self.n                               #update node map with node index
#                             #ax1.scatter(point[1],point[0], marker='o', edgecolor='g', facecolor='none')   #debugging
#                             if merge == False:                                              #if this conduit has not merged with an existing conduit
#                                 self.edges[self.e] = [self.n-1, self.n]                      #add and edge connecting the previous node to the current node
#                                 self.e = self.e + 1                                            #increment edge counter up by one
#                                 #ax1.plot((self.nodes[self.n][1], self.nodes[self.n-1][1]),(self.nodes[self.n][0], self.nodes[self.n-1][0]), c='gold', marker=None)
#                             else:                                                           #if this conduit HAS merged with an existing conduit
#                                 [[fromiz, fromiy, fromix], error] = self.fastMarching.IndexFromPoint(path[:, p-1]) #get xy indices of previous point in current conduit path
#                                 n_from = self.maps['nodes'][fromix, fromiy, fromiz]                   #get node index of the node already in the cell where the previous point was
#                                 self.edges[self.e] = [n_from, self.n]                        #add an edge connecting existing conduit node to current node
#                                 self.e = self.e + 1                                            #increment edge counter up by one
#                                 merge = False                                                #reset merge indicator to show that current conduit has left                                                              #if this is the first point in current path
#                         else:                                                                #if this is the first point in the current path (counter <= 0, therefore it is an inlet)
#                             self.nodes[self.n] = [point[2], point[1], point[0], 'inlet']               #add an inlet node here (with the node type for SWMM)
#                             self.maps['nodes'][ix, iy, iz] = self.n                               #update node map with node index
#                             #ax1.scatter(point[1],point[0], marker='o', edgecolor='sienna', facecolor='none')
#                         self.n = self.n + 1                                                   #increment node counter up by one
#                 elif ~np.isnan(self.maps['nodes'][ix, iy, iz]):                                 #if there is already a node in this cell (either because there is a conduit here, or because there are two nodes in the same cell)
#                     n_existing = self.maps['nodes'][ix, iy, iz]                                  #get index of node already present in current cell
#                     if merge == True:                                                       #if this conduit has already merged into an existing conduit
#                         pass                                                                 #skip this node (there is already a node here)
#                     elif n_existing == self.n-1:                                            #if existing index is only one less than next node to be added index, this is a duplicate node and can be skipped
#                         pass                                                                 #skip this node (duplicate)
#                     else:                                                                   #if existing node index is >1 less than next node to be added index
#                         if p > 0:                                                           #if this is not the first point in the current path
#                             self.edges[self.e] = [self.n-1, n_existing]                      #add an edge connecting most recently added node and existing node in cell
#                             self.e = self.e + 1                                                #increment edge counter up by one
#                             merge = True                                                     #add a flag indicating that this conduit has merged into an existing one
#                             #ax1.plot((self.nodes[self.n-1][1], self.nodes[n_existing][1]),(self.nodes[self.n-1][0], self.nodes[n_existing][0]), c='r', marker=None)
#                         else:                                                                #if this is the first point in the current path (i.e. the inlet is on an exising conduit)
#                             self.nodes[self.n] = [point[2], point[1], point[0], 'inlet']                #add a node here (with the node type for SWMM)- this will cause there to be two nodes in the same cell
#                             self.maps['nodes'][ix, iy, iz] = self.n                                #update node map with node index
#                             self.n = self.n + 1                                                 #increment node counter by 1
#                             #ax1.scatter(point[1],point[0], marker='o', edgecolor='g', facecolor='none')  #debugging
#                 self.maps['karst'][iteration][ix, iy, iz] = 1                               #update karst map to put a conduit in current cell


#         ### Debugging plot
#         # # Display inlets and outlets
#         # ax1.scatter(self._outlets.x, self._outlets.y, c='cyan',   s=100)
#         # ax1.scatter(self._inlets.x,  self._inlets.y,  c='orange', s=100)
#         # ax1.scatter(self._outlets[self._outlets.iteration==iteration].x, self._outlets[self._outlets.iteration==iteration].y, c='cyan',   s=100)
#         # ax1.scatter(self._inlets[self._inlets.iteration==iteration].x,   self._inlets[self._inlets.iteration==iteration].y,   c='orange', s=100)
#         # # Display karst network
#         # ax1.imshow(np.transpose(self.maps['karst'][iteration], (1,0,2)), origin='lower', extent=self.GRID.extent, cmap='gray_r')
#         # ax1.imshow(np.transpose(self.maps['nodes'], (1,0,2)), origin='lower', extent=self.GRID.extent, cmap='gray_r')
#         return None

#     ######################
#     ### Karst Analysis ###
#     ######################

#     def compare_stats(self, mean=False):
#         """
#         Compare statistics between reference indicators and calculated networks.

#         TODO
#         param 'iteration=0'

#         """
#         indicators = ['cpd', 'cv degree', 'cv length', 'orientation entropy', 'length entropy', 'aspl', 'mean degree', 'mean length', 'correlation vertex degree']
#         stats = pd.DataFrame(columns=indicators)

#         for (i, karst_network) in enumerate(self.SIMULATIONS):
#             stats.loc[i] = karst_network.stats

#         # Apply style
#         def _bg_color(x, min_val, max_val):
#             if (x < min_val) or (x > max_val):
#                 return 'background-color: red'
#             else:
#                 return 'background-color: green'

#         display(stats.style.applymap(_bg_color, min_val = self.reference_statistics['cpd']['min'], max_val = self.reference_statistics['cpd']['max'], subset = ['cpd'])\
#         .applymap(_bg_color, min_val = self.reference_statistics['cv degree']['min'], max_val = self.reference_statistics['cv degree']['max'], subset = ['cv degree'])\
#         .applymap(_bg_color, min_val = self.reference_statistics['cv length']['min'], max_val = self.reference_statistics['cv length']['max'], subset = ['cv length'])\
#         .applymap(_bg_color, min_val = self.reference_statistics['orientation entropy']['min'], max_val = self.reference_statistics['orientation entropy']['max'], subset = ['orientation entropy'])\
#         .applymap(_bg_color, min_val = self.reference_statistics['length entropy']['min'], max_val = self.reference_statistics['length entropy']['max'], subset = ['length entropy'])\
#         .applymap(_bg_color, min_val = self.reference_statistics['aspl']['min'], max_val = self.reference_statistics['aspl']['max'], subset = ['aspl'])\
#         .applymap(_bg_color, min_val = self.reference_statistics['mean degree']['min'], max_val = self.reference_statistics['mean degree']['max'], subset = ['mean degree'])\
#         .applymap(_bg_color, min_val = self.reference_statistics['mean length']['min'], max_val = self.reference_statistics['mean length']['max'], subset = ['mean length'])\
#         .applymap(_bg_color, min_val = self.reference_statistics['correlation vertex degree']['min'], max_val = self.reference_statistics['correlation vertex degree']['max'], subset = ['correlation vertex degree']))
#         return None


#     def compute_average_paths(self, mask=0):
#         """
#         TODO
#         Compute the mean of all the simulations.
#         """

#         # Calculate the average from all the simulations
#         karst_maps = []
#         for karst_simulation in self.SIMULATIONS:
#             data = karst_simulation.maps['karst'][-1]
#             karst_maps.append(data)
#         karst_prob = sum(karst_maps)/len(karst_maps)

#         self.karst_prob = karst_prob
#         return karst_prob

#     def show_average_paths(self):
#         """
#         todo
#         """
#         ### Call the plotter
#         p = pv.Plotter(notebook=False)

#         ### Construct the grid
#         vtk = pv.UniformGrid()
#         vtk.dimensions = np.array((self.GRID.nx, self.GRID.ny, self.GRID.nz)) + 1
#         vtk.origin     = (self.GRID.x0 - self.GRID.dx/2, self.GRID.y0 - self.GRID.dy/2, self.GRID.z0 - self.GRID.dz/2)
#         vtk.spacing    = (self.GRID.dx, self.GRID.dy, self.GRID.dz)

#         vtk['values'] = self.karst_prob.flatten(order="F")

#         mesh = vtk.cast_to_unstructured_grid()
#         ghosts = np.argwhere(vtk['values'] < 1.0)
#         mesh.remove_cells(ghosts)
#         p.add_mesh(mesh, show_edges=False)

#         ### Plotting
#         # p.add_title(feature)
#         p.add_axes()
#         bounds = p.show_bounds(mesh=vtk)
#         p.add_actor(bounds)
#         p.show(cpos='xy')

#         return None




#     ###################
#     ### DEBUG PLOTS ###
#     ###################

#     def debug_plot(self, feature):
#         """
#         TODO
#         """
#         ### Call the plotter
#         p = pv.Plotter(notebook=False)

#         ### Construct the grid
#         vtk = pv.UniformGrid()
#         vtk.dimensions = np.array((self.GRID.nx, self.GRID.ny, self.GRID.nz)) + 1
#         vtk.origin     = (self.GRID.x0 - self.GRID.dx/2, self.GRID.y0 - self.GRID.dy/2, self.GRID.z0 - self.GRID.dz/2)
#         vtk.spacing    = (self.GRID.dx, self.GRID.dy, self.GRID.dz)

#         ### Parameters according to feature
#         if feature == 'grid':
#             pass
#         if feature == 'mask':
#             vtk['values'] = self.mask.mask.flatten(order="F")
#         if (feature == 'geology') or (feature == 'faults') or (feature == 'fractures') or (feature == 'field'):
#             vtk['values'] = self.geology.geologic_features[feature].data.flatten(order="F")
#         if feature == 'points':
#             try:
#                 vtk['values'] = self.mask.mask.flatten(order="F")
#             except:
#                 try:
#                     vtk['values'] = self.geology.geologic_features[feature].data.flatten(order="F")
#                 except:
#                     pass
#             inlets  = self.points.point_features['inlets'] .points[['x', 'y', 'z']].values
#             outlets = self.points.point_features['outlets'].points[['x', 'y', 'z']].values
#             cloud_i = pv.wrap(inlets)
#             cloud_o = pv.wrap(outlets)
#             p.add_points(cloud_i, render_points_as_spheres=True, point_size=10, color='r')
#             p.add_points(cloud_o, render_points_as_spheres=True, point_size=10, color='b')

#         p.add_mesh(vtk, show_edges=False)

#         ### Plotting
#         p.add_title(feature)
#         p.add_axes()
#         bounds = p.show_bounds(mesh=vtk)
#         p.add_actor(bounds)
#         p.show(cpos='xy')

#         # TODO : Creates a multiplot

#         # p = pv.Plotter(notebook=False)
#         # p.add_title(feature)
#         # slices = vtk.slice_orthogonal()
#         # p.add_mesh(slices, show_edges=False)
#         # p.show(cpos='xy')
#         return None

#     def debug_plot_initialize(self):
#         """
#         TODO
#         """
#         ### Call the plotter
#         p = pv.Plotter(shape=(2, 4), border=True, notebook=False)

#         ### Construct the grid
#         vtk = pv.UniformGrid()
#         vtk.dimensions = np.array((self.GRID.nx, self.GRID.ny, self.GRID.nz)) + 1
#         vtk.origin     = (self.GRID.x0 - self.GRID.dx/2, self.GRID.y0 - self.GRID.dy/2, self.GRID.z0 - self.GRID.dz/2)
#         vtk.spacing    = (self.GRID.dx, self.GRID.dy, self.GRID.dz)

#         features = ['geology', 'faults', 'fractures', 'field']

#         for (i, feature) in enumerate(features):
#             vtk = vtk.copy()
#             p.subplot(0, i)
#             p.add_text(feature, font_size=24)
#             try:
#                 vtk['values'] = self.geology.geologic_features[feature].data.flatten(order="F")
#             except:
#                 vtk['values'] = np.full((self.GRID.nx, self.GRID.ny, self.GRID.nz), 0).flatten(order="F")
#             p.add_mesh(vtk, show_edges=False)
#             bounds = p.show_bounds(mesh=vtk)
#             p.add_actor(bounds)

#         for (i, feature) in enumerate(features):
#             vtk = vtk.copy()
#             p.subplot(1, i)
#             p.add_text(feature, font_size=24)
#             try:
#                 vtk['values'] = self.geology.geologic_features[feature].data.flatten(order="F")
#             except:
#                 vtk['values'] = np.full((self.GRID.nx, self.GRID.ny, self.GRID.nz), 0).flatten(order="F")
#             slices = vtk.slice_orthogonal()
#             p.add_mesh(slices, show_edges=False)
#             bounds = p.show_bounds(mesh=vtk)
#             p.add_actor(bounds)

#         p.link_views()
#         p.show(cpos='xy')
#         return None

#     def debug_plot_compute(self, iteration=-1):
#         """
#         TODO
#         """
#         ### Call the plotter
#         p = pv.Plotter(shape=(1, 3), border=True, notebook=False)

#         ### Construct the grid
#         vtk = pv.UniformGrid()
#         vtk.dimensions = np.array((self.GRID.nx, self.GRID.ny, self.GRID.nz)) + 1
#         vtk.origin     = (self.GRID.x0 - self.GRID.dx/2, self.GRID.y0 - self.GRID.dy/2, self.GRID.z0 - self.GRID.dz/2)
#         vtk.spacing    = (self.GRID.dx, self.GRID.dy, self.GRID.dz)

#         features = ['cost', 'time', 'karst']

#         for (i, feature) in enumerate(features):
#             vtk = vtk.copy()
#             p.subplot(0, i)
#             p.add_text(feature, font_size=24)
#             vtk['values'] = self.maps[feature][iteration].flatten(order="F")
#             try:
#                 vtk['values'] = self.maps[feature][iteration].flatten(order="F")
#             except:
#                 vtk['values'] = np.full((self.GRID.nx, self.GRID.ny, self.GRID.nz), 0).flatten(order="F")

#             bounds = p.show_bounds(mesh=vtk)
#             p.add_actor(bounds)

#             if (feature == 'cost'):
#                 mesh = vtk.cast_to_unstructured_grid()
#                 ghosts = np.argwhere(vtk['values'] > 1.0)
#                 mesh.remove_cells(ghosts)
#                 p.add_mesh(mesh, show_edges=False)
#                 # p.add_mesh(vtk, show_edges=False, clim=[0,1])

#             if (feature == 'time'):
#                 p.add_mesh(vtk, show_edges=False)

#             if (feature == 'karst'):
#                 mesh = vtk.cast_to_unstructured_grid()
#                 ghosts = np.argwhere(vtk['values'] < 1.0)
#                 mesh.remove_cells(ghosts)
#                 p.add_mesh(mesh, show_edges=False)

#                 inlets  = self.points.point_features['inlets'] .points[['x', 'y', 'z']].values
#                 outlets = self.points.point_features['outlets'].points[['x', 'y', 'z']].values
#                 cloud_i = pv.wrap(inlets)
#                 cloud_o = pv.wrap(outlets)
#                 p.add_points(cloud_i, render_points_as_spheres=True, point_size=10, color='r')
#                 p.add_points(cloud_o, render_points_as_spheres=True, point_size=10, color='b')

#             # p.subplot(1, i)
#             # p.add_text(feature, font_size=24)
#             # slices = vtk.slice_orthogonal()
#             # bounds = p.show_bounds(mesh=vtk)
#             # p.add_actor(bounds)

#         p.link_views()
#         p.show(cpos='xy')
#         return None

#     def debug_plot_fmm_feature(self, feature, iteration=-1):
#         """
#         TODO
#         """
#         ### Construct the grid
#         vtk = pv.UniformGrid()
#         vtk.dimensions = np.array((self.GRID.nx, self.GRID.ny, self.GRID.nz)) + 1
#         vtk.origin     = (self.GRID.x0 - self.GRID.dx/2, self.GRID.y0 - self.GRID.dy/2, self.GRID.z0 - self.GRID.dz/2)
#         vtk.spacing    = (self.GRID.dx, self.GRID.dy, self.GRID.dz)

#         ### Parameters according to feature
#         if (feature == 'costmap') or (feature == 'timemap'):
#             p = pv.Plotter(shape=(1, 2), border=True, notebook=False)
#             vtk['values'] = self.maps[feature[:4]][iteration].flatten(order="F")

#         if feature == 'karstmap':
#             p = pv.Plotter(notebook=False)
#             vtk['values'] = self.maps['karst'][iteration].flatten(order="F")
#             mesh = vtk.cast_to_unstructured_grid()
#             ghosts = np.argwhere(vtk['values'] < 1.0)
#             mesh.remove_cells(ghosts)
#             vtk = mesh

#             inlets  = self.points.point_features['inlets'] .points[['x', 'y', 'z']].values
#             outlets = self.points.point_features['outlets'].points[['x', 'y', 'z']].values
#             cloud_i = pv.wrap(inlets)
#             cloud_o = pv.wrap(outlets)
#             p.add_points(cloud_i, render_points_as_spheres=True, point_size=10, color='r')
#             p.add_points(cloud_o, render_points_as_spheres=True, point_size=10, color='b')

#         p.add_mesh(vtk, show_edges=False)

#         if feature != 'karstmap':
#             p.subplot(0, 1)
#             # slices = vtk.slice_orthogonal()
#             slices = vtk.slice_along_axis(n=5, axis="y")
#             p.add_mesh(slices, show_edges=False)

#             p.link_views()

#         ### Plotting
#         # p.add_title(feature + '/ iteration:' + str(iteration))
#         # p.add_axes()
#         # bounds = p.show_bounds(mesh=vtk)
#         # p.add_actor(bounds)


#         p.show(cpos='xy')

#         return None

# #############################
# ### Visualization methods ###
# #############################

#     # def show():
#     #     fig, ax = plt.subplots(figsize=(10, 10))
#     #     ax.imshow(np.transpose(, (1,0,2)), origin="lower", extent=self.GRID.extent)

#     # def show(self):
#     #     """
#     #     Displays the apparent last computed karst network in 3 directions.
#     #     """
#     #
#     #     karst = self.maps['karst'][-1]
#     #
#     #     karst_xy = np.sum(karst, axis=2)
#     #     karst_zx = np.sum(karst, axis=0)
#     #     karst_zy = np.sum(karst, axis=1)
#     #
#     #     fig, ax = plt.subplots(1, 3, figsize=(20, 10))



#     # def show(self):
#     #     """
#     #     Displays the last computed karst network.
#     #     """
#     #     import pyvista as pv
#     #     vtk = pv.UniformGrid()
#     #     vtk.dimensions = np.array((self.GRID.nx, self.GRID.ny, self.GRID.nz)) + 1
#     #     vtk.origin     = (self.GRID.x0 - self.GRID.dx/2, self.GRID.y0 - self.GRID.dy/2, self.GRID.z0 - self.GRID.dz/2)
#     #     vtk.spacing    = (self.GRID.dx, self.GRID.dy, self.GRID.dz)
#     #     vtk.values     = self.maps['karst'][-1].flatten(order="F")
#     #     vtk.
#     #     # pv.set_jupyter_backend(option)
#     #     vtk.plot(show_edges=True)
#     #     return None





# # #     #############################
# # #     ### Visualization methods ###
# # #     #############################
# # #
# # #     def show_catchment(self, label='geology', title=None, cmap='binary'):
# # #         """
# # #         Show the entire study domain.
# # #
# # #         Parameters
# # #         ----------
# # #         label : str, optional
# # #             Data to show : 'geology', 'topography', 'orientationx', 'orientationy' 'faults' or 'fractures'.
# # #             By default : 'geology'.
# # #         title : str, optional
# # #             Title of the plot. If 'None', 'data' becomes the label.
# # #         cmap : str, optional
# # #             Color map, 'binary' by default.
# # #         """
# # #         import matplotlib.patches as mp
# # #         fig, ax1 = plt.subplots()
# # #         #if title is None:
# # #         #    title = label
# # #         #fig.suptitle(title, fontsize=16)
# # #
# # #         # Load data to show
# # #         try:
# # #             data = [data.data[:,:,0] for data in self.geology.data if data.label==label][-1]
# # #         except:
# # #             print('no data for indicated label parameter')
# # #
# # #         im1 = ax1.imshow(data.T, origin="lower", extent=self.grid.extent, cmap=cmap)
# # #
# # #         fig.colorbar(im1, ax=ax1)
# # #         if self.settings['data_has_mask']:
# # #             import matplotlib.patches as mp
# # #             p1 = mp.PathPatch(self.mask.polygon, lw=2, fill=0, edgecolor='red', label='mask')
# # #             ax1.add_patch(p1)
# # #
# # #         # Points
# # #         for pts in self.points.points:
# # #             x, y = zip(*pts.points)
# # #             ax1.plot(x, y, 'o', label=pts.points_key)
# # #
# # #         ax1.set_aspect('equal', 'box')
# # #         plt.legend(loc='upper right')
# # #         plt.show()
# # #         return fig
# # #
# # #     def _show_maps(self, sim=-1, iteration=-1, cmap='binary'):
# # #         """
# # #         Show the simulated karst network as an image.
# # #         """
# # #         karst_network = self.karst_simulations[sim]
# # #
# # #         fig, ([ax1,ax2],[ax3,ax4]) = plt.subplots(2, 2, sharex=True, sharey=True)
# # #         fig.suptitle('Karst Network', fontsize=16)
# # #
# # #         ax1.imshow(karst_network.maps['outlets'], extent=self.grid.extent, origin='lower', cmap=cmap)
# # #         ax1.set_title('Outlets')
# # #
# # #         ax2.imshow(karst_network.maps['cost'][iteration], extent=self.grid.extent, origin='lower', cmap=cmap)
# # #         ax2.set_title('Cost')
# # #
# # #         ax3.imshow(karst_network.maps['time'][iteration], extent=self.grid.extent, origin='lower', cmap=cmap)
# # #         ax3.set_title('Time')
# # #
# # #         ax4.imshow(karst_network.maps['karst'][iteration], extent=self.grid.extent, origin='lower', cmap=cmap)
# # #         ax4.set_title('Karst')
# # #
# # #         fig.subplots_adjust(hspace=0.5)
# # #         plt.show()
# # #         return None
# # #
# # #
# # #     def show(self, data=None, title=None):
# # #         """
# # #         Show the entire study domain (defaults to showing most recent simulation).
# # #         """
# # #         if data is None:
# # #             data = self.karst_simulations[-1]
# # #
# # #         fig = plt.figure(figsize=(20,10))
# # #
# # #         # Cost map
# # #         fig.add_subplot(131, aspect='equal')
# # #         d = data.maps['cost'][-1]
# # #         plt.xlabel('Cost array'+str(d.shape))
# # #         d = np.transpose(d, (1,0)) # imshow read MxN and we have NxM
# # #         plt.imshow(d, extent=self.grid.extent, origin='lower', cmap='gray') #darker=slower
# # #         plt.colorbar(fraction=0.046, pad=0.04)
# # #
# # #         # Travel time map
# # #         fig.add_subplot(132, aspect='equal')
# # #         d = data.maps['time'][-1]
# # #         plt.xlabel('Travel time array'+str(d.shape))
# # #         d = np.transpose(d, (1,0)) # imshow read MxN and we have NxM
# # #         plt.imshow(d, extent=self.grid.extent, origin='lower', cmap='cividis') #darker=faster
# # #         plt.colorbar(fraction=0.046, pad=0.04)
# # #
# # #         # Karst map
# # #         fig.add_subplot(133, aspect='equal')
# # #         d = data.maps['karst'][-1]
# # #         plt.xlabel('Karst array'+str(d.shape))
# # #         d = np.transpose(d, (1,0)) # imshow read MxN and we have NxM
# # #         plt.imshow(d, extent=self.grid.extent, origin='lower', cmap='gray_r') #darker=conduits
# # #         plt.colorbar(fraction=0.046, pad=0.04)
# # #         i = plt.scatter(data.points['inlets'].x,  data.points['inlets'].y,  c='orange')
# # #         o = plt.scatter(data.points['outlets'].x, data.points['outlets'].y, c='steelblue')
# # #         p = matplotlib.patches.Rectangle((0,0),0,0, ec='r', fc='none')
# # #         if self.settings['data_has_mask']:
# # #             closed_polygon = self.mask.vertices[:]
# # #             closed_polygon.append(closed_polygon[0])
# # #             x,y = zip(*closed_polygon)
# # #             plt.plot(x,y, color='red', label='mask')
# # #         #plt.legend([i,o,p], ['inlets', 'outlets', 'catchment'], loc='upper right')
# # #         plt.legend([i,o], ['inlets', 'outlets'], loc='upper right')
# # #
# # #         if title is not None:
# # #             fig.suptitle(title, fontsize=16)
# # #         plt.show()
# # #         return fig
# # #
# # #     def show_network(self, data=None, simplify=False, ax=None, plot_nodes=True, mask=True, labels=['inlets', 'outlets'], title=None, cmap=None, color='k', alpha=1, legend=True):
# # #         """
# # #         #Chloe: This is a new function that I use to create all the figures for the paper.
# # #         Show the karst network as a graph with nodes and edges. Defaults to showing latest iteration.
# # #
# # #         Parameters
# # #         ----------
# # #         data:
# # #             karst simulation object containing nodes, edges, points, etc. Can be obtained from self.karst_simulations[i]
# # #         ax :
# # #             axis to plot on
# # #         label :
# # #             None or list of strings ['nodes','edges','inlets','outlets'], indicating which components to label
# # #         title : str
# # #             title of plot
# # #         cmap : str
# # #             colormap to use when plotting
# # #         color : str
# # #             single color to use when plotting (cannot have both cmap and color)
# # #         alpha : float
# # #             opacity to plot with (1=opaque, 0=transparent)
# # #         legend : bool
# # #             whether to display legend
# # #         plot_nodes : bool
# # #             whether to display nodes
# # #         polygon : bool
# # #             whether to display the bounding polygon
# # #         """
# # #
# # #         if ax == None:
# # #             fig,ax = plt.subplots(figsize=(10,10))
# # #             ax.set_aspect('equal')
# # #
# # #         if data == None:
# # #             data = self.karst_simulations[-1]
# # #
# # #         if mask == True:
# # #             if self.settings['data_has_mask']:
# # #                 closed_polygon = self.mask.vertices[:]
# # #                 closed_polygon.append(closed_polygon[0])
# # #                 x,y = zip(*closed_polygon)
# # #                 ax.plot(x,y, color='maroon')
# # #                 p = matplotlib.lines.Line2D([0],[0], color='k')
# # #
# # #         if simplify == True:
# # #             nodes = data.network['nodes']   #get all nodes
# # #             nodes_simple = data.network['karstnet'].graph_simpl.nodes  #get indices of only the nodes in the simplified graph
# # #             nodes_simple = {key: nodes[key] for key in nodes_simple}   #make df of only the nodes in the simplified graph, for plotting
# # #             edges = data.network['edges']   #get all edges
# # #             edges_simple = data.network['karstnet'].graph_simpl.edges  #get only the edges in the simplified graph
# # #             edges_simple = {i: edge for i,edge in enumerate(edges_simple)}   #make df of only the edges in the simplified graph, for p
# # #             nodes = pd.DataFrame.from_dict(nodes_simple, orient='index', columns=['x','y','type']) #convert to pandas for easier plotting
# # #             edges = pd.DataFrame.from_dict(edges_simple, orient='index', columns=['inNode','outNode'])
# # #         else:
# # #             nodes = pd.DataFrame.from_dict(data.network['nodes'], orient='index', columns=['x','y','type']) #convert to pandas for easier plotting
# # #             edges = pd.DataFrame.from_dict(data.network['edges'], orient='index', columns=['inNode','outNode'])
# # #
# # #         #Set up data for plotting:
# # #         fromX = nodes.x.loc[edges.inNode]      #calculate coordinates for link start and end points
# # #         fromY = nodes.y.loc[edges.inNode]
# # #         toX   = nodes.x.loc[edges.outNode]
# # #         toY   = nodes.y.loc[edges.outNode]
# # #
# # #         #Plot nodes and edges:
# # #         if plot_nodes:
# # #             n = ax.scatter(nodes.x,              nodes.y,                  c='k',         alpha=alpha, s=5)  #scatterplot nodes
# # #         i = ax.scatter(data.points['inlets'].x,  data.points['inlets'].y,  c='orange',    s=30) #scatterplot inlets
# # #         o = ax.scatter(data.points['outlets'].x, data.points['outlets'].y, c='steelblue', s=30) #scatterplot outlets
# # #         e = matplotlib.lines.Line2D([0],[0])                                                  #line artist for legend
# # #         for ind in edges.index:                                                               #loop over edge indices
# # #             if cmap is not None:
# # #                 ax.plot((fromX.iloc[ind], toX.iloc[ind]), (fromY.iloc[ind], toY.iloc[ind]), c=plt.cm.get_cmap(cmap)(ind/len(edges)), alpha=alpha)  #plot each edge, moving along color gradient to show order
# # #             elif color is not None:
# # #                 ax.plot((fromX.iloc[ind], toX.iloc[ind]), (fromY.iloc[ind], toY.iloc[ind]), c=color, alpha=alpha)  #plot each edge in same color
# # #
# # #         #Add labels:
# # #         if labels == None:
# # #             pass
# # #         else:
# # #             if 'nodes' in labels:                                         #label node indices
# # #                 for ind in nodes.index:                                   #loop over node indices
# # #                     ax.annotate(str(ind), xy=(nodes.y[ind]-10, nodes.x[ind]))  #annotate slightly to left of each node
# # #             if 'edges' in labels:
# # #                 for ind in edges.index:
# # #                     ax.annotate(str(ind), xy=(edges.y[ind]-10, edges.x[ind]))  #annotate slightly to left of each edge
# # #             if 'inlets' in labels:
# # #                 for index,inlet in data.points['inlets'].iterrows():
# # #                     ax.annotate(str(int(inlet.outlet))+'-'+str(int(inlet.inlet_iteration)),  xy=(inlet.x-(6*self.grid.dx),  inlet.y))
# # #             if 'outlets' in labels:
# # #                 for index,outlet in data.points['outlets'].iterrows():
# # #                     ax.annotate(str(int(outlet.name)), xy=(outlet.x-(4*self.grid.dx), outlet.y))
# # #
# # #         #Add legend & title:
# # #         if legend:
# # #             if plot_nodes:
# # #                 if plot_polygon:
# # #                     ax.legend([i,o,n,e,p],['inlets','outlets','nodes','edges','mask'])
# # #                 else:
# # #                     ax.legend([i,o,n,e],['inlets','outlets','nodes','edges'])
# # #             else:
# # #                 if plot_polygon:
# # #                     ax.legend([i,o,e,p],['inlets','outlets','edges','mask'])
# # #                 else:
# # #                     ax.legend([i,o,e],['inlets','outlets','edges','mask'])
# # #         if title is not None:
# # #             ax.set_title(title, fontsize=16)
# # #
# # #         return None
# # #
# # #

# # #





# ################################################################################################
# # #
# # #     def get_fractures_numbers(self):
# # #         """
# # #         TODO
# # #         Gets the number of fractures.
# # #
# # #         Returns
# # #         -------
# # #         result : dict
# # #             A dictionary with the number of fractures given by user, and with the number of fractures calculated by the model.
# # #
# # #         Examples
# # #         --------
# # #         >>> frac_nbr = catchment.get_fractures_numbers()
# # #         """
# # #         user = self.fractures_numbers
# # #         model = []
# # #         for frac_family in self.geology.fractures:
# # #             model.append(self.geology.fractures[frac_family]['frac_nbr'])
# # #         fracs = {'user' : user, 'model' : model}
# # #         return fracs
# # #
# # #     def increment_rand_seed(self):
# # #         """
# # #         Increment by one the value of the random seed.
# # #         """
# # #         self.settings['rand_seed'] += 1
# # #         np.random.seed(self.settings['rand_seed'])
# # #         return None




# # TO DELETE ?

# # #      ####################
# # #     ### UPDATE methods ###
# # #      ####################
# # #
# # #     def update_feature(self, feature):
# # #         """
# # #         TODO
# # #         Update the feature according to the value of the parameter.
# # #
# # #         Parameters
# # #         ----------
# # #         feature : str
# # #             'mask', 'inlets', 'outlets', 'geology', 'topography', 'orientation', 'faults', 'fractures', 'all'.
# # #
# # #         Examples
# # #         --------
# # #         >>> catchment.update_feature('geology')
# # #         """
# # #         assert feature in self.UPDATE_FUNCTIONS
# # #         FUNC = self.UPDATE_FUNCTIONS[feature]
# # #         FUNC(self, feature)
# # #         return None
# # #
# # #     def _update_mask(self, feature):
# # #         """
# # #         Updates the mask settings.
# # #         Applies to the model the new parameters given with the 'set_' methods.
# # #         """
# # #         if self.settings['data_has_mask']:
# # #             self.hasMask = True
# # #             self.points.hasMask = True
# # #             self.mask.set_mask(self.settings['mask_data'])
# # #             self.mask_array = self.mask.mask
# # #
# # #             self.points.inspect_points()
# # #             self.inlets  = [pts.points for pts in self.points.points if pts.points_key == 'inlets'][0]
# # #             self.outlets = [pts.points for pts in self.points.points if pts.points_key == 'outlets'][0]
# # #         else:
# # #             self.hasMask = False
# # #             self.points.hasMask = False
# # #
# # #             self.points.inspect_points()
# # #             self.inlets  = [pts.points for pts in self.points.points if pts.points_key == 'inlets'][0]
# # #             self.outlets = [pts.points for pts in self.points.points if pts.points_key == 'outlets'][0]
# # #         return None
# # #
# # #     def _update_points(self, feature): # TODO update a set of points based on ID
# # #         """
# # #         Update inlets or outlets settings.
# # #         """
# # #         points = [pts for pts in self.points.points if pts.points_key==feature][-1]
# # #         n_points = self.points._set_points(feature, self.settings[feature+'_mode'], self.settings[feature+'_data'], self.settings[feature+'_number'], func="_update_points()")
# # #         points.points = n_points
# # #         points.mode   = self.settings[feature+'_mode']
# # #         points.nbr    = len(n_points)
# # #         self.points.inspect_points()
# # #
# # #         if feature == 'inlets':
# # #             self.inlets = points.points
# # #             if self.settings['inlets_shuffle']:
# # #                 np.random.shuffle(self.inlets)
# # #         if feature == 'outlets':
# # #             self.outlets = points.points
# # #             if self.settings['outlets_shuffle']:
# # #                 np.random.shuffle(self.outlets)
# # #         return None
# # #
# # #     def _update_data(self, feature):
# # #         """
# # #         Update the geologic feature settings.
# # #         """
# # #         if feature == "orientation":
# # #             if self.settings['orientation_mode'] == "null":
# # #                 n_orientationx = self.geology._set_data(feature, self.settings['orientation_mode'], self.settings['orientation_datafile'], func="_update_data()")
# # #                 n_orientationy = self.geology._set_data(feature, self.settings['orientation_mode'], self.settings['orientation_datafile'], func="_update_data()")
# # #             else:
# # #                 n_orientationx, n_orientationy = self.geology._set_data(feature, self.settings['orientation_mode'], self.settings['orientation_datafile'], func="_update_data()")
# # #             orientationx = [d for d in self.geology.data if d.label=="orientationx"][-1]
# # #             orientationx.data = n_orientationx
# # #             orientationx.mode = self.settings['orientation_mode']
# # #             orientationy = [d for d in self.geology.data if d.label=="orientationy"][-1]
# # #             orientationy.data = n_orientationy
# # #             orientationy.mode = self.settings['orientation_mode']
# # #             if self.hasMask:
# # #                 self.geology_masked['orientationx'] = ma.MaskedArray(orientationx.data, mask=self.mask_array[:,:,0])
# # #                 self.geology_masked['orientationy'] = ma.MaskedArray(orientationy.data, mask=self.mask_array[:,:,0])
# # #         else:
# # #             data = [d for d in self.geology.data if d.label==feature][-1]
# # #
# # #             if feature == 'fractures' and self.settings['fractures_mode'] == 'random':
# # #                 frac_parameters = [self.settings['fractures_densities'], self.settings['fractures_alpha'], self.settings['fractures_min_orientation'], self.settings['fractures_max_orientation'],
# # #                                self.settings['fractures_min_dip'], self.settings['fractures_max_dip'], self.settings['fractures_min_length'], self.settings['fractures_max_length']]
# # #                 n_data, fractures = self.geology._set_data(feature, self.settings[feature+'_mode'], self.settings[feature+'_datafile'], func="_update_data()", frac_parameters=frac_parameters)
# # #                 data.fractures = fractures
# # #             else:
# # #                 n_data = self.geology._set_data(feature, self.settings[feature+'_mode'], self.settings[feature+'_datafile'], func="_update_data()")
# # #
# # #             data.data = n_data
# # #             data.mode = self.settings[feature+'_mode']
# # #             self.geology.compute_stats_on_data(data)
# # #
# # #             if self.hasMask:
# # #                 self.geology_masked[feature] = ma.MaskedArray(data.data, mask=self.mask_array)
# # #         return None
# # #
# # #     def _update_all(self, feature):
# # #         """
# # #         Update all the parameters.
# # #         """
# # #         self._update_mask('')
# # #         self._update_points('inlets')
# # #         self._update_points('outlets')
# # #         self._update_data('geology')
# # #         self._update_data('topography')
# # #         self._update_data('orientation')
# # #         self._update_data('faults')
# # #         self._update_data('fractures')
# # #         return None
# # #
























# # #     #######
# # #     # DOC #
# # #     #######
# # #
# # #     UPDATE_FUNCTIONS = {
# # #         'mask'        : _update_mask,
# # #         'inlets'      : _update_points,
# # #         'outlets'     : _update_points,
# # #         'geology'     : _update_data,
# # #         'topography'  : _update_data,
# # #         'orientation' : _update_data,
# # #         'faults'      : _update_data,
# # #         'fractures'   : _update_data,
# # #         'all'         : _update_all
# # #         }
# # #
# # #
# # #     PARAMETERS_DOC = """
# # #     'data_has_mask' : bool
# # #         Defines if a mask is used or not.
# # #         If true, a mask must be defined with 'mask_data' parameter.
# # #     'mask_data' : str || array
# # #         Defines the mask vertices.
# # #         Mask datafile path or list of vertices coordinates.
# # #         Useful only when 'data_has_mask' is true.
# # #     'outlets_mode' : str
# # #         Defines the outlets mode.
# # #         'random'    - Full random points
# # #         'import'    - Import points
# # #         'composite' - Add n random points to imported points
# # #     'outlets_data' : str || list
# # #         Defines the outlets.
# # #         Outlets datafile path or list of outlets coordinates.
# # #         Useful only when 'outlets_mode' parameter is on 'import' or 'composite'.
# # #     'outlets_number' : str
# # #         Defines the number of outlets to generate.
# # #         Useful only when 'outlets_mode' parameter is on 'random' or 'composite'.
# # #     'outlets_shuffle' : bool
# # #         Defines whether to shuffle the order of the outlets randomly.
# # #         False - don't shuffle, True - shuffle randomly
# # #         Useful only when iterating over outlets.
# # #     'outlets_importance' : list
# # #         Defines the proportion of outlets to be distributed across each iteration.
# # #         Length of array indicates number of outlet iterations,
# # #         each integer indicates number of outlets to run in that iteration,
# # #         sum of integers = total number of outlets
# # #         [1] - a single iteration with all outlets,
# # #         [1,1,1] - three iterations with one outlet in each,
# # #         [1,2,3] - three iterations with one outlet in the first, 2 outlets in the second, and 3 outlets in the third.
# # #         Useful only when iterating over outlets.
# # #     'inlets_mode' : str
# # #         Defines the inlets mode.
# # #         'random'    - Full random points
# # #         'import'    - Import points
# # #         'composite' - Add n random points to imported points
# # #     'inlets_data' : str || list
# # #         Defines the inlets.
# # #         Inlets datafile path or list of inlets coordinates.
# # #         Useful only when 'inlets_mode' parameter is on 'import' or 'composite'.
# # #     'inlets_number' : str
# # #         Defines the number of inlets to generate.
# # #         Useful only when 'inlets_mode' parameter is on 'random' or 'composite'.
# # #     'inlets_shuffle' : bool
# # #         Defines whether to shuffle the order of the inlets randomly.
# # #         False - don't shuffle, True - shuffle randomly
# # #         Useful only when iterating over inlets.
# # #     'inlets_per_outlet' : list
# # #         Defines the proportion of inlets to be distributed across each outlet.
# # #         Length of array indicates number of outlets,
# # #         each integer indicates number of inlets to assign to that outlet,
# # #         sum of integers = total number of inlets
# # #         [1] - a single iteration with all inlets to one outlet,
# # #         [1,1,1] - three outlets with one inlet in each,
# # #         [1,2,3] - three outlets with one inlet in the first, 2 inlets in the second, and 3 inlets in the third.
# # #         Useful only when iterating over inlets and outlets.
# # #     'inlets_importance' : list
# # #         Defines the proportion of inlets to be distributed across each iteration.
# # #         Length of array indicates number of inlet iterations,
# # #         each integer indicates number of inlets to run in that iteration,
# # #         sum of integers = total number of inlets
# # #         [1] - a single iteration with all inlets,
# # #         [1,1,1] - three iterations with one inlet in each,
# # #         [1,2,3] - three iterations with one inlet in the first, 2 inlets in the second, and 3 inlets in the third.
# # #         Useful only when iterating over inlets.
# # #     'geology_mode' : str
# # #         Defines the geological mode.
# # #         'null'  - No geology
# # #         'gslib' - Import geology via gslib
# # #         'csv'   - Import geology via csv
# # #         'image' - Import geology via image
# # #     'geology_datafile' : str
# # #         Defines the geological datafile path.
# # #         Useful only when 'geology_mode' parameter is not 'null'.
# # #     'topography_mode' : str
# # #         Defines the topography mode.
# # #         'null'  - No topography
# # #         'gslib' - Import topography via gslib
# # #         'csv'   - Import topography via csv
# # #     'topography_datafile' : str
# # #         Defines the topography datafile path.
# # #         Useful only when 'topography_mode' parameter is not 'null'.
# # #     'orientation_mode' : str
# # #         Defines the orientation mode.
# # #         'null'    - No orientation
# # #         'topo'    - Calculate from topography
# # #         'surface' - Calculate from csv file of a surface (useful if using lower surface of karst unit)
# # #     'orientation_datafile' : str
# # #         Defines the orientation datafile path.
# # #         Useful only when 'orientation_mode' parameter is not 'null'.
# # #     'faults_mode' : str
# # #         Defines the mode for the faults.
# # #         'null'  - No faults
# # #         'gslib' - Import faults via gslib
# # #         'csv'   - Import faults via csv
# # #         'image' - Import faults via image
# # #     'faults_datafile' : str
# # #         Defines the faults datafile path.
# # #         Useful only when the 'faults_mode' parameter is on 'gslib', 'csv' or 'image'.
# # #     'fractures_mode' : str
# # #         Defines the mode for the fractures.
# # #         'null'   - No fractures
# # #         'gslib'  - Import fractures via gslib
# # #         'csv'    - Import fractures via csv
# # #         'image'  - Import fractures via image
# # #         'random' - Generate fractures randomly
# # #     'fracture_datafile' : str
# # #         Defines the fractures datafile path.
# # #         Useful only when the 'fractures_mode' parameter is on 'gslib', 'csv' or 'image'.
# # #     'fractures_densities' : list
# # #         Defines the fractures densitiy for each fractures family.
# # #         Useful only when the 'fractures_mode' parameter is on 'random'.
# # #     'fractures_min_orientation' : list
# # #         Defines the minimum orientation of the fracturation for each fractures family.
# # #         Useful only when the 'fractures_mode' parameter is on 'random'.
# # #     'fractures_max_orientation' : list
# # #         Defines the maximum orientation of the fracturation for each fractures family.
# # #         Useful only when the 'fractures_mode' parameter is on 'random'.
# # #     'fractures_min_dip' : list
# # #         Defines the minimum dip of the fracturation for each fractures family.
# # #         Useful only when the 'fractures_mode' parameter is on 'random'.
# # #     'fractures_max_dip' : list
# # #         Defines the maximum dip of the fracturation for each fractures family.
# # #         Useful only when the 'fractures_mode' parameter is on 'random'.
# # #     'fractures_alpha' : list
# # #         Defines alpha, a parameter in the fracturation law.
# # #         Useful only when the 'fractures_mode' parameter is on 'random'.
# # #     'fractures_min_length' : list
# # #         Defines the minimum lenght for all the fractures.
# # #         Useful only when the 'fractures_mode' parameter is on 'random'.
# # #     'fractures_max_length' : list
# # #         Defines the maximum lenght for all the fractures.
# # #         Useful only when the 'fractures_mode' parameter is on 'random'.
# # #     'algorithm' : str
# # #         Defines the algorithm to use when calculating travel time to spring.
# # #         'Isotropic2' - isotropic 2D
# # #         'Isotropic3' - isotropic 3D
# # #         'Riemann2'   - anisotropic 2D
# # #         'Riemann3'   - anisotropic 3D
# # #         See AGD-HFM documentation for full list of options.
# # #     'cost_out' : float, (default: 0.999)
# # #         Defines the fast-marching value for the outside of the study area.
# # #         The value must be between 0 and 1 and should be high to avoid unrealistic conduits.
# # #     'cost_aquifer' : float, (default: 0.3)
# # #         Defines the fast-marching value for the aquifer cells.
# # #         Should be between 0 and 1 and lower than aquiclude but higher than conduits.
# # #     'cost_aquiclude' : float, (default: 0.8)
# # #         Defines the fast-marching value for the aquiclude cells.
# # #         Should be between 0 and 1 and higher than aquiclude but lower than cost_out
# # #     'cost_faults' : float, (default: 0.2)
# # #         Defines the fast-marching value for the faults cells.
# # #         Should be between 0 and 1 and between conduits and cost_out. Higher = conduit will avoid faults, lower = conduits will follow faults
# # #     'cost_fractures' : float, (default: 0.2)
# # #         Defines the fast-marching value for the fractures cells.
# # #         Should be between 0 and 1 and between conduits and cost_out. Higher = conduit will avoid fractures, lower = conduits will follow fractures
# # #     'cost_conduits' : float, (default: 0.01)
# # #         Defines the fast-marching value for the conduits cells.
# # #         Should be between 0 and 1 but lower than aquifer (for conduits to preferentially follow each other)
# # #     'cost_ratio' : float, (default: 0.25)
# # #         Defines the fast-marching ratio of travel cost parallel to gradient / travel cost prependicular to gradient.
# # #         Should be between 0 and 0.5.
# # #     'geology_id' : list
# # #         Defines the geology id (from geology datafile) to consider in the simulation.
# # #         Useful only when the 'geology_mode' parameter is on 'gslib' or 'csv'.
# # #     'rand_seed' : int
# # #         Defines the random seed.
# # #         May help for reproduicity.
# # #     'verbosity' : int
# # #         Define the verbosity (how much output to print during runs).
# # #         0 - print minimal output, 1 - print medium output, 2 - print max output
# # #     """
