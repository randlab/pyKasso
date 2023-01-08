############
### TODO ###
############
# 002 (+)   define default fmm cost
# 003 (---) define a load_project function
# 007 (---) add initial karstic system feature
# 008 (---) add initial field feature
# 009 (---) define the tracers
# 010 - define z function for inlets and outlets
# 012 - Rename logger sections
# 013 - 
# 014 - 
# 015 - 
# 016 - 
# 017 - 
# 018 - 
# 019 - 
# 020 - 

### Others
# Methods / Fonctions / Classes documentation
# better 'debug'/'verbosity' management


####################
### Dependencies ###
####################
import os
import sys
import shutil
import logging
import datetime

### Local dependencies
from ._validations import validate_settings_structure, validate_grid_settings, validate_mask_settings, validate_geologic_feature_settings, validate_points_feature_settings
from .grid import Grid
from .mask import Mask
from .geologic_feature import Geology, Topography, Karst, Field, Faults, Fractures
from .points import generate_random_points, is_point_valid

### External dependencies
import yaml
import numpy as np
import pandas as pd
import karstnet as kn

### Fast-Marching package
import agd
from agd import Eikonal
from agd.Metrics import Riemann

########################
### Module variables ###
########################

this = sys.modules[__name__]
this.ACTIVE_PROJECT = None

# Reference indicators for statistical karstic network analysis
usecols = None
# usecols = "B:I,K"
this.STATISTICS = pd.read_excel(os.path.dirname(os.path.abspath(__file__)) + '/../../misc/' + 'statistics.xlsx', usecols=usecols).describe()

# TODO 002 (+)
# Defines default fast-marching costs
this.cost = {
    'out'       : 0.999,
    'aquifer'   : 0.4,
    'aquiclude' : 0.8,
    'faults'    : 0.2,
    'fractures' : 0.2,
    'conduits'  : 0.1,
    'ratio'     : 0.5,
}

def create_project(project_directory:str):
    """
    TODO
    """
    # Uncomment this to tell apart the development version and the main version
    print('CAUTION: You are using the development version of this package.')

    # Creates the project directories
    core_settings = {
        'project_directory'  : project_directory,
        'inputs_directory'   : project_directory + '/inputs/',
        'outputs_directory'  : project_directory + '/outputs/',
        'settings_directory' : project_directory + '/settings/',
    }
    [os.makedirs(directory, exist_ok=True) for directory in core_settings.values()]

    ### Populates the settings directory with initial yaml settings file
    misc_directory = os.path.dirname(os.path.abspath(__file__)) + '/../../misc/'
    core_settings.update({
        'core_settings_filename' : 'CORE_SETTINGS.yaml',
        'sks_settings_filename'  : 'SKS_SETTINGS.yaml',
        'sim_settings_filename'  : 'SIM_SETTINGS.yaml',
        'log_filename'           : 'pyKasso.log',
    })

    # CORE_SETTINGS
    core_settings_filename = core_settings['settings_directory'] + core_settings['core_settings_filename']
    with open(core_settings_filename, 'w') as f:
        text = "# pyKasso CORE SETTINGS \n---\n"
        yaml_text = yaml.safe_dump(core_settings)
        text = text + yaml_text + "..."
        f.write(text)

    # SKS_SETTINGS
    shutil.copy2(misc_directory + core_settings['sks_settings_filename'], project_directory + '/settings/')

    # SIM_SETTINGS
    shutil.copy2(misc_directory + core_settings['sim_settings_filename'], project_directory + '/settings/')

    # TODO
    this.ACTIVE_PROJECT = this.ACTIVE_PROJECT = {
        'project_directory' : project_directory,
        'n_simulation'      : 0,
        'settings'          : {
            # static features
            'grid'       : None,
            'mask'       : None,
            'topography' : None,
            'geology'    : None,
            'faults'     : None,
            # dynamic features
            'fractures' : None,
            'outlets'   : None,
            'inlets'    : None,
        },
        'model'             : {}
    }

    return None

# TODO 003 (---)
# define a load_project function
def load_project():
    """
    TODO
    """
    return None

#########################################################################################################################################
### SKS CLASS ###
#################

class SKS():
    """
    Class storing all the parameters, the data and the resulting simulated stochastic karst networks.
    TODO
    """

    def __init__(self, core_settings=None, sks_settings=None, sim_settings=None):
        """
        TODO
        """
        ### Loads core, sks and simulation settings
        settings_list  = [core_settings, sks_settings, sim_settings]
        settings_names = ['CORE_SETTINGS', 'SKS_SETTINGS', 'SIM_SETTINGS']
        for settings, settings_name in zip(settings_list, settings_names):
            if settings is None:
                settings = this.ACTIVE_PROJECT['project_directory'] + '/settings/{}.YAML'.format(settings_name)
            if isinstance(settings, str):
                with open(settings, 'r') as f:
                    setattr(self, settings_name.upper(), yaml.safe_load(f))
            if isinstance(settings, dict):
                setattr(self, settings_name.upper(), settings)
        # TODO - 
        # test validity of arguments inputs

        ### Creates simulation directory
        this.ACTIVE_PROJECT['n_simulation'] += 1
        self.CORE_SETTINGS['simulation_directory'] = self.CORE_SETTINGS['outputs_directory'] + 'simulation_{}/'.format(this.ACTIVE_PROJECT['n_simulation'])
        os.makedirs(self.CORE_SETTINGS['simulation_directory'], exist_ok=True)

        ### Resets and creates log file
        if this.ACTIVE_PROJECT['n_simulation'] == 1:

            # Resets logging module
            root = logging.getLogger()
            list(map(root.removeHandler, root.handlers))
            list(map(root.removeFilter,  root.filters))

            # Sets new log file
            log_file = self.CORE_SETTINGS['outputs_directory'] + self.CORE_SETTINGS['log_filename']
            logging.basicConfig(
                filename=log_file, 
                encoding='utf-8', 
                level=logging.INFO, 
                filemode="w",
                format=' %(name)-21s | %(levelname)-8s | %(message)s'
            )
            # this.logger.setLevel(logging.INFO)
            this.logger = logging.getLogger("sks")

        # Prints current simulation number
        l = len(str(this.ACTIVE_PROJECT['n_simulation']))
        this.logger.info('***********************{}****'.format('*' * l))
        this.logger.info('*** pyKasso simulation {} ***'.format(this.ACTIVE_PROJECT['n_simulation']))
        this.logger.info('***********************{}****'.format('*' * l))
        this.logger.info(datetime.datetime.now())
        this.logger.info('---')

        ### Validates sks and simulation dictionary settings structure
        self.SKS_SETTINGS = validate_settings_structure(self.SKS_SETTINGS, 'SKS')
        self.SIM_SETTINGS = validate_settings_structure(self.SIM_SETTINGS, 'SIM')
        

##########################################################################################################################################################################
### WRAPPERS ###
################

    def _decorator_logging(feature_name):
        def _(function):
            # print(function.__name__)
            # @wraps(function)
            def _wrapper_logging(*args, **kwargs):
                # TODO - this logger (feature_name)
                try:
                    result = function(*args, **kwargs)
                except:
                    this.logger.error(function.__name__)
                    raise
                else:
                    this.logger.info(function.__name__)
                    return result
                finally:
                    pass
            return _wrapper_logging
        return _

##########################################################################################################################################################################
### CONSTRUCTING STATIC AND DYNAMIC FEATURES ###
################################################

    def build_model(self):
        """
        TODO
        """
        ##################################
        ### Constructs static features ###
        ##################################


        # 01 - Constructs the grid
        if self.SKS_SETTINGS['grid'] != this.ACTIVE_PROJECT['settings']['grid']:
            self.SKS_SETTINGS['grid'] = validate_grid_settings(self.SKS_SETTINGS['grid'])
            self._construct_feature_grid()
            this.ACTIVE_PROJECT['settings']['grid'] = self.SKS_SETTINGS['grid']
            this.ACTIVE_PROJECT['model']['grid']    = self.GRID
        else:
            self.GRID = this.ACTIVE_PROJECT['model']['grid']


        # 02 - Constructs the mask
        # TODO - Si la grille change, alors le mask doit être recontrôlé
        if self.SKS_SETTINGS['mask'] != this.ACTIVE_PROJECT['settings']['mask']:
            self.SKS_SETTINGS['mask'] = validate_mask_settings(self.SKS_SETTINGS['mask'])
            self._construct_feature_mask()
            this.ACTIVE_PROJECT['settings']['mask'] = self.SKS_SETTINGS['mask']
            this.ACTIVE_PROJECT['model']['mask']    = self.MASK
        else:
            self.MASK = this.ACTIVE_PROJECT['model']['mask']


        # 03 - Constructs the topography
        if self.SKS_SETTINGS['topography'] != this.ACTIVE_PROJECT['settings']['topography']:
            self.SKS_SETTINGS['topography'] = validate_geologic_feature_settings('topography', self.SKS_SETTINGS['topography'], self.GRID)
            self._construct_feature_topography()
            this.ACTIVE_PROJECT['settings']['topography'] = self.SKS_SETTINGS['topography']
            this.ACTIVE_PROJECT['model']['topography']    = self.TOPOGRAPHY
        else:
            self.TOPOGRAPHY = this.ACTIVE_PROJECT['model']['topography']


        # 04 - Constructs the geology
        if self.SKS_SETTINGS['geology'] != this.ACTIVE_PROJECT['settings']['geology']:
            self.SKS_SETTINGS['geology'] = validate_geologic_feature_settings('geology', self.SKS_SETTINGS['geology'], self.GRID)
            self._construct_feature_geology()
            this.ACTIVE_PROJECT['settings']['geology'] = self.SKS_SETTINGS['geology']
            this.ACTIVE_PROJECT['model']['geology']    = self.GEOLOGY
        else:
            self.GEOLOGY = this.ACTIVE_PROJECT['model']['geology']

        
        # 05 - Constructs the faults
        if self.SKS_SETTINGS['faults'] != this.ACTIVE_PROJECT['settings']['faults']:
            self.SKS_SETTINGS['faults'] = validate_geologic_feature_settings('faults', self.SKS_SETTINGS['faults'], self.GRID)
            self._construct_feature_faults()
            this.ACTIVE_PROJECT['settings']['faults'] = self.SKS_SETTINGS['faults']
            this.ACTIVE_PROJECT['model']['faults']    = self.FAULTS
        else:
            self.FAULTS = this.ACTIVE_PROJECT['model']['faults']
 

        # TODO 007
        # add initial karstic system feature
        # 06 - Constructs the initial karst network
        # self._construct_feature_karst()


        # TODO 008
        # add initial field feature
        # 07 - Constructs the field
        # self._construct_feature_field()


        ###################################
        ### Constructs dynamic features ###
        ###################################

        
        # 08 - Sets the seeds
        self._set_rng()
    

        # 09 - Constructs the outlets and shuffles
        if self.SIM_SETTINGS['outlets'] != this.ACTIVE_PROJECT['settings']['outlets']:
            self.SIM_SETTINGS['outlets'] = validate_points_feature_settings('outlets', self.SIM_SETTINGS['outlets'])
            self.OUTLETS = self._construct_feature_points('outlets')
            if self.SIM_SETTINGS['outlets']['shuffle']:
                self.OUTLETS = self.OUTLETS.sample(frac=1, random_state=self.RNG['master']).reset_index(drop=True)
            this.ACTIVE_PROJECT['settings']['outlets'] = self.SIM_SETTINGS['outlets']
            this.ACTIVE_PROJECT['model']['outlets']    = self.OUTLETS
        else:
            self.OUTLETS = this.ACTIVE_PROJECT['model']['outlets']

        
        # 10 - Constructs the inlets and shuffles
        if self.SIM_SETTINGS['inlets'] != this.ACTIVE_PROJECT['settings']['inlets']:
            self.SIM_SETTINGS['inlets'] = validate_points_feature_settings('inlets', self.SIM_SETTINGS['inlets'])
            self.INLETS = self._construct_feature_points('inlets')
            if self.SIM_SETTINGS['inlets']['shuffle']:
                self.INLETS = self.INLETS.sample(frac=1, random_state=self.RNG['master']).reset_index(drop=True)
            this.ACTIVE_PROJECT['settings']['inlets'] = self.SIM_SETTINGS['inlets']
            this.ACTIVE_PROJECT['model']['inlets']    = self.INLETS
        else:
            self.INLETS = this.ACTIVE_PROJECT['model']['inlets']


        # TODO 009
        # define the tracers
        # 11 - Constructs the tracers
        # self._construct_feature_tracers()


        # 12 - Constructs the fractures
        # if self.SIM_SETTINGS['fractures'] != this.ACTIVE_PROJECT['settings']['fractures']:
        #     self._construct_feature_fractures()
        #     this.ACTIVE_PROJECT['settings']['fractures'] = self.SIM_SETTINGS['fractures']
        #     this.ACTIVE_PROJECT['model']['fractures']    = self.FRACTURES
        # else:
        #     self.FRACTURES = this.ACTIVE_PROJECT['model']['fractures']

        return None

##########################################################################################################################################################################
### STATIC FEATURES ###
#######################


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


    @_decorator_logging('topography')
    def _construct_feature_topography(self):
        """
        Constructs the topography.
        """
        if self.SKS_SETTINGS['topography']['data'] != '':
            self.TOPOGRAPHY = Topography(**self.SKS_SETTINGS['topography'], grid=self.GRID)
            self.TOPOGRAPHY._compute_topographic_surface()
            self.TOPOGRAPHY._compute_statistics(self.GRID)
        else:
            self.TOPOGRAPHY = None
        return None


    @_decorator_logging('geology')
    def _construct_feature_geology(self):
        """
        Constructs the geology.
        """
        self.GEOLOGY = Geology(**self.SKS_SETTINGS['geology'], grid=self.GRID)

        if self.TOPOGRAPHY is None:
            topography = np.where(self.GEOLOGY.data > 0, 1, 0)
            self.TOPOGRAPHY = Topography(data=topography, grid=self.GRID)
            self.TOPOGRAPHY._compute_topographic_surface()
            self.TOPOGRAPHY._compute_statistics(self.GRID)
        else:
            self.GEOLOGY.data = np.where(self.TOPOGRAPHY.data == 1, self.GEOLOGY.data, 0)
        
        self.GEOLOGY._compute_surface(self.TOPOGRAPHY.data)
        self.GEOLOGY._compute_statistics(self.GRID)
        return None


    @_decorator_logging('faults')
    def _construct_feature_faults(self):
        """
        TODO
        """
        if self.SKS_SETTINGS['faults']['data'] != '':
            self.FAULTS = Faults(**self.SKS_SETTINGS['faults'], grid=self.GRID)
            self.FAULTS.data = np.where(self.TOPOGRAPHY.data == 1, self.FAULTS.data, 0)
            self.FAULTS._compute_surface(self.TOPOGRAPHY.data)
            self.FAULTS._compute_statistics(self.GRID)
        else:
            self.FAULTS = None

        return None


    # TODO 007
    # @_decorator_logging('karst')
    # gestion des réseaux incomplets
    # https://scipy-lectures.org/packages/scikit-image/auto_examples/plot_labels.html
    # def _construct_feature_karst(self):
    #     """
    #     Constructs the initial karst conduits network.
    #     """
    #     if self.SKS_SETTINGS['karst']['data'] != '':
    #         self.KARST = Karst(**self.SKS_SETTINGS['karst'], grid=self.GRID)
    #         self.KARST.data = np.where(self.TOPOGRAPHY.data == 1, self.KARST.data, 0)
    #         self.KARST._compute_surface(self.TOPOGRAPHY.data)
    #         self.KARST._compute_statistics(self.GRID)
    #     else:
    #         self.KARST = None
    #     return None


    # TODO 008
    # @_decorator_logging('field')
    # def _construct_feature_field(self):
    #     """
    #     Constructs the field.
    #     """
    #     if self.SKS_SETTINGS['field']['data'] != '':
    #         self.FIELD = Field(**self.SKS_SETTINGS['field'], grid=self.GRID)
    #     else:
    #         self.FIELD = None
    #     return None


##########################################################################################################################################################################
### DYNAMIC FEATURES ###
########################


    @_decorator_logging('seed')
    def _set_rng(self):
        """
        Sets the seed(s).
        TODO
        """

        if ('seed' not in self.SIM_SETTINGS) or (self.SIM_SETTINGS['seed'] == 0):
            self.SIM_SETTINGS['seed'] = np.random.default_rng().integers(low=0, high=10**6)

        attributes = ['inlets', 'outlets', 'fractures']
        for attribute in attributes:
            if ('seed' not in self.SIM_SETTINGS[attribute]) or (self.SIM_SETTINGS[attribute]['seed'] == 0):
                self.SIM_SETTINGS[attribute]['seed'] = np.random.default_rng().integers(low=0, high=10**6)

        self.RNG = {
            'master'    : np.random.default_rng(self.SIM_SETTINGS['seed']),
            'inlets'    : np.random.default_rng(self.SIM_SETTINGS['inlets']['seed']),
            'outlets'   : np.random.default_rng(self.SIM_SETTINGS['outlets']['seed']),
            'fractures' : np.random.default_rng(self.SIM_SETTINGS['fractures']['seed']),
        }

        return None


    @_decorator_logging('points')
    def _construct_feature_points(self, kind):
        """
        TODO
        Constructs the inlets / outlets.

        Four cases

        1. No points declared
        2. More points required than provided : Generates additional random points
        3. Less points required than provided : Pick random points among provided ones
        4. Points required equals points declared
        """

        ### Get existing points

        # Loads points if needed
        if isinstance(self.SIM_SETTINGS[kind]['data'], (str)):
            self.SIM_SETTINGS[kind]['data'] = np.genfromtxt(self.SIM_SETTINGS[kind]['data'])

        # Inspects validity of points
        points = self.SIM_SETTINGS[kind]['data']
        tests  = map(lambda point : is_point_valid(point, self.GRID, self.MASK), points)
        validated_points = [point for (point, test) in zip(points, tests) if test == True]
        diff = len(self.SIM_SETTINGS[kind]['data']) - len(validated_points)
        if diff > 0:
            # TODO - LOG - VERBOSITY
            this.logger.warning('{}/{} {} have been discarded because out of domain.'.format(diff, len(self.SIM_SETTINGS[kind]['data']), kind))
        self.SIM_SETTINGS[kind]['data'] = validated_points


        ### Get new points according to the right case

        # Case 1 - No points declared
        if (self.SIM_SETTINGS[kind]['data'] == '') or (self.SIM_SETTINGS[kind]['data'] == []):
            points = generate_random_points(self.SIM_SETTINGS[kind]['number'], self.RNG[kind], self.GRID, self.MASK, self.GEOLOGY)

        # Case 2 - More points required than provided
        elif (self.SIM_SETTINGS[kind]['number'] > len(self.SIM_SETTINGS[kind]['data'])):
            n_points = self.SIM_SETTINGS[kind]['number'] - len(self.SIM_SETTINGS[kind]['data'])
            points = self.SIM_SETTINGS[kind]['data'] + generate_random_points(n_points, self.RNG[kind], self.GRID, self.MASK, self.GEOLOGY)

        # Case 3 - Less points required than provided
        elif (self.SIM_SETTINGS[kind]['number'] < len(self.SIM_SETTINGS[kind]['data'])):
            points = self.RNG['master'].choice(self.SIM_SETTINGS[kind]['data'], self.SIM_SETTINGS[kind]['number'], replace=False)

        # Case 4 - Points required equals points declared
        else:
            points = self.SIM_SETTINGS[kind]['data']


        ### Populates the DataFrame
        if len(points[0]) == 2:
            x, y = zip(*points)

            ### TODO 010
            # Define z function for inlets and outlets
            # z_func ?
            i, j = self.GRID.get_i(x), self.GRID.get_j(y)
            k = self.TOPOGRAPHY.surface_indices[i,j]

            if kind == 'inlets':
                z = self.GRID.get_z(k)
                z = self.GRID.zmax
                # z_surf = z + self.GRID.dz / 2
            if kind == 'outlets':
                z = self.GRID.z0
                z = self.GRID.zmin
                # z_surf = z - self.GRID.dz / 2

        if len(points[0]) == 3:
            x, y, z = zip(*points)

        data = {
            'x' : x,
            'y' : y,
            'z' : z,
        }
        return pd.DataFrame(data=data)


    # TODO 009
    @_decorator_logging('tracers')
    def _construct_feature_tracers(self):
        """
        Constructs the tracers.
        """
        self.TRACERS = {}
        return None


    @_decorator_logging('fractures')
    def _construct_feature_fractures(self):
        """
        TODO
        """
        if self.SIM_SETTINGS['fractures']['data'] != '':
            self.FRACTURES = Fractures(**self.SIM_SETTINGS['fractures'], grid=self.GRID, rng=self.RNG['fractures'])
            self.FRACTURES.data = np.where(self.TOPOGRAPHY.data == 1, self.FRACTURES.data, 0)
            self.FRACTURES._compute_surface(self.TOPOGRAPHY.data)
            self.FRACTURES._compute_statistics(self.GRID)
        else:
            self.FRACTURES = None
        return None


##########################################################################################################################################################################
### KARST NETWORK SIMULATION ###
################################

    def compute_karst_networks(self):
        """
        TODO
        """

        # 1 - Initializes parameters from farst-marching method
        self._initialize_karst_network_parameters()

        # # TODO
        # if 'fmm' in self.SKS_SETTINGS['debug']:
        #     if self.SKS_SETTINGS['debug']['fmm']:
        #         # TODO
        #         import pykasso.visualization as pkv
        #         pkv.debug_plot_initialize(self)

        # 2 - Compute karst network
        self._compute_karst_network()

        return None


    # 1 - Initializes parameters from farst-marching method
    def _initialize_karst_network_parameters(self):
        """
        Initialize the karst network parameters.
        """
        ### Inlets - Outlets - Iterations

        # Defining some variables
        outlets_nbr = len(self.OUTLETS)
        inlets_nbr  = len(self.INLETS)
        self.outlets_importance = self.SIM_SETTINGS['outlets']['importance']
        self.inlets_importance  = self.SIM_SETTINGS['inlets']['importance']
        inlets_per_outlet       = self.SIM_SETTINGS['inlets']['per_outlet']

        # Calculating inlets and outlets repartitions
        self.nbr_iteration  = len(self.outlets_importance) * len(self.inlets_importance)        # total number of iterations that will occur
        outlets_repartition = self._repartition_points(outlets_nbr, self.outlets_importance)    # correct for outlets_importance not summing to correct number of actual outlets
        inlets_repartition  = self._repartition_points(inlets_nbr , inlets_per_outlet)          # correct for inlets_per_outlet not summing to correct number of actual inlets

        # print('outlets_repartition')
        # print(outlets_repartition)
        # print('inlets_repartition')
        # print(inlets_repartition)
        
        # Distributing inlets and outlets iterations
        outlets_distribution = pd.Series([k for (k, n) in enumerate(outlets_repartition) for j in range(n)], name='outlet_iteration')
        inlets_distribution  = pd.Series([k for (k, n) in enumerate(inlets_repartition)  for j in range(n)], name='outlet_key')
        self._outlets = pd.concat([self.OUTLETS, outlets_distribution], axis=1) # store as a semi-private variable for internal use only
        self._inlets  = pd.concat([self.INLETS , inlets_distribution] , axis=1) # store as a semi-private variable for internal use only

        # Distributing iterations for each inlet
        for (outlet_key, row) in self._outlets.iterrows():
            inlets_test    = self._inlets['outlet_key']==outlet_key
            inlets_current = self._inlets[inlets_test]
            inlets_nbr     = len(inlets_current)
            repartition  = self._repartition_points(inlets_nbr, self.inlets_importance)
            distribution = pd.Series([k for (k, n) in enumerate(repartition) for j in range(n)], name='inlet_iteration', index=inlets_current.index)
            self._inlets.loc[inlets_test, 'inlet_iteration'] = distribution

        # print('outlets')
        # print(self._outlets)
        # print('inlets')
        # print(self._inlets)

        ### Raster maps
        self.maps = {
            'outlets' : np.full((self.GRID.nx, self.GRID.ny, self.GRID.nz), np.nan), # map of null values where each cell with an outlet will have the index of that outlet
            'nodes'   : np.full((self.GRID.nx, self.GRID.ny, self.GRID.nz), np.nan), # map of null values where each cell that has a node will be updated with that node index
            'cost'    : np.zeros((self.nbr_iteration, self.GRID.nx, self.GRID.ny, self.GRID.nz)), # cost of travel through each cell
            'alpha'   : np.zeros((self.nbr_iteration, self.GRID.nx, self.GRID.ny, self.GRID.nz)), # cost of travel along gradient through each cell
            'beta'    : np.zeros((self.nbr_iteration, self.GRID.nx, self.GRID.ny, self.GRID.nz)), # cost of travel perpendicular to gradient through each cell
            'time'    : np.zeros((self.nbr_iteration, self.GRID.nx, self.GRID.ny, self.GRID.nz)), # travel time to outlet from each cell
            'karst'   : np.zeros((self.nbr_iteration, self.GRID.nx, self.GRID.ny, self.GRID.nz)), # presence/absence of karst conduit in each cell
        }

        ### Vector maps
        # TODO - self.vector ??
        self.nodes     = {} #empty dic to store nodes (key: nodeID, val: [x, y, z, type])
        self.edges     = {} #empty dic to store edges (key: edgeID, val: [inNode, outNode])
        self.n         = 0  #start node counter at zero
        self.e         = 0  #start edge counter at zero
        self.geodesics = [] #empty list to store raw fast-marching path output

        ### Set up fast-marching:
        # Note: AGD-HFM library has different indexing, so model dimensions must be [nz, ny, nx],
        # and model extent must be [zmin,zmax, ymin,ymax, xmin,xmax] (NOT x0,y0,z0)
        # TODO
        # if self.FIELD is None:
        #     pass
        # else:
        #     pass
        self.algorithm = self.SKS_SETTINGS['fmm']['algorithm']
        self.riemannMetric = []                   # this changes at every iteration, but cannot be stored?
        self.fastMarching = agd.Eikonal.dictIn({
            'model'             : self.algorithm,      # set algorithm from settings file ('Isotropic2', 'Isotropic3', 'Riemann2', 'Riemann3')
            'order'             : 2,              # recommended setting: 2 (replace by variable)
            'exportValues'      : 1,              # export the travel time field
            'exportGeodesicFlow': 1               # export the walker paths (i.e. the conduits)
        })
        self.fastMarching.SetRect(                              # give the fast-marching algorithm the model grid
            sides=[[self.GRID.xmin, self.GRID.xmax],            # bottom edge,   top edge (NOT centerpoint)
                   [self.GRID.ymin, self.GRID.ymax],            # leftmost edge, rightmost edge (NOT centerpoint)
                   [self.GRID.zmin, self.GRID.zmax]],           
            dims=[self.GRID.nx, self.GRID.ny, self.GRID.nz])    # number of cells, number of cells, number of cells

        return None


    def _repartition_points(self, nbr_points, importance):
        """
        Correct for integers in importance factors list not summing correctly to total number of points.
        """
        total_importance = float(sum(importance))                                               # total number of points as assigned (this may not be the actual total)
        proportion       = [float(i)/total_importance       for i in importance]                # percent of points to use this iteration
        repartition      = [round(proportion[i]*nbr_points) for i in range(len(proportion))]    # number of points to use this iteration (need to round bc percentage will not result in whole number)
        repartition[-1]  = nbr_points-sum(repartition[0:-1])                                    # leftover points are assignd to last iteration
        return repartition



    # 2 - Compute karst network
    def _compute_karst_network(self):
        """
        Compute the karst network according to the parameters.

        Append the results to the `karst_simulations` attribute as an element of an array.

        Parameters
        ----------
        TODO
        """
        ### 2.1 - Computes conduits for each generation & store nodes and edges for network
        self._compute_iterations_karst_network()

        # if self.SETTINGS['debug']['costmap']:
        #     self.debug_plot_fmm_feature('costmap')
        # if self.SETTINGS['debug']['timemap']:
        #     self.debug_plot_fmm_feature('timemap')
        # if self.SETTINGS['debug']['karstmap']:
        #     self.debug_plot_fmm_feature('karstmap')

#         ### ii. Calculates the karst network statistics indicators with karstnet and save karst network
#         # karstnet_edges = list(self.edges.values()) #convert to format needed by karstnet (list)
#         # karstnet_nodes = copy.deepcopy(self.nodes) #convert to format needed by karstnet (dic with only coordinates) - make sure to only modify a copy!
#         # for key, value in karstnet_nodes.items():  #drop last item in list (the node type) for each dictionary entry
#         #     value.pop()

#         # # Computes karstnet indicators
#         # k = kn.KGraph(karstnet_edges, karstnet_nodes)  #make graph - edges must be a list, and nodes must be a dic of format {nodeindex: [x,y]}
#         # stats = k.characterize_graph(verbose)

#         # ### iii. Store all the relevant data for this network in dictionaries:
#         # maps = copy.deepcopy(self.maps)                  #store copy of maps
#         # points = {}
#         # points['inlets']  = copy.deepcopy(self._inlets)  #store copy of inlets in pandas df format with info about iteration and inlet/outlet assignment (this will be used in plotting later)
#         # points['outlets'] = copy.deepcopy(self._outlets) #store copy of outlets in pandas df format with info about iteration and inlet/outlet assignment (this will be used in plotting later)
#         # network = {}
#         # network['edges'] = copy.deepcopy(self.edges) #store copy of edges list
#         # network['nodes'] = copy.deepcopy(self.nodes) #store copy of nodes list
#         # network['karstnet'] = copy.deepcopy(k)       #store copy of karstnet network object (including graph)
#         # config = copy.deepcopy(self.SETTINGS)        #store copy of settings for the run being stored
#         # try:
#         #     del config['debug']
#         # except:
#         #     pass
#         # self.SIMULATIONS.append(KarstNetwork(maps, points, network, stats, config))

#         # TODO ??????????
#         ### iv. - Return inlets and outlets to their original format:
#         # #Chloe: this is to convert inlets and outlets back from pandas dataframes
#         # #if the private self._inlets and self._outlets variables are working correctly, this step may be removed
#         # self.inlets  = np.asarray(self._inlets)[:,0:2]    #return inlets to array format without info about iteration and inlet/outlet assignment (important for later iterations)
#         # self.outlets = np.asarray(self._outlets)[:,0:2]   #return outlets to array format without info about iteration and inlet/outlet assignment (important for later iterations)
#         # return None




    # 2.1 - Computes conduits for each generation
    def _compute_iterations_karst_network(self):
        """
        Compute each generation of karst conduits.

        TODO : à terminer
        """
        # TODO
        # if self.SKS_SETTINGS['debug']['verbosity'] > 0:
        #     print('-START-')

        # Define outlets map according to outlets emplacements
        for (i, outlet) in self._outlets.iterrows(): #assign outlet indices. Compute the outlets map (array indicating location of outlets as their index and everywhere else as nan).
            X = self.GRID.get_i(outlet.x)
            Y = self.GRID.get_j(outlet.y)
            Z = self.GRID.get_k(outlet.z)
            self.maps['outlets'][X][Y][Z] = i

        # TODO
        # Plot for debugging:
        # if self.SETTINGS['verbosity'] > 1:
            # f, ax = plt.subplots(1, 1, figsize=(10, 10))
            # ax.imshow(np.transpose(self.DATA['geology']['geology'].data[:,:,0], (1,0)), extent=self.GRID.extent, origin='lower', cmap='gray_r', alpha=0.5)

        # Set up iteration structure:
        iteration = 0

        # Iterate in outlets groups
        for outlet_iteration in range(len(self.outlets_importance)):

            # if self.SETTINGS['verbosity'] > 2:
                # print('Total Iteration:', iteration, 'Outlet iteration:', outlet_iteration)

            outlets_current = self._outlets[self._outlets['outlet_iteration'] == outlet_iteration]
            # if self.SETTINGS['verbosity'] > 1:
                # ax.scatter(outlets_current.x, outlets_current.y, c='c') # Debugging

            lst = outlets_current.index.values.tolist()
            inlets_current = self._inlets[self._inlets['outlet_key'].isin(lst)]

            # Iterate in inlets groups
            for inlet_iteration in range(len(self.inlets_importance)):
                inlets_current_iteration = inlets_current[inlets_current['inlet_iteration'] == inlet_iteration]
                self._outlets.loc[outlets_current.index, 'iteration'] = iteration          # Store total iteration counter
                self._inlets.loc [inlets_current_iteration.index, 'iteration'] = iteration # Store total iteration counter

                #Compute travel time maps and conduit network
                if self.algorithm == 'Isotropic3':
                    # 2.1.1
                    self._compute_cost_map(iteration)
                    # 2.1.2
                    self._compute_time_map(iteration)
                    # 2.1.3
                    self._compute_karst_map(iteration)
                elif self.algorithm == 'Riemann3':
                    # 2.1.1
                    self._compute_cost_map(iteration)
                    # 2.1.4
                    self._compute_alpha_map(iteration)
                    # 2.1.5
                    self._compute_beta_map(iteration)
                    # TODO
                    # notebooks mirebeau
                    # 2.1.6
                    self._compute_riemann_metric(iteration)
                    # 2.1.7
                    self._compute_time_map(iteration)
                    # 2.1.3
                    self._compute_karst_map(iteration)

                # TODO - this part in the verification functions
                else:
                    print('Unrecognized algorithm:', self.algorithm)
                
                iteration = iteration + 1   #increment total iteration number by 1
                print(iteration)

                # TODO
                # if self.SETTINGS['verbosity'] > 0:
                    # print('iteration:{}/{}'.format(iteration, self.nbr_iteration))
        # if self.SETTINGS['verbosity'] > 0:
            # print('- END -')

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


    ### 2.1.1 Iso- and anisotropic case
    def _compute_cost_map(self, iteration):
        """
        Compute the cost map (how difficult it is to traverse each cell).

        TODO
        """
        # If it's the first iteration, iniatialize the cost map according to the geological settings.
        if iteration == 0:
            
            ### Geology
            for key, cost in self.GEOLOGY.cost.items():
                self.maps['cost'][0] = np.where(self.GEOLOGY.data == key, cost, self.maps['cost'][0])

            ### Faults
            for key, cost in self.FAULTS.cost.items():
                self.maps['cost'][0] = np.where(self.FAULTS.data == key, cost, self.maps['cost'][0])

            ### TODO
            ### Karst
            # for key, cost in self.KARST.cost.items():
            #     self.maps['cost'][0] = np.where(self.KARST.data == key, cost, self.maps['cost'][0])

            ### TODO
            ### Fractures
            # for key, cost in self.FRACTURES.cost.items():
                # self.maps['cost'][0] = np.where(self.FRACTURES.data == key, cost, self.maps['cost'][0])

            ### If out of mask
            if self.MASK is not None:
                self.maps['cost'][0] = np.where(self.MASK.mask == 1, this.cost['out'], self.maps['cost'][0])

        # If it's not the first iteration
        else: 
            self.maps['cost'][iteration] = self.maps['cost'][iteration-1] # set cost map to be same as previous iteration
            self.maps['cost'][iteration] = np.where(self.maps['karst'][iteration-1] > 0, this.cost['conduits'], self.maps['cost'][iteration]) # where karst conduits are present from previous iteration, set cost to conduit cost, elsewhere, leave unchanged
        
        return None


    ### 2.1.4 Anisotropic case
    def _compute_alpha_map(self, iteration):
        """
        Compute the alpha map: travel cost in the same direction as the gradient.
        Cost map * topography map, so that the cost is higher at higher elevations, encouraging conduits to go downgradient.

        TODO : à terminer
        """
        self.maps['alpha'][iteration] = self.maps['cost'][iteration] * self.GRID.Z

        # TODO
        # if self.settings['topography_mode'] != 'null':
        #     self.maps['alpha'][iteration] = self.maps['cost'][iteration] * self.data['topography'].data
        # elif self.settings['orientation_mode'] == 'surface':
        #     self.maps['alpha'][iteration] = self.maps['cost'][iteration] * self.data['surface'].data
        # else:
        #     self.maps['alpha'][iteration] = self.maps['cost'][iteration]

        return None


    ### 2.1.5 Anisotropic case
    # TODO - _set_beta_map()
    # anisotropie avec pendage
    # objet Surface
    def _compute_beta_map(self, iteration):
        """
        Compute the beta map: travel cost perpendicular to the gradient.
        If beta is higher than alpha, conduits will follow the steepest gradient.
        If beta is lower than alpha, conduits will follow contours.
        """
        self.maps['beta'][iteration] = self.maps['alpha'][iteration] / self.SKS_SETTINGS['fmm']['cost']['ratio']

        # TODO - alpha = beta dans zone saturée
        return None


    ### 2.1.6 Anisotropic case
    def _compute_riemann_metric(self, iteration):
        """
        Compute the riemann metric: Define the Riemannian metric needed as input for the anisotropic fast marching.

        TODO : à terminer
        """
        x, y, z = self.GRID.Z.shape
        # orientationx = np.transpose(np.gradient(self.GRID.Z, axis=0), (2,1,0))
        # orientationy = np.transpose(np.gradient(self.GRID.Z, axis=1), (2,1,0))
        # orientationz = np.transpose(np.gradient(self.GRID.Z, axis=2), (2,1,0))
        orientationx = np.gradient(self.GRID.Z, axis=0)
        orientationy = np.gradient(self.GRID.Z, axis=1)
        orientationz = np.gradient(self.GRID.Z, axis=2)

        # TODO
        # if (z == 1):
        #     orientationz = np.transpose(self.GRID.Z, (2,1,0)) # TODO ??!!
        # else:
        #     orientationz = np.transpose(np.gradient(self.GRID.Z, axis=2), (2,1,0))

        # alpha = np.transpose(self.maps['alpha'][iteration], (2,1,0))
        # beta  = np.transpose(self.maps['beta'][iteration],  (2,1,0))
        alpha = self.maps['alpha'][iteration]
        beta  = self.maps['beta'][iteration]
        # self.riemannMetric = agd.Metrics.Riemann.needle([orientationz, orientationy, orientationx], alpha, beta) 
        self.riemannMetric = agd.Metrics.Riemann.needle([orientationx, orientationy, orientationz], alpha, beta) 
        # TODO gamma ??? 
        return None


    ### 2.1.2
    def _compute_time_map(self, iteration):
        """
        Compute the travel time map (how long it takes to get to the outlet from each cell),
        using the ani- or isotropic agd-hfm fast-marching algorithm, and store travel time map.
        Note: the AGD-HFM library uses different indexing, so x and y indices are reversed for inlets and outlets.
        TODO
        """
        # Set the outlets for this iteration
        seeds_x = self._outlets[self._outlets['iteration']==iteration].x
        seeds_y = self._outlets[self._outlets['iteration']==iteration].y
        seeds_z = self._outlets[self._outlets['iteration']==iteration].z
        seeds   = list(zip(seeds_x, seeds_y, seeds_z))
        self.fastMarching['seeds'] = seeds

        # Select inlets for current iteration
        tips_x = self._inlets[self._inlets['iteration']==iteration].x
        tips_y = self._inlets[self._inlets['iteration']==iteration].y
        tips_z = self._inlets[self._inlets['iteration']==iteration].z
        tips   = list(zip(tips_x, tips_y, tips_z))
        self.fastMarching['tips'] = tips

        # Set the travel cost through each cell
        if self.algorithm == 'Isotropic3':
            self.fastMarching['cost'] = self.maps['cost'][iteration]
        if self.algorithm == 'Riemann3':
            self.fastMarching['metric'] = self.riemannMetric

        # Set verbosity of hfm run
        self.fastMarching['verbosity'] = self.SKS_SETTINGS['verbosity']['agd']

        # Run the fast marching algorithm and store the outputs
        self.fastMarchingOutput = self.fastMarching.Run()

        # Store travel time maps
        self.maps['time'][iteration] = self.fastMarchingOutput['values']

        # Store fastest travel paths
        self.geodesics.append(self.fastMarchingOutput['geodesics'])
        return None


    ### 2.1.3 
    def _compute_karst_map(self, iteration):
        """
        Compute the karst map based on the paths from agd-hfm.
        Array of all zeros, with ones in cells containing a karst conduit.
        """
        # Get karst map from previous iteration (except for the very first iteration)
        if iteration > 0:
            self.maps['karst'][iteration] = self.maps['karst'][iteration-1]

        # Debugging plot:
        # Chloe: this should stay in since it is very useful if there are problems
        # f1, ax1 = plt.subplots(1, 1, figsize=(10, 10))

        ### Loop over conduit paths generated by fast marching:
        for path in self.fastMarchingOutput['geodesics']:   #loop over conduit paths in this iteration (there is one path from each inlet)
            merge = False                                   #reset indicator for whether this conduit has merged with an existing conduit
            for p in range(path.shape[1]):                  #loop over points making up this conduit path
                point = path[:,p]                           #get coordinates of current point
                # print('p', point)
                [[ix, iy, iz], error] = self.fastMarching.IndexFromPoint(point) #convert to coordinates to indices, /!\ returning iy first then ix
                # print(ix, iy, iz)
                # ax1.scatter(point[1], point[0], c='g',s=5)  #debugging

                ############## WHAT IS HAPPENING HERE?
                # if ix < 0 or iy < 0 or iz < 0:
                #     print(ix,iy,iz)
                #     continue

                #Place nodes and links:
                if np.isnan(self.maps['nodes'][ix, iy, iz]):                                    #if there is no existing conduit node here
                    if ~np.isnan(self.maps['outlets'][ix, iy, iz]):                              #if there is an outlet here (cell value is not nan)
                        outlet = self._outlets.iloc[int(self.maps['outlets'][ix, iy, iz])]         #get the outlet coordinates using the ID in the outlets map
                        self.nodes[self.n]             = [outlet.x, outlet.y, outlet.z, 'outfall']           #add a node at the outlet coordinates (with the node type for SWMM)
                        self.maps['nodes'][ix, iy, iz] = self.n                                   #update node map with node index
                        # ax1.scatter(outlet.x, outlet.y, marker='o', c='b')                   #debugging
                        if p > 0:                                                           #if this is not the first point (i.e. the inlet) in the current path
                            if merge == False:                                               #if this conduit has not merged with an existing conduit
                                self.edges[self.e] = [self.n-1, self.n]                       #add an edge connecting the previous node to the current node
                                self.e = self.e + 1                                             #increment edge counter up by one
                                # ax1.plot((self.nodes[self.n][0], self.nodes[self.n-1][0]), (self.nodes[self.n][1], self.nodes[self.n-1][1]))
                            else:                                                          #if this conduit HAS merged with an existing conduit
                                [[fromix, fromiy, fromiz], error] = self.fastMarching.IndexFromPoint(path[:, p-1]) #get xyz indices of previous point in current conduit path
                                n_from = self.maps['nodes'][fromix, fromiy, fromiz]           #get node index of the node already in the cell where the previous point was
                                self.edges[self.e] = [n_from, self.n]                         #add an edge connecting existing conduit node to current node
                                self.e = self.e + 1                                             #increment edge counter up by one
                                # ax1.plot((self.nodes[self.n].x, self.nodes[n_from].x), (self.nodes[self.n].y, self.nodes[n_from].y))
                        self.n = self.n + 1                                                   #increment node counter up by one
                    else:                                                                  #if there is NOT an outlet here
                        if p > 0:                                                           #if this is not the first point in the current path
                            #possible improvement: if the next point on the path is on an existing point, skip the current point.
                            self.nodes[self.n] = [point[0], point[1], point[2], 'junction']            #add a junction node here (with the node type for SWMM)
                            self.maps['nodes'][ix, iy, iz] = self.n                               #update node map with node index
                            #ax1.scatter(point[1],point[0], marker='o', edgecolor='g', facecolor='none')   #debugging
                            if merge == False:                                              #if this conduit has not merged with an existing conduit
                                self.edges[self.e] = [self.n-1, self.n]                      #add and edge connecting the previous node to the current node
                                self.e = self.e + 1                                            #increment edge counter up by one
                                #ax1.plot((self.nodes[self.n][1], self.nodes[self.n-1][1]),(self.nodes[self.n][0], self.nodes[self.n-1][0]), c='gold', marker=None)
                            else:                                                           #if this conduit HAS merged with an existing conduit
                                [[fromix, fromiy, fromiz], error] = self.fastMarching.IndexFromPoint(path[:, p-1]) #get xy indices of previous point in current conduit path
                                n_from = self.maps['nodes'][fromix, fromiy, fromiz]                   #get node index of the node already in the cell where the previous point was
                                self.edges[self.e] = [n_from, self.n]                        #add an edge connecting existing conduit node to current node
                                self.e = self.e + 1                                            #increment edge counter up by one
                                merge = False                                                #reset merge indicator to show that current conduit has left                                                              #if this is the first point in current path
                        else:                                                                #if this is the first point in the current path (counter <= 0, therefore it is an inlet)
                            self.nodes[self.n] = [point[0], point[1], point[2], 'inlet']               #add an inlet node here (with the node type for SWMM)
                            self.maps['nodes'][ix, iy, iz] = self.n                               #update node map with node index
                            #ax1.scatter(point[1],point[0], marker='o', edgecolor='sienna', facecolor='none')
                        self.n = self.n + 1                                                   #increment node counter up by one
                elif ~np.isnan(self.maps['nodes'][ix, iy, iz]):                                 #if there is already a node in this cell (either because there is a conduit here, or because there are two nodes in the same cell)
                    n_existing = self.maps['nodes'][ix, iy, iz]                                  #get index of node already present in current cell
                    if merge == True:                                                       #if this conduit has already merged into an existing conduit
                        pass                                                                 #skip this node (there is already a node here)
                    elif n_existing == self.n-1:                                            #if existing index is only one less than next node to be added index, this is a duplicate node and can be skipped
                        pass                                                                 #skip this node (duplicate)
                    else:                                                                   #if existing node index is >1 less than next node to be added index
                        if p > 0:                                                           #if this is not the first point in the current path
                            self.edges[self.e] = [self.n-1, n_existing]                      #add an edge connecting most recently added node and existing node in cell
                            self.e = self.e + 1                                                #increment edge counter up by one
                            merge = True                                                     #add a flag indicating that this conduit has merged into an existing one
                            #ax1.plot((self.nodes[self.n-1][1], self.nodes[n_existing][1]),(self.nodes[self.n-1][0], self.nodes[n_existing][0]), c='r', marker=None)
                        else:                                                                #if this is the first point in the current path (i.e. the inlet is on an exising conduit)
                            self.nodes[self.n] = [point[0], point[1], point[2], 'inlet']                #add a node here (with the node type for SWMM)- this will cause there to be two nodes in the same cell
                            self.maps['nodes'][ix, iy, iz] = self.n                                #update node map with node index
                            self.n = self.n + 1                                                 #increment node counter by 1
                            #ax1.scatter(point[1],point[0], marker='o', edgecolor='g', facecolor='none')  #debugging
                self.maps['karst'][iteration][ix, iy, iz] = 1                               #update karst map to put a conduit in current cell


        ### Debugging plot
        # # Display inlets and outlets
        # ax1.scatter(self._outlets.x, self._outlets.y, c='cyan',   s=100)
        # ax1.scatter(self._inlets.x,  self._inlets.y,  c='orange', s=100)
        # ax1.scatter(self._outlets[self._outlets.iteration==iteration].x, self._outlets[self._outlets.iteration==iteration].y, c='cyan',   s=100)
        # ax1.scatter(self._inlets[self._inlets.iteration==iteration].x,   self._inlets[self._inlets.iteration==iteration].y,   c='orange', s=100)
        # # Display karst network
        # ax1.imshow(np.transpose(self.maps['karst'][iteration], (1,0,2)), origin='lower', extent=self.GRID.extent, cmap='gray_r')
        # ax1.imshow(np.transpose(self.maps['nodes'], (1,0,2)), origin='lower', extent=self.GRID.extent, cmap='gray_r')
        return None

# #     ######################
# #     ### Karst Analysis ###
# #     ######################

# #     def compare_stats(self, mean=False):
# #         """
# #         Compare statistics between reference indicators and calculated networks.

# #         TODO
# #         param 'iteration=0'

# #         """
# #         indicators = ['cpd', 'cv degree', 'cv length', 'orientation entropy', 'length entropy', 'aspl', 'mean degree', 'mean length', 'correlation vertex degree']
# #         stats = pd.DataFrame(columns=indicators)

# #         for (i, karst_network) in enumerate(self.SIMULATIONS):
# #             stats.loc[i] = karst_network.stats

# #         # Apply style
# #         def _bg_color(x, min_val, max_val):
# #             if (x < min_val) or (x > max_val):
# #                 return 'background-color: red'
# #             else:
# #                 return 'background-color: green'

# #         display(stats.style.applymap(_bg_color, min_val = self.reference_statistics['cpd']['min'], max_val = self.reference_statistics['cpd']['max'], subset = ['cpd'])\
# #         .applymap(_bg_color, min_val = self.reference_statistics['cv degree']['min'], max_val = self.reference_statistics['cv degree']['max'], subset = ['cv degree'])\
# #         .applymap(_bg_color, min_val = self.reference_statistics['cv length']['min'], max_val = self.reference_statistics['cv length']['max'], subset = ['cv length'])\
# #         .applymap(_bg_color, min_val = self.reference_statistics['orientation entropy']['min'], max_val = self.reference_statistics['orientation entropy']['max'], subset = ['orientation entropy'])\
# #         .applymap(_bg_color, min_val = self.reference_statistics['length entropy']['min'], max_val = self.reference_statistics['length entropy']['max'], subset = ['length entropy'])\
# #         .applymap(_bg_color, min_val = self.reference_statistics['aspl']['min'], max_val = self.reference_statistics['aspl']['max'], subset = ['aspl'])\
# #         .applymap(_bg_color, min_val = self.reference_statistics['mean degree']['min'], max_val = self.reference_statistics['mean degree']['max'], subset = ['mean degree'])\
# #         .applymap(_bg_color, min_val = self.reference_statistics['mean length']['min'], max_val = self.reference_statistics['mean length']['max'], subset = ['mean length'])\
# #         .applymap(_bg_color, min_val = self.reference_statistics['correlation vertex degree']['min'], max_val = self.reference_statistics['correlation vertex degree']['max'], subset = ['correlation vertex degree']))
# #         return None


# #     def compute_average_paths(self, mask=0):
# #         """
# #         TODO
# #         Compute the mean of all the simulations.
# #         """

# #         # Calculate the average from all the simulations
# #         karst_maps = []
# #         for karst_simulation in self.SIMULATIONS:
# #             data = karst_simulation.maps['karst'][-1]
# #             karst_maps.append(data)
# #         karst_prob = sum(karst_maps)/len(karst_maps)

# #         self.karst_prob = karst_prob
# #         return karst_prob

# #     def show_average_paths(self):
# #         """
# #         todo
# #         """
# #         ### Call the plotter
# #         p = pv.Plotter(notebook=False)

# #         ### Construct the grid
# #         vtk = pv.UniformGrid()
# #         vtk.dimensions = np.array((self.GRID.nx, self.GRID.ny, self.GRID.nz)) + 1
# #         vtk.origin     = (self.GRID.x0 - self.GRID.dx/2, self.GRID.y0 - self.GRID.dy/2, self.GRID.z0 - self.GRID.dz/2)
# #         vtk.spacing    = (self.GRID.dx, self.GRID.dy, self.GRID.dz)

# #         vtk['values'] = self.karst_prob.flatten(order="F")

# #         mesh = vtk.cast_to_unstructured_grid()
# #         ghosts = np.argwhere(vtk['values'] < 1.0)
# #         mesh.remove_cells(ghosts)
# #         p.add_mesh(mesh, show_edges=False)

# #         ### Plotting
# #         # p.add_title(feature)
# #         p.add_axes()
# #         bounds = p.show_bounds(mesh=vtk)
# #         p.add_actor(bounds)
# #         p.show(cpos='xy')

# #         return None




# ######################################################################################################################################################
# #
# #     #######
# #     # DOC #
# #     #######
# #
# #     PARAMETERS_DOC = """
# #     'data_has_mask' : bool
# #         Defines if a mask is used or not.
# #         If true, a mask must be defined with 'mask_data' parameter.
# #     'mask_data' : str || array
# #         Defines the mask vertices.
# #         Mask datafile path or list of vertices coordinates.
# #         Useful only when 'data_has_mask' is true.
# #     'outlets_mode' : str
# #         Defines the outlets mode.
# #         'random'    - Full random points
# #         'import'    - Import points
# #         'composite' - Add n random points to imported points
# #     'outlets_data' : str || list
# #         Defines the outlets.
# #         Outlets datafile path or list of outlets coordinates.
# #         Useful only when 'outlets_mode' parameter is on 'import' or 'composite'.
# #     'outlets_number' : str
# #         Defines the number of outlets to generate.
# #         Useful only when 'outlets_mode' parameter is on 'random' or 'composite'.
# #     'outlets_shuffle' : bool
# #         Defines whether to shuffle the order of the outlets randomly.
# #         False - don't shuffle, True - shuffle randomly
# #         Useful only when iterating over outlets.
# #     'outlets_importance' : list
# #         Defines the proportion of outlets to be distributed across each iteration.
# #         Length of array indicates number of outlet iterations,
# #         each integer indicates number of outlets to run in that iteration,
# #         sum of integers = total number of outlets
# #         [1] - a single iteration with all outlets,
# #         [1,1,1] - three iterations with one outlet in each,
# #         [1,2,3] - three iterations with one outlet in the first, 2 outlets in the second, and 3 outlets in the third.
# #         Useful only when iterating over outlets.
# #     'inlets_mode' : str
# #         Defines the inlets mode.
# #         'random'    - Full random points
# #         'import'    - Import points
# #         'composite' - Add n random points to imported points
# #     'inlets_data' : str || list
# #         Defines the inlets.
# #         Inlets datafile path or list of inlets coordinates.
# #         Useful only when 'inlets_mode' parameter is on 'import' or 'composite'.
# #     'inlets_number' : str
# #         Defines the number of inlets to generate.
# #         Useful only when 'inlets_mode' parameter is on 'random' or 'composite'.
# #     'inlets_shuffle' : bool
# #         Defines whether to shuffle the order of the inlets randomly.
# #         False - don't shuffle, True - shuffle randomly
# #         Useful only when iterating over inlets.
# #     'inlets_per_outlet' : list
# #         Defines the proportion of inlets to be distributed across each outlet.
# #         Length of array indicates number of outlets,
# #         each integer indicates number of inlets to assign to that outlet,
# #         sum of integers = total number of inlets
# #         [1] - a single iteration with all inlets to one outlet,
# #         [1,1,1] - three outlets with one inlet in each,
# #         [1,2,3] - three outlets with one inlet in the first, 2 inlets in the second, and 3 inlets in the third.
# #         Useful only when iterating over inlets and outlets.
# #     'inlets_importance' : list
# #         Defines the proportion of inlets to be distributed across each iteration.
# #         Length of array indicates number of inlet iterations,
# #         each integer indicates number of inlets to run in that iteration,
# #         sum of integers = total number of inlets
# #         [1] - a single iteration with all inlets,
# #         [1,1,1] - three iterations with one inlet in each,
# #         [1,2,3] - three iterations with one inlet in the first, 2 inlets in the second, and 3 inlets in the third.
# #         Useful only when iterating over inlets.
# #     'geology_mode' : str
# #         Defines the geological mode.
# #         'null'  - No geology
# #         'gslib' - Import geology via gslib
# #         'csv'   - Import geology via csv
# #         'image' - Import geology via image
# #     'geology_datafile' : str
# #         Defines the geological datafile path.
# #         Useful only when 'geology_mode' parameter is not 'null'.
# #     'topography_mode' : str
# #         Defines the topography mode.
# #         'null'  - No topography
# #         'gslib' - Import topography via gslib
# #         'csv'   - Import topography via csv
# #     'topography_datafile' : str
# #         Defines the topography datafile path.
# #         Useful only when 'topography_mode' parameter is not 'null'.
# #     'orientation_mode' : str
# #         Defines the orientation mode.
# #         'null'    - No orientation
# #         'topo'    - Calculate from topography
# #         'surface' - Calculate from csv file of a surface (useful if using lower surface of karst unit)
# #     'orientation_datafile' : str
# #         Defines the orientation datafile path.
# #         Useful only when 'orientation_mode' parameter is not 'null'.
# #     'faults_mode' : str
# #         Defines the mode for the faults.
# #         'null'  - No faults
# #         'gslib' - Import faults via gslib
# #         'csv'   - Import faults via csv
# #         'image' - Import faults via image
# #     'faults_datafile' : str
# #         Defines the faults datafile path.
# #         Useful only when the 'faults_mode' parameter is on 'gslib', 'csv' or 'image'.
# #     'fractures_mode' : str
# #         Defines the mode for the fractures.
# #         'null'   - No fractures
# #         'gslib'  - Import fractures via gslib
# #         'csv'    - Import fractures via csv
# #         'image'  - Import fractures via image
# #         'random' - Generate fractures randomly
# #     'fracture_datafile' : str
# #         Defines the fractures datafile path.
# #         Useful only when the 'fractures_mode' parameter is on 'gslib', 'csv' or 'image'.
# #     'fractures_densities' : list
# #         Defines the fractures densitiy for each fractures family.
# #         Useful only when the 'fractures_mode' parameter is on 'random'.
# #     'fractures_min_orientation' : list
# #         Defines the minimum orientation of the fracturation for each fractures family.
# #         Useful only when the 'fractures_mode' parameter is on 'random'.
# #     'fractures_max_orientation' : list
# #         Defines the maximum orientation of the fracturation for each fractures family.
# #         Useful only when the 'fractures_mode' parameter is on 'random'.
# #     'fractures_min_dip' : list
# #         Defines the minimum dip of the fracturation for each fractures family.
# #         Useful only when the 'fractures_mode' parameter is on 'random'.
# #     'fractures_max_dip' : list
# #         Defines the maximum dip of the fracturation for each fractures family.
# #         Useful only when the 'fractures_mode' parameter is on 'random'.
# #     'fractures_alpha' : list
# #         Defines alpha, a parameter in the fracturation law.
# #         Useful only when the 'fractures_mode' parameter is on 'random'.
# #     'fractures_min_length' : list
# #         Defines the minimum lenght for all the fractures.
# #         Useful only when the 'fractures_mode' parameter is on 'random'.
# #     'fractures_max_length' : list
# #         Defines the maximum lenght for all the fractures.
# #         Useful only when the 'fractures_mode' parameter is on 'random'.
# #     'algorithm' : str
# #         Defines the algorithm to use when calculating travel time to spring.
# #         'Isotropic2' - isotropic 2D
# #         'Isotropic3' - isotropic 3D
# #         'Riemann2'   - anisotropic 2D
# #         'Riemann3'   - anisotropic 3D
# #         See AGD-HFM documentation for full list of options.
# #     'cost_out' : float, (default: 0.999)
# #         Defines the fast-marching value for the outside of the study area.
# #         The value must be between 0 and 1 and should be high to avoid unrealistic conduits.
# #     'cost_aquifer' : float, (default: 0.3)
# #         Defines the fast-marching value for the aquifer cells.
# #         Should be between 0 and 1 and lower than aquiclude but higher than conduits.
# #     'cost_aquiclude' : float, (default: 0.8)
# #         Defines the fast-marching value for the aquiclude cells.
# #         Should be between 0 and 1 and higher than aquiclude but lower than cost_out
# #     'cost_faults' : float, (default: 0.2)
# #         Defines the fast-marching value for the faults cells.
# #         Should be between 0 and 1 and between conduits and cost_out. Higher = conduit will avoid faults, lower = conduits will follow faults
# #     'cost_fractures' : float, (default: 0.2)
# #         Defines the fast-marching value for the fractures cells.
# #         Should be between 0 and 1 and between conduits and cost_out. Higher = conduit will avoid fractures, lower = conduits will follow fractures
# #     'cost_conduits' : float, (default: 0.01)
# #         Defines the fast-marching value for the conduits cells.
# #         Should be between 0 and 1 but lower than aquifer (for conduits to preferentially follow each other)
# #     'cost_ratio' : float, (default: 0.25)
# #         Defines the fast-marching ratio of travel cost parallel to gradient / travel cost prependicular to gradient.
# #         Should be between 0 and 0.5.
# #     'geology_id' : list
# #         Defines the geology id (from geology datafile) to consider in the simulation.
# #         Useful only when the 'geology_mode' parameter is on 'gslib' or 'csv'.
# #     'rand_seed' : int
# #         Defines the random seed.
# #         May help for reproduicity.
# #     'verbosity' : int
# #         Define the verbosity (how much output to print during runs).
# #         0 - print minimal output, 1 - print medium output, 2 - print max output
# #     """
