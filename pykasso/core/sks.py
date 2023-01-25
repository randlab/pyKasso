############
### TODO ###
############
# 002 (+)   define default fmm cost
# 003 (---) define a load_project function
# 007 (---) add initial karstic system feature
# 008 (---) add initial field feature
# 009 (---) define the tracers
# 015 - 
# 016 - 
# 017 - 
# 018 - 
# 019 - 
# 020 - 

### Others
# Methods / Fonctions / Classes documentation

####################
### Dependencies ###
####################
import os
import sys
import shutil
import logging
import datetime

### Local dependencies
import pykasso.core._wrappers as wp
from .grid import Grid
from .mask import Mask
from .geologic_feature import Geology, Topography, Karst, Field, Faults, Fractures
from .points import generate_coordinate, is_point_valid

### External dependencies
import yaml
import numpy as np
import pandas as pd
import karstnet as kn

### Fast-Marching package
import agd
from agd import Eikonal
from agd.Metrics import Riemann

################################################################################################################################################
### Module variables ###
########################

this = sys.modules[__name__]
this.ACTIVE_PROJECT = None
this.feature_kind = None

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

################################################################################################################################################
### Module functions ###
########################

def create_project(project_directory):
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

    ### Creates a dictionary used for settings comparison between actual and previous simulation 
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
        'model' : {}
    }

    return None

# TODO 003 (---) - define a load_project function
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

    def __init__(self, sks_settings=None, core_settings=None, debug_level=None):
        """
        TODO
        """
        # TODO
        # test validity of arguments inputs
        # -> 'verbosity'

        ### Loads core and sks settings
        settings_list  = [core_settings, sks_settings]
        settings_names = ['CORE_SETTINGS', 'SKS_SETTINGS']
        for settings, settings_name in zip(settings_list, settings_names):
            if settings is None:
                settings = this.ACTIVE_PROJECT['project_directory'] + '/settings/{}.YAML'.format(settings_name)
            if isinstance(settings, str):
                with open(settings, 'r') as f:
                    setattr(self, settings_name.upper(), yaml.safe_load(f))
            if isinstance(settings, dict):
                setattr(self, settings_name.upper(), settings)
        
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
            
            # Remove PIL package log entries
            logging.getLogger("PIL.PngImagePlugin").setLevel(logging.CRITICAL + 1)

            # Sets logging level
            levels = {
                0 : logging.DEBUG,
                1 : logging.INFO,
                2 : logging.WARNING,
                3 : logging.ERROR,
                4 : logging.CRITICAL,
            }
            level = levels[self.SKS_SETTINGS['verbosity']['logging']]

            # Sets new logging file
            log_file = self.CORE_SETTINGS['outputs_directory'] + self.CORE_SETTINGS['log_filename']
            logging.basicConfig(
                filename=log_file, 
                encoding='utf-8', 
                level=level,
                filemode="w",
                format=' %(name)-23s | %(levelname)-8s | %(message)s'
            )
            
            # Prints basic information in the log
            # - pyKasso + version
            # - dossier du projet
            # - logging level
            
            
        # Prints current simulation number
        this.logger = logging.getLogger(".")
        l = len(str(this.ACTIVE_PROJECT['n_simulation']))
        this.logger.info('***********************{}****'.format('*' * l))
        this.logger.info('*** pyKasso simulation {} ***'.format(this.ACTIVE_PROJECT['n_simulation']))
        this.logger.info('***********************{}****'.format('*' * l))
        this.logger.info('Date : ' + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        this.logger.info('---')

        # Sets debug level
        if debug_level is None:
            self.debug_mode  = False
            self.debug_level = {
                'model'      : 15,
                'simulation' : 10,
                'iteration'  : 1000 # TODO
            }
        else:
            self.debug_mode  = True
            self.debug_level = debug_level
            if 'model' not in self.debug_level:
                self.debug_level['model'] = 15
            if 'simulation' not in self.debug_level:
                self.debug_level['simulation'] = 10
            if 'iteration' not in self.debug_level:
                self.debug_level['iteration'] = 1000
            this.logger.info('Debug mode : {}'.format(self.debug_mode))
            this.logger.info('Debug level :')
            this.logger.info('-> model          : {}'.format(self.debug_level['model']))
            this.logger.info('-> simulation     : {}'.format(self.debug_level['simulation']))
            this.logger.info('-> max. iteration : {}'.format(self.debug_level['iteration']))
            this.logger.info('---')
        

##########################################################################################################################################################################
### BUILDING THE MODEL ###
##########################

    def build_model(self):
        """
        TODO
        """
        ### Constructs geologic static features
        self._build_model_grid()        # 01 - Builds the grid
        self._build_model_mask()        # 02 - Builds the mask
        self._build_model_topography()  # 03 - Builds the topography
        self._build_model_geology()     # 04 - Builds the geology
        self._build_model_faults()      # 05 - Builds the faults
        self._build_model_karst()       # 06 - Builds the initial karst network
        self._build_model_field()       # 07 - Builds the field

        ### Constructs geologic dynamic features 
        self._build_model_seed()        # 08 - Sets the main seed
        self._build_model_outlets()     # 09 - Builds the outlets and shuffles them
        self._build_model_inlets()      # 10 - Builds the inlets and shuffles them
        self._build_model_tracers()     # 11 - Builds the tracers
        self._build_model_fractures()   # 12 - Builds the fractures

        ### Constructs fmm features 
        self._construct_fmm_variables() # 13 - Constructs fmm variables
        return None
    

    # TODO - Grille inchangeable ???
    # TODO - Si la grille change, alors le mask doit être recontrôlé, etc pour les autres features
    # 1
    @wp._debug_level(1)
    @wp._parameters_validation('grid', 'required')
    @wp._memoize('grid')
    @wp._logging()
    def _build_model_grid(self):
        """
        Builds the grid.
        """
        self.grid = Grid(**self.SKS_SETTINGS['grid'])
        return None
    
    # 2
    @wp._debug_level(2)
    @wp._parameters_validation('mask', 'optional')
    @wp._memoize('mask')
    @wp._logging()
    def _build_model_mask(self):
        """
        Builds the mask.
        """
        if self.SKS_SETTINGS['mask']['data'] != '':
            self.mask = Mask(**self.SKS_SETTINGS['mask'])
            self.mask.validate_vertices(self.grid)
        else:
            self.mask = None
        return None
    
    # 3
    @wp._debug_level(3)
    @wp._parameters_validation('topography', 'optional')
    @wp._memoize('topography')
    @wp._logging()
    def _build_model_topography(self):
        """
        Builds the topography.
        """
        if self.SKS_SETTINGS['topography']['data'] != '':
            self.topography = Topography(**self.SKS_SETTINGS['topography'], grid=self.grid)
            self.topography._compute_topographic_surface()
            self.topography._compute_statistics(self.grid)
        else:
            self.topography = None
        return None
        
    # 4
    @wp._debug_level(4)
    @wp._parameters_validation('geology', 'optional')
    @wp._memoize('geology')
    @wp._logging()
    def _build_model_geology(self):
        """
        Builds the geology.
        """
        self.geology = Geology(**self.SKS_SETTINGS['geology'], grid=self.grid)

        if self.topography is None:
            topography = np.where(self.geology.data > 0, 1, 0)
            self.topography = Topography(data=topography, grid=self.grid)
            self.topography._compute_topographic_surface()
            self.topography._compute_statistics(self.grid)
        else:
            self.geology.data = np.where(self.topography.data == 1, self.geology.data, 0)
        
        self.geology._compute_surface(self.topography.data)
        self.geology._compute_statistics(self.grid)
        return None
    
    # 5
    @wp._debug_level(5)
    @wp._parameters_validation('faults', 'optional')
    @wp._memoize('faults')
    @wp._logging()
    def _build_model_faults(self):
        """
        Builds the faults.
        """
        if self.SKS_SETTINGS['faults']['data'] != '':
            self.faults = Faults(**self.SKS_SETTINGS['faults'], grid=self.grid)
            self.faults.data = np.where(self.topography.data == 1, self.faults.data, 0)
            self.faults._compute_surface(self.topography.data)
            self.faults._compute_statistics(self.grid)
        else:
            self.faults = None
        return None
    
    # 6
    # TODO
    # gestion des réseaux incomplets
    # https://scipy-lectures.org/packages/scikit-image/auto_examples/plot_labels.html
    @wp._debug_level(6)
    # @wp._parameters_validation('karst', 'optional')
    # @wp._memoize('karst')
    # @wp._logging()
    def _build_model_karst(self):
        """
        Builds the initial karst conduits network.
        """
        # if self.SKS_SETTINGS['karst']['data'] != '':
        #     self.karst = Karst(**self.SKS_SETTINGS['karst'], grid=self.grid)
        #     self.karst.data = np.where(self.topography.data == 1, self.karst.data, 0)
        #     self.karst._compute_surface(self.topography.data)
        #     self.karst._compute_statistics(self.grid)
        # else:
        #     self.karst = None
        return None
    
    # 7
    # TODO
    @wp._debug_level(7)
    # @wp._parameters_validation('field', 'optional')
    # @wp._memoize('field')
    # @wp._logging()
    def _build_model_field(self):
        """
        Builds the field.
        """
        # if self.SKS_SETTINGS['field']['data'] != '':
        #     self.field = Field(**self.SKS_SETTINGS['field'], grid=self.grid)
        # else:
        #     self.field = None
        return None
    
    # 8
    @wp._debug_level(8)
    @wp._parameters_validation('sks', 'required') # todo - rename ?
    def _build_model_seed(self):
        """ 
        TODO
        """
        if self.SKS_SETTINGS['sks']['seed'] == 0:
            self.SKS_SETTINGS['sks']['seed'] = np.random.default_rng().integers(low=0, high=10**6)
        
        self.RNG = {
            'master' : np.random.default_rng(self.SKS_SETTINGS['sks']['seed']),
        }
        return None

    # 9
    @wp._debug_level(9)
    @wp._parameters_validation('outlets', 'required')
    @wp._memoize('outlets')
    @wp._logging()
    def _build_model_outlets(self):
        """ 
        TODO
        """
        self._set_rng('outlets')
        self.outlets = self._construct_feature_points('outlets')
        if self.SKS_SETTINGS['outlets']['shuffle']:
            self.outlets = self.outlets.sample(frac=1, random_state=self.RNG['master']).reset_index(drop=True)
        return None

    # 10
    @wp._debug_level(10)
    @wp._parameters_validation('inlets', 'required')
    @wp._memoize('inlets')
    @wp._logging()
    def _build_model_inlets(self):
        """ 
        TODO
        """
        self._set_rng('inlets')
        self.inlets = self._construct_feature_points('inlets')
        if self.SKS_SETTINGS['inlets']['shuffle']:
            self.inlets = self.inlets.sample(frac=1, random_state=self.RNG['master']).reset_index(drop=True)
        return None

    # 11
    @wp._debug_level(11)
    # @wp._parameters_validation('tracers', 'optional')
    # @wp._memoize('tracers')
    # @wp._logging()
    def _build_model_tracers(self):
        """
        Builds the tracers.
        """
        # self.tracers = {}
        return None
    
    # 12
    @wp._debug_level(12)
    @wp._parameters_validation('fractures', 'optional')
    @wp._memoize('fractures')
    @wp._logging()
    def _build_model_fractures(self):
        """ 
        Builds the fractures
        """
        # TODO ????? d'ou vient la clé 'axis' ??
        self._set_rng('fractures')
        if self.SKS_SETTINGS['fractures']['data'] != '':
            self.fractures = Fractures(**self.SKS_SETTINGS['fractures'], grid=self.grid, rng=self.RNG['fractures'])
            self.fractures.data = np.where(self.topography.data == 1, self.fractures.data, 0)
            self.fractures._compute_surface(self.topography.data)
            self.fractures._compute_statistics(self.grid)
        else:
            self.fractures = None
        return None
    
    # 13
    @wp._debug_level(13)
    @wp._parameters_validation('fmm', 'required')
    @wp._logging('fmm', 'construction')
    def _construct_fmm_variables(self):
        """
        TODO
        """
        self._initialize_fmm_iterations()   # Distributes inlets and outlets among karstic generations
        self._construct_fmm_iterations()    # Distributes inlets and outlets among computing iterations
        self._initialize_fmm_variables()    # Initializes useful variables for the farst-marching method
        return None
    
    
    ########################
    ### BUILDING METHODS ###
    ########################

    def _set_rng(self, attribute):
        """
        Sets the corresponding seed.
        TODO
        """
        if self.SKS_SETTINGS[attribute]['seed'] == 0:
            self.SKS_SETTINGS[attribute]['seed'] = np.random.default_rng().integers(low=0, high=10**6)
        self.RNG[attribute] = np.random.default_rng(self.SKS_SETTINGS[attribute]['seed'])
        return None

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

        # TODO - propositions de mot clé pour les fonctions lambda
        
        lambda_functions = {
            'x' : self.SKS_SETTINGS[kind]['x'],
            'y' : self.SKS_SETTINGS[kind]['y'],
            'z' : self.SKS_SETTINGS[kind]['z'],
        }
        geologic_ids = self.SKS_SETTINGS[kind]['geology']

        ### Get existing points

        # Loads points if needed
        if isinstance(self.SKS_SETTINGS[kind]['data'], (str)) and not (self.SKS_SETTINGS[kind]['data'] == ''):
            self.SKS_SETTINGS[kind]['data'] = np.genfromtxt(self.SKS_SETTINGS[kind]['data'])

        # Inspects validity of points
        points = self.SKS_SETTINGS[kind]['data']
        validated_points = []
        for point in points:
            if is_point_valid(point, self.grid, self.mask):

                # Inspects completeness of point
                if len(point) == 2:
                    x, y = point
                    point = generate_coordinate(lambda_functions, self.RNG[kind], self.grid, self.mask, self.geology, xi=x, yi=y, geologic_ids=geologic_ids)
                    # TODO - Log

                validated_points.append(point)

        diff = len(points) - len(validated_points)
        if diff > 0:
            # TODO - LOG - VERBOSITY
            # TODO - log name is not correct
            msg = '{}/{} {} have been discarded because out of domain.'.format(diff, len(self.SKS_SETTINGS[kind]['data']), kind)
            this.logger.warning(msg)
        self.SKS_SETTINGS[kind]['data'] = validated_points


        ### Get new points according to the right case

        # Case 1 - No points declared
        if (self.SKS_SETTINGS[kind]['data'] == '') or (self.SKS_SETTINGS[kind]['data'] == []):
            points = [generate_coordinate(lambda_functions, self.RNG[kind], self.grid, self.mask, self.geology, geologic_ids=geologic_ids) for _ in range(self.SKS_SETTINGS[kind]['number'])]

        # Case 2 - More points required than provided
        elif (self.SKS_SETTINGS[kind]['number'] > len(self.SKS_SETTINGS[kind]['data'])):
            n_points = self.SKS_SETTINGS[kind]['number'] - len(self.SKS_SETTINGS[kind]['data'])
            points = self.SKS_SETTINGS[kind]['data'] + [generate_coordinate(lambda_functions, self.RNG[kind], self.grid, self.mask, self.geology, geologic_ids=geologic_ids) for _ in range(n_points)]

        # Case 3 - Less points required than provided
        elif (self.SKS_SETTINGS[kind]['number'] < len(self.SKS_SETTINGS[kind]['data'])):
            points = self.RNG['master'].choice(self.SKS_SETTINGS[kind]['data'], self.SKS_SETTINGS[kind]['number'], replace=False)

        # Case 4 - Points required equals points declared
        else:
            points = self.SKS_SETTINGS[kind]['data']

        ### Populates the DataFrame
        x, y, z = zip(*points)
        data = {
            'x' : x,
            'y' : y,
            'z' : z,
        }
        return pd.DataFrame(data=data)


    ####################
    ### FMM FEATURES ###
    ####################

    def _initialize_fmm_iterations(self):
        """
        TODO
        """
        # Defining some variables
        outlets_nbr = len(self.outlets)
        inlets_nbr  = len(self.inlets)
        self.outlets_importance = self.SKS_SETTINGS['outlets']['importance']
        self.inlets_importance  = self.SKS_SETTINGS['inlets']['importance']
        inlets_per_outlet       = self.SKS_SETTINGS['inlets']['per_outlet']

        # Calculating inlets and outlets repartitions
        self.nbr_iteration  = len(self.outlets_importance) * len(self.inlets_importance)        # total number of iterations that will occur
        outlets_repartition = self._repartition_points(outlets_nbr, self.outlets_importance)    # correct for outlets_importance not summing to correct number of actual outlets
        inlets_repartition  = self._repartition_points(inlets_nbr , inlets_per_outlet)          # correct for inlets_per_outlet not summing to correct number of actual inlets

        # Distributing inlets and outlets iterations
        outlets_distribution = pd.Series([k for (k, n) in enumerate(outlets_repartition) for j in range(n)], name='outlet_iteration')
        inlets_distribution  = pd.Series([k for (k, n) in enumerate(inlets_repartition)  for j in range(n)], name='outlet_key')
        self._outlets = pd.concat([self.outlets, outlets_distribution], axis=1) # store as a semi-private variable for internal use only
        self._inlets  = pd.concat([self.inlets , inlets_distribution] , axis=1) # store as a semi-private variable for internal use only

        # Distributing iterations for each inlet
        for (outlet_key, row) in self._outlets.iterrows():
            inlets_test    = self._inlets['outlet_key']==outlet_key
            inlets_current = self._inlets[inlets_test]
            inlets_nbr     = len(inlets_current)
            repartition  = self._repartition_points(inlets_nbr, self.inlets_importance)
            distribution = pd.Series([k for (k, n) in enumerate(repartition) for j in range(n)], name='inlet_iteration', index=inlets_current.index, dtype='object')
            self._inlets.loc[inlets_test, 'inlet_iteration'] = distribution

        return None
    
    def _construct_fmm_iterations(self):
        """
        TODO
        """
        # Set up iteration structure:
        iteration = 0
        self.iteration_outlets = pd.DataFrame()
        self.iteration_inlets  = pd.DataFrame()

        # Loops over outlet iterations
        for outlet_iteration in range(len(self.outlets_importance)):

            # Loops over inlet iterations
            for inlet_iteration in range(len(self.inlets_importance)):
                
                # Gets the outlets assigned to the current outlet iteration
                outlets_current = self._outlets[self._outlets['outlet_iteration'] == outlet_iteration]
                outlets_current = outlets_current.assign(iteration=iteration)
                self.iteration_outlets = pd.concat([self.iteration_outlets, outlets_current[['iteration', 'x', 'y', 'z']]])

                # Gets the inlets assigned to the current inlet iteration
                inlets_current = self._inlets[self._inlets['outlet_key'] == outlet_iteration]
                inlets_current_iteration = inlets_current[inlets_current['inlet_iteration'] == inlet_iteration]
                inlets_current_iteration = inlets_current_iteration.assign(iteration=iteration)
                self.iteration_inlets = pd.concat([self.iteration_inlets, inlets_current_iteration[['iteration', 'x', 'y', 'z']]])

                # Increments total iteration number by 1
                iteration = iteration + 1   
                
        return None
    
    def _repartition_points(self, nbr_points, importance):
        """
        Correct for integers in importance factors list not summing correctly to total number of points. # TODO ?
        """
        total_importance = float(sum(importance))                                               # total number of points as assigned (this may not be the actual total)
        proportion       = [float(i)/total_importance       for i in importance]                # percent of points to use this iteration
        repartition      = [round(proportion[i]*nbr_points) for i in range(len(proportion))]    # number of points to use this iteration (need to round bc percentage will not result in whole number)
        repartition[-1]  = nbr_points-sum(repartition[0:-1])                                    # leftover points are assignd to last iteration
        return repartition

    def _initialize_fmm_variables(self):
        """
        TODO
        """
        ### Raster maps
        self.maps = {
            'outlets' : np.full((self.grid.nx, self.grid.ny, self.grid.nz), np.nan), # map of null values where each cell with an outlet will have the index of that outlet
            'nodes'   : np.full((self.grid.nx, self.grid.ny, self.grid.nz), np.nan), # map of null values where each cell that has a node will be updated with that node index
            'cost'    : [], # cost of travel through each cell
            'alpha'   : [], # cost of travel along gradient through each cell
            'beta'    : [], # cost of travel perpendicular to gradient through each cell
            'time'    : [], # travel time to outlet from each cell
            'karst'   : [], # presence/absence of karst conduit in each cell
        }
        
        ### Defines outlets map according to outlets emplacements
        for (i, outlet) in self._outlets.iterrows(): #assign outlet indices. Compute the outlets map (array indicating location of outlets as their index and everywhere else as nan).
            X = self.grid.get_i(outlet.x)
            Y = self.grid.get_j(outlet.y)
            Z = self.grid.get_k(outlet.z)
            self.maps['outlets'][X][Y][Z] = i

        ### Vector maps
        self.vectors = {
            'nodes'     : {}, # empty dict to store nodes (key: nodeID, val: [x, y, z, type])
            'edges'     : {}, # empty dict to store edges (key: edgeID, val: [inNode, outNode])
            'n'         : 0,  # start node counter at zero
            'e'         : 0,  # start edge counter at zero
            'geodesics' : [], # empty list to store raw fast-marching path output
        }

        ### Set up fast-marching
        self.fmm = {
            'algorithm'     : self.SKS_SETTINGS['fmm']['algorithm'],
            'riemannMetric' : [],                                            # this changes at every iteration, but cannot be stored?
            'fastMarching'  : agd.Eikonal.dictIn({
                'model'             : self.SKS_SETTINGS['fmm']['algorithm'], # set algorithm from settings file ('Isotropic2', 'Isotropic3', 'Riemann2', 'Riemann3')
                'order'             : 2,                                     # recommended setting: 2 (replace by variable)
                'exportValues'      : 1,                                     # export the travel time field
                'exportGeodesicFlow': 1                                      # export the walker paths (i.e. the conduits)
            })
        }
        self.fmm['fastMarching'].SetRect(                                    # give the fast-marching algorithm the model grid
            sides=[[self.grid.xmin, self.grid.xmax],                         # bottom edge,   top edge (NOT centerpoint)
                [self.grid.ymin, self.grid.ymax],                            # leftmost edge, rightmost edge (NOT centerpoint)
                [self.grid.zmin, self.grid.zmax]],           
            dims=[self.grid.nx, self.grid.ny, self.grid.nz]                  # number of cells, number of cells, number of cells
        )
        return None


##########################################################################################################################################################################
### KARST NETWORK SIMULATION ###
################################

    # Computes karst network
    def compute_karst_network(self):
        """
        Compute the karst network according to the parameters.

        Append the results to the `karst_simulations` attribute as an element of an array.

        Parameters
        ----------
        TODO
        """
        ### Computes conduits for each generation & store nodes and edges for network
        self._compute_karst_network()

        ### ii. Calculates the karst network statistics indicators with karstnet and save karst network
        # karstnet_edges = list(self.edges.values()) #convert to format needed by karstnet (list)
        # karstnet_nodes = copy.deepcopy(self.nodes) #convert to format needed by karstnet (dic with only coordinates) - make sure to only modify a copy!
        # for key, value in karstnet_nodes.items():  #drop last item in list (the node type) for each dictionary entry
        #     value.pop()

        # Computes karstnet indicators
        # k = kn.KGraph(karstnet_edges, karstnet_nodes)  #make graph - edges must be a list, and nodes must be a dic of format {nodeindex: [x,y]}
        # stats = k.characterize_graph(verbose)

        ### iii. Store all the relevant data for this network in dictionaries:
        # maps = copy.deepcopy(self.maps)                  #store copy of maps
        # points = {}
        # points['inlets']  = copy.deepcopy(self._inlets)  #store copy of inlets in pandas df format with info about iteration and inlet/outlet assignment (this will be used in plotting later)
        # points['outlets'] = copy.deepcopy(self._outlets) #store copy of outlets in pandas df format with info about iteration and inlet/outlet assignment (this will be used in plotting later)
        # network = {}
        # network['edges'] = copy.deepcopy(self.edges) #store copy of edges list
        # network['nodes'] = copy.deepcopy(self.nodes) #store copy of nodes list
        # network['karstnet'] = copy.deepcopy(k)       #store copy of karstnet network object (including graph)
        # config = copy.deepcopy(self.SETTINGS)        #store copy of settings for the run being stored
        # try:
        #     del config['debug']
        # except:
        #     pass
        # self.SIMULATIONS.append(KarstNetwork(maps, points, network, stats, config))

        # TODO ??????????
        ### iv. - Return inlets and outlets to their original format:
        # #Chloe: this is to convert inlets and outlets back from pandas dataframes
        # #if the private self._inlets and self._outlets variables are working correctly, this step may be removed
        # self.inlets  = np.asarray(self._inlets)[:,0:2]    #return inlets to array format without info about iteration and inlet/outlet assignment (important for later iterations)
        # self.outlets = np.asarray(self._outlets)[:,0:2]   #return outlets to array format without info about iteration and inlet/outlet assignment (important for later iterations)
        return None
    
    
    def _compute_karst_network(self):
        """
        TODO
        """
        this.logger = logging.getLogger("fmm.modelisation")
        this.logger.info("Computing karst network")
        
        self.iteration = 0
        
        for self.iteration in range(self.nbr_iteration):
        
            # Compute travel time maps and conduit network
            if self.fmm['algorithm'] == 'Isotropic3':
                self._compute_cost_map()    # 2.1.1
                self._compute_time_map()    # 2.1.2
                self._compute_karst_map()   # 2.1.3

            elif self.fmm['algorithm'] == 'Riemann3':
                # TODO - check notebooks mirebeau
                self._compute_cost_map()        # 2.1.3
                self._compute_alpha_map()       # 2.1.4
                self._compute_beta_map()        # 2.1.5
                self._compute_riemann_metric()  # 2.1.6
                self._compute_time_map()        # 2.1.7
                self._compute_karst_map()       # 2.1.3

        # self.iteration = self.iteration + 1   #increment total iteration number by 1
            this.logger.info("iteration : {}/{}".format(self.iteration+1, self.nbr_iteration))
        return None


    ### 2.1.1 Iso- and anisotropic case
    @wp._debug_level(1, True)
    @wp._logging()
    def _compute_cost_map(self):
        """
        Compute the cost map (how difficult it is to traverse each cell).

        TODO
        """
        # If it's the first iteration, iniatialize the cost map according to the geological settings.
        if self.iteration == 0:

            cost_map_0 = np.zeros((self.grid.nx, self.grid.ny, self.grid.nz))
            self.maps['cost'].append(cost_map_0)
            
            ### Geology
            for key, cost in self.geology.cost.items():
                cost_map_0 = np.where(self.geology.data == key, cost, cost_map_0)
                self.maps['cost'][0] = cost_map_0

            ### Faults
            if self.faults is not None:
                for key, cost in self.faults.cost.items():
                    cost_map_0 = np.where(self.faults.data == key, cost, cost_map_0)
                    self.maps['cost'][0] = cost_map_0

            ### TODO
            ### Karst
            # if self.karst is not None:
            # for key, cost in self.karst.cost.items():
            #     self.maps['cost'][0] = np.where(self.karst.data == key, cost, self.maps['cost'][0])

            ### TODO
            ### Fractures
            # if self.fractures is not None:
            # for key, cost in self.fractures.cost.items():
                # self.maps['cost'][0] = np.where(self.fractures.data == key, cost, self.maps['cost'][0])

            ### If out of mask
            if self.mask is not None:
                cost_map_0 = np.where(self.mask.mask == 0, this.cost['out'], cost_map_0)
                self.maps['cost'][0] = cost_map_0

        # If it's not the first iteration
        else: 
            self.maps['cost'].append(self.maps['cost'][self.iteration-1])
            self.maps['cost'][self.iteration] = np.where(self.maps['karst'][self.iteration-1] > 0, this.cost['conduits'], self.maps['cost'][self.iteration])

        return None


    ### 2.1.4 Anisotropic case
    # @wp._decorator_logging_fmm()
    @wp._debug_level(2, True)
    @wp._logging()
    def _compute_alpha_map(self):
        """
        Compute the alpha map: travel cost in the same direction as the gradient.
        Cost map * topography map, so that the cost is higher at higher elevations, encouraging conduits to go downgradient.

        TODO : à terminer
        """
        alpha_map = self.maps['cost'][self.iteration] * self.grid.Z
        self.maps['alpha'].append(alpha_map)

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
    @wp._debug_level(3, True)
    @wp._logging()
    def _compute_beta_map(self):
        """
        Compute the beta map: travel cost perpendicular to the gradient.
        If beta is higher than alpha, conduits will follow the steepest gradient.
        If beta is lower than alpha, conduits will follow contours.
        """
        beta_map = self.maps['alpha'][self.iteration] / self.SKS_SETTINGS['fmm']['cost']['ratio']
        self.maps['beta'].append(beta_map)

        # TODO - alpha = beta dans zone saturée
        return None


    ### 2.1.6 Anisotropic case
    @wp._debug_level(4, True)
    @wp._logging()
    def _compute_riemann_metric(self):
        """
        Compute the riemann metric: Define the Riemannian metric needed as input for the anisotropic fast marching.

        TODO : à terminer
        """
        x, y, z = self.grid.Z.shape
        # orientationx = np.transpose(np.gradient(self.grid.Z, axis=0), (2,1,0))
        # orientationy = np.transpose(np.gradient(self.grid.Z, axis=1), (2,1,0))
        # orientationz = np.transpose(np.gradient(self.grid.Z, axis=2), (2,1,0))
        orientationx = np.gradient(self.grid.Z, axis=0)
        orientationy = np.gradient(self.grid.Z, axis=1)
        orientationz = np.gradient(self.grid.Z, axis=2)

        # TODO
        # if (z == 1):
        #     orientationz = np.transpose(self.grid.Z, (2,1,0)) # TODO ??!!
        # else:
        #     orientationz = np.transpose(np.gradient(self.grid.Z, axis=2), (2,1,0))

        # alpha = np.transpose(self.maps['alpha'][iteration], (2,1,0))
        # beta  = np.transpose(self.maps['beta'][iteration],  (2,1,0))
        alpha = self.maps['alpha'][self.iteration]
        beta  = self.maps['beta'][self.iteration]
        self.fmm['riemannMetric'] = agd.Metrics.Riemann.needle([orientationx, orientationy, orientationz], alpha, beta)  
        return None


    ### 2.1.2
    @wp._debug_level(5, True)
    @wp._logging()
    def _compute_time_map(self):
        """
        Compute the travel time map (how long it takes to get to the outlet from each cell),
        using the ani- or isotropic agd-hfm fast-marching algorithm, and store travel time map.
        Note: the AGD-HFM library uses different indexing, so x and y indices are reversed for inlets and outlets.
        TODO
        """
        # Set the outlets for this iteration
        seeds_x = self.iteration_outlets[self.iteration_outlets['iteration']==self.iteration].x
        seeds_y = self.iteration_outlets[self.iteration_outlets['iteration']==self.iteration].y
        seeds_z = self.iteration_outlets[self.iteration_outlets['iteration']==self.iteration].z
        seeds   = list(zip(seeds_x, seeds_y, seeds_z))
        self.fmm['fastMarching']['seeds'] = seeds

        # Select inlets for current iteration
        tips_x = self.iteration_inlets[self.iteration_inlets['iteration']==self.iteration].x
        tips_y = self.iteration_inlets[self.iteration_inlets['iteration']==self.iteration].y
        tips_z = self.iteration_inlets[self.iteration_inlets['iteration']==self.iteration].z
        tips   = list(zip(tips_x, tips_y, tips_z))
        self.fmm['fastMarching']['tips'] = tips

        # Set the travel cost through each cell
        if self.fmm['algorithm'] == 'Isotropic3':
            self.fmm['fastMarching']['cost'] = self.maps['cost'][self.iteration]
        elif self.fmm['algorithm'] == 'Riemann3':
            self.fmm['fastMarching']['metric'] = self.fmm['riemannMetric']

        # Set verbosity of hfm run
        self.fmm['fastMarching']['verbosity'] = self.SKS_SETTINGS['verbosity']['agd']

        # Run the fast marching algorithm and store the outputs
        self.fmm['fastMarchingOutput'] = self.fmm['fastMarching'].Run()

        # Store travel time maps
        self.maps['time'].append(self.fmm['fastMarchingOutput']['values'])

        # Store fastest travel paths
        self.vectors['geodesics'].append(self.fmm['fastMarchingOutput']['geodesics'])
        return None


    ### 2.1.3 
    @wp._debug_level(6, True)
    @wp._logging()
    def _compute_karst_map(self):
        """
        Compute the karst map based on the paths from agd-hfm.
        Array of all zeros, with ones in cells containing a karst conduit.
        """
        # Get karst map from previous iteration (except for the very first iteration)
        if self.iteration == 0:
            karst_map_0 = np.zeros((self.grid.nx, self.grid.ny, self.grid.nz))
            self.maps['karst'].append(karst_map_0)
        else:
            self.maps['karst'].append(self.maps['karst'][self.iteration-1])

        # Debugging plot:
        # Chloe: this should stay in since it is very useful if there are problems
        # f1, ax1 = plt.subplots(1, 1, figsize=(10, 10))

        ### Loop over conduit paths generated by fast marching:
        for path in self.fmm['fastMarchingOutput']['geodesics']: #loop over conduit paths in this iteration (there is one path from each inlet)
            merge = False                                        #reset indicator for whether this conduit has merged with an existing conduit
            for p in range(path.shape[1]):                       #loop over points making up this conduit path
                point = path[:,p]                                #get coordinates of current point
                [[ix, iy, iz], error] = self.fmm['fastMarching'].IndexFromPoint(point) #convert to coordinates to indices, /!\ returning iy first then ix TODO ???
                # ax1.scatter(point[1], point[0], c='g',s=5)  #debugging
                #Place nodes and links:
                if np.isnan(self.maps['nodes'][ix, iy, iz]):                                    #if there is no existing conduit node here
                    if ~np.isnan(self.maps['outlets'][ix, iy, iz]):                              #if there is an outlet here (cell value is not nan)
                        outlet = self._outlets.iloc[int(self.maps['outlets'][ix, iy, iz])]         #get the outlet coordinates using the ID in the outlets map
                        self.vectors['nodes'][self.vectors['n']] = [outlet.x, outlet.y, outlet.z, 'outfall']           #add a node at the outlet coordinates (with the node type for SWMM)
                        self.maps['nodes'][ix, iy, iz] = self.vectors['n']                                   #update node map with node index
                        # ax1.scatter(outlet.x, outlet.y, marker='o', c='b')                   #debugging
                        if p > 0:                                                           #if this is not the first point (i.e. the inlet) in the current path
                            if merge == False:                                               #if this conduit has not merged with an existing conduit
                                self.vectors['edges'][self.vectors['e']] = [self.vectors['n']-1, self.vectors['n']]                       #add an edge connecting the previous node to the current node
                                self.vectors['e'] = self.vectors['e'] + 1                                             #increment edge counter up by one
                                # ax1.plot((self.vectors['nodes'][self.vectors['n']][0], self.vectors['nodes'][self.vectors['n']-1][0]), (self.vectors['nodes'][self.vectors['n']][1], self.vectors['nodes'][self.vectors['n']-1][1]))
                            else:                                                          #if this conduit HAS merged with an existing conduit
                                [[fromix, fromiy, fromiz], error] = self.fmm['fastMarching'].IndexFromPoint(path[:, p-1]) #get xyz indices of previous point in current conduit path
                                n_from = self.maps['nodes'][fromix, fromiy, fromiz]           #get node index of the node already in the cell where the previous point was
                                self.vectors['edges'][self.vectors['e']] = [n_from, self.vectors['n']]                         #add an edge connecting existing conduit node to current node
                                self.vectors['e'] = self.vectors['e'] + 1                                             #increment edge counter up by one
                                # ax1.plot((self.vectors['nodes'][self.vectors['n']].x, self.vectors['nodes'][n_from].x), (self.vectors['nodes'][self.vectors['n']].y, self.vectors['nodes'][n_from].y))
                        self.vectors['n'] = self.vectors['n'] + 1                                                   #increment node counter up by one
                    else:                                                                  #if there is NOT an outlet here
                        if p > 0:                                                           #if this is not the first point in the current path
                            #possible improvement: if the next point on the path is on an existing point, skip the current point.
                            self.vectors['nodes'][self.vectors['n']] = [point[0], point[1], point[2], 'junction']            #add a junction node here (with the node type for SWMM)
                            self.maps['nodes'][ix, iy, iz] = self.vectors['n']                               #update node map with node index
                            #ax1.scatter(point[1],point[0], marker='o', edgecolor='g', facecolor='none')   #debugging
                            if merge == False:                                              #if this conduit has not merged with an existing conduit
                                self.vectors['edges'][self.vectors['e']] = [self.vectors['n']-1, self.vectors['n']]                      #add and edge connecting the previous node to the current node
                                self.vectors['e'] = self.vectors['e'] + 1                                            #increment edge counter up by one
                                #ax1.plot((self.vectors['nodes'][self.vectors['n']][1], self.vectors['nodes'][self.vectors['n']-1][1]),(self.vectors['nodes'][self.vectors['n']][0], self.vectors['nodes'][self.vectors['n']-1][0]), c='gold', marker=None)
                            else:                                                           #if this conduit HAS merged with an existing conduit
                                [[fromix, fromiy, fromiz], error] = self.fmm['fastMarching'].IndexFromPoint(path[:, p-1]) #get xy indices of previous point in current conduit path
                                n_from = self.maps['nodes'][fromix, fromiy, fromiz]                   #get node index of the node already in the cell where the previous point was
                                self.vectors['edges'][self.vectors['e']] = [n_from, self.vectors['n']]                        #add an edge connecting existing conduit node to current node
                                self.vectors['e'] = self.vectors['e'] + 1                                            #increment edge counter up by one
                                merge = False                                                #reset merge indicator to show that current conduit has left                                                              #if this is the first point in current path
                        else:                                                                #if this is the first point in the current path (counter <= 0, therefore it is an inlet)
                            self.vectors['nodes'][self.vectors['n']] = [point[0], point[1], point[2], 'inlet']               #add an inlet node here (with the node type for SWMM)
                            self.maps['nodes'][ix, iy, iz] = self.vectors['n']                               #update node map with node index
                            #ax1.scatter(point[1],point[0], marker='o', edgecolor='sienna', facecolor='none')
                        self.vectors['n'] = self.vectors['n'] + 1                                                   #increment node counter up by one
                elif ~np.isnan(self.maps['nodes'][ix, iy, iz]):                                 #if there is already a node in this cell (either because there is a conduit here, or because there are two nodes in the same cell)
                    n_existing = self.maps['nodes'][ix, iy, iz]                                  #get index of node already present in current cell
                    if merge == True:                                                       #if this conduit has already merged into an existing conduit
                        pass                                                                 #skip this node (there is already a node here)
                    elif n_existing == self.vectors['n']-1:                                            #if existing index is only one less than next node to be added index, this is a duplicate node and can be skipped
                        pass                                                                 #skip this node (duplicate)
                    else:                                                                   #if existing node index is >1 less than next node to be added index
                        if p > 0:                                                           #if this is not the first point in the current path
                            self.vectors['edges'][self.vectors['e']] = [self.vectors['n']-1, n_existing]                      #add an edge connecting most recently added node and existing node in cell
                            self.vectors['e'] = self.vectors['e'] + 1                                                #increment edge counter up by one
                            merge = True                                                     #add a flag indicating that this conduit has merged into an existing one
                            #ax1.plot((self.vectors['nodes'][self.vectors['n']-1][1], self.vectors['nodes'][n_existing][1]),(self.vectors['nodes'][self.vectors['n']-1][0], self.vectors['nodes'][n_existing][0]), c='r', marker=None)
                        else:                                                                #if this is the first point in the current path (i.e. the inlet is on an exising conduit)
                            self.vectors['nodes'][self.vectors['n']] = [point[0], point[1], point[2], 'inlet']                #add a node here (with the node type for SWMM)- this will cause there to be two nodes in the same cell
                            self.maps['nodes'][ix, iy, iz] = self.vectors['n']                                #update node map with node index
                            self.vectors['n'] = self.vectors['n'] + 1                                                 #increment node counter by 1
                            #ax1.scatter(point[1],point[0], marker='o', edgecolor='g', facecolor='none')  #debugging
                self.maps['karst'][self.iteration][ix, iy, iz] = 1                               #update karst map to put a conduit in current cell
                

        ### Debugging plot
        # # Display inlets and outlets
        # ax1.scatter(self._outlets.x, self._outlets.y, c='cyan',   s=100)
        # ax1.scatter(self._inlets.x,  self._inlets.y,  c='orange', s=100)
        # ax1.scatter(self._outlets[self._outlets.iteration==iteration].x, self._outlets[self._outlets.iteration==iteration].y, c='cyan',   s=100)
        # ax1.scatter(self._inlets[self._inlets.iteration==iteration].x,   self._inlets[self._inlets.iteration==iteration].y,   c='orange', s=100)
        # # Display karst network
        # ax1.imshow(np.transpose(self.maps['karst'][iteration], (1,0,2)), origin='lower', extent=self.grid.extent, cmap='gray_r')
        # ax1.imshow(np.transpose(self.maps['nodes'], (1,0,2)), origin='lower', extent=self.grid.extent, cmap='gray_r')
        return None