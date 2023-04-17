"""
This module contains a class modeling the karstic network generator tool.
"""

### Internal dependencies
import os
import sys
import copy
import pickle
import shutil
import logging
import datetime

### External dependencies
import yaml
import numpy as np
import pandas as pd

### Fast-Marching package
import agd
from agd import Eikonal
from agd.Metrics import Riemann

### Local dependencies
import pykasso.core._wrappers as wp
from .grid import Grid
from .domain import Domain, Delimitation, Topography, Bedrock, WaterLevel
from .geologic_features import Geology, Faults
from .fracturation import Fractures
from .points import PointGenerator

### Typing
from pykasso._typing import DataFrame

### Module variables
this = sys.modules[__name__]

# Dictionary used for memory optimization/memoizing operations
this.ACTIVE_PROJECT = None

# Defines default fast-marching costs
this.default_fmm_costs = {
    # 'out': 0.999,
    'out': 10,
    'aquifer': 0.4,
    'aquiclude': 0.8,
    'beddings': 0.35,
    'faults': 0.2,
    'fractures': 0.2,
    'karst': 0.1,
    'conduits': 0.1,
    'ratio': 0.5,
}

# Defines misc directory path
this.misc_relative_path = '/../_misc/'

########################
### Module functions ###
########################


def create_project(project_directory: str, example: str = None) -> None:
    """
    Creates a new pyKasso's project within the provided `project_directory`
    path.
    Three subdirectories will be created :
     - ``inputs/``, for input model data files
     - ``outputs/``, for model results storage
     - ``settings/``, for pyKasso's settings file
    
    .. warning:: pyKasso doesn't check if the directory already exists and if
    yes, its actual content.

    Parameters
    ----------
    project_directory : str
        Location, path of the new project directory.
    example : str
        __doc__ :
        -"betteraz"
        
    Examples
    --------
    >>> import pykasso as pk
    >>> pk.create_project('examples/betteraz/')
    """
    # Uncomment this to tell apart the development version and the main version
    # print('CAUTION: You are using the development version of this package.')
    
    ### Creates the project subdirectories
    this.core_settings = {
        'project_directory': project_directory,
        'inputs_directory': project_directory + '/inputs/',
        'outputs_directory': project_directory + '/outputs/',
        'settings_directory': project_directory + '/settings/',
    }
    [os.makedirs(directory, exist_ok=True) for directory in core_settings.values()]

    ### Populates the settings directory with initial yaml settings file
    location = os.path.dirname(os.path.abspath(__file__))
    misc_directory = location + this.misc_relative_path
    this.core_settings.update({
        'core_settings_filename': 'CORE_SETTINGS.yaml',
        'sks_settings_filename': 'SKS_SETTINGS.yaml',
        'log_filename': 'pyKasso.log',
    })
    
    # Copies directly the requested example from the _misc directory
    if example is not None:
        example = example.strip('/')
        # path = 
        return None

    ### Creates the CORE_SETTINGS.yaml file in the project directory
    # core_settings_filename = (core_settings['settings_directory']
    #                           + core_settings['core_settings_filename'])
    # with open(core_settings_filename, 'w') as f:
    #     text = "# pyKasso CORE SETTINGS \n---\n"
    #     yaml_text = yaml.safe_dump(core_settings)
    #     text = text + yaml_text + "..."
    #     f.write(text)

    ### Copies the default SKS_SETTINGS.yaml file from the misc directory to
    ### the project directory
    path = misc_directory + this.core_settings['sks_settings_filename']
    shutil.copy2(path, project_directory + '/settings/')

    ### Constructs the dictionary used for settings comparison between actual
    ### and previous simulation, used for the memoization operation
    this.ACTIVE_PROJECT = this.ACTIVE_PROJECT = {
        'project_directory': project_directory,
        'n_simulation': 0,
        'simulation_locations': [],
        'settings': {
            'grid': None,
            'domain': None,
            'geology': None,
            'beddings': None,
            'faults': None,
            # static features
            'fractures': None,
            'outlets': None,
            'inlets': None,
            # dynamic features
        },
        'model': {}
    }
    return None


def save_project(filename: str) -> None:
    """
    Exports a pyKasso's project in a python pickle and saves it to the
    provided `filename`.

    Parameters
    ----------
    filename : str
        Filename to which the project is saved.
        A ``.pickle`` extension will be appended to the filename if it does
        not already have one.
        
    Examples
    --------
    >>> import pykasso as pk
    >>> pk.create_project('examples/betteraz/')
    >>> pk.save_project('examples/betteraz/project.pickle')
    
    Notes
    -----
    For a description of the ``.pickle`` format, see https://docs.python.org/3/library/pickle.html.
    """
    # Controls filename extension
    from numpy.compat import os_fspath
    filename = os_fspath(filename)
    if not filename.endswith('.pickle'):
        filename = filename + '.pickle'
    
    # Saves the project in a pickle
    with open(filename, 'wb') as handle:
        pickle.dump(this.ACTIVE_PROJECT,
                    handle,
                    protocol=pickle.HIGHEST_PROTOCOL)
        
    return None


def load_project(filename: str) -> None:
    """
    Loads a pyKasso's project from a python pickle file located to the
    provided `filename`.
    
    .. warning:: Loading files that contain object arrays uses the `pickle`
                 module, which is not secure against erroneous or maliciously
                 constructed data.

    Parameters
    ----------
    filename : str
        Filename to which the project is stored.

    Examples
    --------
    >>> import pykasso as pk
    >>> pk.load_project('examples/betteraz/project.pickle')
    """
    # Loads the pickle
    with open(filename, 'rb') as handle:
        this.ACTIVE_PROJECT = pickle.load(handle)
        
    return None


#################
### SKS Class ###
#################

class SKS():
    """
    Class storing all the parameters, the data and the resulting simulated
    stochastic karst networks.
    
    Attributes
    ----------
    SKS_SETTINGS : dict
        ...
    grid : _type_
        ...
    domain : _type_
        ...
    geology : _type_
        ...
    faults : _type_
        ...
    fractures : _type_
        ...
    inlets : DataFrame
        ...
    outlets : DataFrame
        ...
    
    .. note:: ``faults`` or ``fractures`` will return None when not defined.
    """

    def __init__(self, sks_settings: dict = {}, export_settings: dict = {},
                 debug_level: dict = {}) -> None:
        """Initializes an instance of SKS class.

        Parameters
        ----------
        sks_settings : dict, optional
            _description_, by default {}
        export_settings : dict, optional
            _description_, by default {}
        debug_level : dict, optional
            _description_, by default {}
            
        Examples
        --------
        >>> import pykasso as pk
        >>> simulation = pk.SKS()
        """
        ### Initialization
        self.grid = None
        self.domain = None
        self.geology = None
        self.faults = None
        self.fractures = None
        self.inlets = None
        self.outlets = None
        self.conceptual_model = None
        self.conceptual_model_table = None
        self.rng = {}
        self.orientation = None  # TODO
        
        ### Loads core and sks settings
        settings_list = [copy.deepcopy(this.core_settings),
                         copy.deepcopy(sks_settings)]
        settings_names = ['CORE_SETTINGS', 'SKS_SETTINGS']
        for settings, settings_name in zip(settings_list, settings_names):
            if settings == {}:
                settings = this.ACTIVE_PROJECT['project_directory'] + '/settings/{}.yaml'.format(settings_name)
            if isinstance(settings, str):
                with open(settings, 'r') as f:
                    setattr(self, settings_name.upper(), yaml.safe_load(f))
            if isinstance(settings, dict):
                setattr(self, settings_name.upper(), settings)
        
        ### Creates simulation directory
        this.ACTIVE_PROJECT['n_simulation'] += 1
        self.CORE_SETTINGS['simulation_directory'] = self.CORE_SETTINGS['outputs_directory'] + 'simulation_{}/'.format(this.ACTIVE_PROJECT['n_simulation'])
        os.makedirs(self.CORE_SETTINGS['simulation_directory'], exist_ok=True)
        this.ACTIVE_PROJECT['simulation_locations'].append(self.CORE_SETTINGS['simulation_directory'])
        
        ### Checks if verbosity levels have been declared
        if 'verbosity' not in self.SKS_SETTINGS:
            self.SKS_SETTINGS['verbosity'] = {
                'logging': 0,
                'agd': 0,
            }
        else:
            if 'logging' not in self.SKS_SETTINGS['verbosity']:
                self.SKS_SETTINGS['logging'] = 0
            if 'agd' not in self.SKS_SETTINGS['verbosity']:
                self.SKS_SETTINGS['agd'] = 0
        
        ### LOGGING ###
        
        ### Sets logging output formats
        # logging_title = ' %(message)s'
        logging_entry = ' %(name)-30s | %(levelname)-8s | %(message)s'

        ### Sets logging level
        levels = {
            0: logging.DEBUG,
            1: logging.INFO,
            2: logging.WARNING,
            3: logging.ERROR,
            4: logging.CRITICAL,
        }
        level = self.SKS_SETTINGS['verbosity']['logging']
        if level > 4:
            level = 4
        if level < 0:
            level = 0
        logging_level = levels[level]
        
        ### Sets logging file
        log_file = self.CORE_SETTINGS['outputs_directory'] + self.CORE_SETTINGS['log_filename']

        ### Sets new logging file
        if this.ACTIVE_PROJECT['n_simulation'] == 1:
            
            # Resets logging module
            root = logging.getLogger()
            list(map(root.removeHandler, root.handlers))
            list(map(root.removeFilter, root.filters))
            
            # Removes PIL package logging entries
            logging.getLogger("PIL.PngImagePlugin").setLevel(logging.CRITICAL + 1)
            
            # Creates new logging file
            logging.basicConfig(
                filename=log_file,
                encoding='utf-8',
                level=logging_level,
                filemode="w",
                format=logging_entry
            )
            
            # Prints pyKasso 'logo'
            this.logger = logging.getLogger("♠")
            this.logger.info("             _                        ")
            this.logger.info("            | |                       ")
            this.logger.info(" _ __  _   _| | ____ _ ___ ___  ___   ")
            this.logger.info("| `_ \| | | | |/ / _` / __/ __|/ _ \  ")
            this.logger.info("| |_) | |_| |   < (_| \__ \__ \ (_) | ")
            this.logger.info("| .__/ \__, |_|\_\__,_|___/___/\___/  ")
            this.logger.info("| |     __/ |                         ")
            this.logger.info("|_|    |___/                          ")
            this.logger.info("                                      ")
      
            # Prints basic information in the log
            this.logger = logging.getLogger("☼")
            this.logger.info("Project directory location : " + this.ACTIVE_PROJECT['project_directory'])
            this.logger.info("Project creation date      : " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
            # TODO - error with binder
            try:
                # pykasso_version = __all__.__version__ # TODO
                pykasso_version = 'TODO'
                this.logger.info("pyKasso version used       : " + pykasso_version)
            except:
                pass
            this.logger.info("Logging level              : " + str(logging_level)) # TODO ????????
            
        ### Prints current simulation number
        logging.basicConfig(
            filename=log_file,
            encoding='utf-8',
            level=logging_level,
            filemode="w",
            format=logging_entry
        )
        this.logger = logging.getLogger("►")
        l = len(str(this.ACTIVE_PROJECT['n_simulation']))
        this.logger.info('***********************{}****'.format('*' * l))
        this.logger.info('*** pyKasso simulation {} ***'.format(this.ACTIVE_PROJECT['n_simulation']))
        this.logger.info('***********************{}****'.format('*' * l))
        this.logger.info('Location : ' + self.CORE_SETTINGS['simulation_directory'])
        this.logger.info('Date     : ' + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        this.logger.info('---')

        ### Sets debug level
        if debug_level == {}:
            self.debug_mode = False
            self.debug_level = {
                'model': 20,
                'simulation': 10,
                'iteration': 1000  # TODO
            }
        else:
            self.debug_mode = True
            self.debug_level = debug_level
            if 'model' not in self.debug_level:
                self.debug_level['model'] = 20
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
        
##########################
### BUILDING THE MODEL ###
##########################

    def build(self) -> None:
        """
        Builds the geological model. Should be called right after ``SKS``
        class instantiation.
        
        Model building procedure:
         1. Constructs the grid : ``self.grid`` ;
         2. Defines the model extension : ``self.domain`` ;
         3. Builds the geological features : ``self.geology``, ``self.faults``,
         ``self.inlets``, ``self.outlets``, ``self.fractures``,
         ``self.conceptual_model``, ``self.conceptual_model_table``;
         4. Prepares the fast-marching method variables.
         
        If a geological feature has not been defined in the settings, calling
        its attribute will returns ``None``.
        
        Examples
        --------
        >>> import pykasso as pk
        >>> simulation = pk.SKS()
        >>> simulation.build()
        """
        self._build_grid()               # 1 - Constructs the grid
        self._build_domain()             # 2 - Constructs the domain
        self._build_model()              # 3 - Constructs the model
        self._construct_fmm_variables()  # 4 - Constructs fmm features
        return None
    
    ##### 1 - GRID #########################
    @wp._debug_level(1)
    @wp._parameters_validation('grid', 'required')
    @wp._memoize('grid')
    @wp._logging()
    def _build_grid(self) -> None:
        """Builds the grid."""
        if isinstance(self.SKS_SETTINGS['grid'], (str)):
            # TODO builds the grid according to a .asc or .grd file
            pass
        else:
            self.grid = Grid(**self.SKS_SETTINGS['grid'])
        return None
    
    ##### 2 - DOMAIN #########################
    @wp._debug_level(2)
    @wp._parameters_validation('domain', 'optional')
    @wp._memoize('domain')
    @wp._logging()
    def _build_domain(self) -> None:
        """Builds the domain."""
        delimitation = self._build_domain_delimitation()
        topography = self._build_domain_topography()
        bedrock = self._build_domain_bedrock()
        water_level = self._build_domain_water_level()
        self.domain = Domain(self.grid,
                             delimitation=delimitation,
                             topography=topography,
                             bedrock=bedrock,
                             water_level=water_level)
        return None
    
    @wp._debug_level(2.1)
    @wp._logging()
    def _build_domain_delimitation(self):
        """Builds the delimitation."""
        delimitation = self.SKS_SETTINGS['domain']['delimitation']
        test_a = isinstance(delimitation, (str))
        test_b = (delimitation == '')
        if not (test_a and test_b):
            if isinstance(delimitation, (str)):
                delimitation = np.genfromtxt(delimitation).tolist()
                self.SKS_SETTINGS['domain']['delimitation'] = delimitation
            instance = Delimitation(vertices=delimitation, grid=self.grid)
            return instance
        else:
            return None

    @wp._debug_level(2.2)
    @wp._logging()
    def _build_domain_topography(self):
        """Builds the topography."""
        topography = self.SKS_SETTINGS['domain']['topography']
        test_a = isinstance(topography, (str))
        test_b = (topography == '')
        if not (test_a and test_b):
            instance = Topography()
            instance.set_data(data=topography, grid=self.grid)
            instance._surface_to_volume('<=', self.grid)
            return instance
        else:
            return None
        
    @wp._debug_level(2.3)
    @wp._logging()
    def _build_domain_bedrock(self):
        """Builds the bedrock elevation."""
        bedrock = self.SKS_SETTINGS['domain']['bedrock']
        test_a = isinstance(bedrock, (str))
        test_b = (bedrock == '')
        if not (test_a and test_b):
            instance = Bedrock()
            instance.set_data(data=bedrock, grid=self.grid)
            instance._surface_to_volume('<=', self.grid)
            return instance
        else:
            return None
        
    @wp._debug_level(2.4)
    @wp._logging()
    def _build_domain_water_level(self):
        """Builds the water level elevation."""
        water_level = self.SKS_SETTINGS['domain']['water_level']
        test_a = isinstance(water_level, (str))
        test_b = (water_level == '')
        if not (test_a and test_b):
            instance = WaterLevel()
            instance.set_data(data=water_level, grid=self.grid)
            instance._surface_to_volume('<=', self.grid)
            return instance
        else:
            return None

    ##### 3 - MODEL #########################
    @wp._debug_level(3)
    @wp._logging()
    def _build_model(self) -> None:
        """Builds the geological features."""
        
        ### Sets the main parameters
        self._build_model_parameters()
        
        ### Constructs the deterministic geological features
        self._build_model_geology()
        self._build_model_faults()

        ### Constructs the stochastic geological features
        self._build_model_outlets()
        self._build_model_inlets()
        # TODO - self._build_model_tracers()
        self._build_model_fractures()
        
        ### Constructs the conceptual model
        self._build_conceptual_model()
        return None
    
    def _set_rng(self, attribute):
        """
        Populates the rng seeds dictionary with the defined seed attribute.
        """
        
        # Generates a random seed when the seed equals 0 (default setting)
        if self.SKS_SETTINGS[attribute]['seed'] == 0:
            seed = np.random.default_rng().integers(low=0, high=10**6)
            self.SKS_SETTINGS[attribute]['seed'] = seed
            
        # Sets the seed
        seed = self.SKS_SETTINGS[attribute]['seed']
        self.rng[attribute] = np.random.default_rng(seed)
        return None
    
    @wp._debug_level(3.1)
    @wp._parameters_validation('sks', 'optional')
    def _build_model_parameters(self) -> None:
        """Builds characteristic model parameters."""
        # Defines the main seed
        self._set_rng('sks')
        return None
    
    @wp._debug_level(3.2)
    @wp._parameters_validation('geology', 'optional')
    @wp._memoize('geology')
    @wp._logging()
    def _build_model_geology(self) -> None:
        """Builds the geology."""
        geology = self.SKS_SETTINGS['geology']['data']
        test_a = isinstance(geology, (str))
        test_b = (geology == '')
        self.geology = Geology()
        if not (test_a and test_b):
            axis = self.SKS_SETTINGS['geology']['axis']
            self.geology.set_data(data=geology, grid=self.grid, axis=axis)
        else:
            self.geology.data_volume = (
                self.geology._set_data_full_3D(grid=self.grid, value=1)
            )
        costs = self.SKS_SETTINGS['geology']['costs']
        self.geology._set_costs(costs)
        self.geology._compute_statistics(self.grid)
        return None
    
    @wp._debug_level(3.3)
    @wp._parameters_validation('faults', 'optional')
    @wp._memoize('faults')
    @wp._logging()
    def _build_model_faults(self) -> None:
        """Builds the faults."""
        faults = self.SKS_SETTINGS['faults']['data']
        test_a = isinstance(faults, (str))
        test_b = (faults == '')
        if not (test_a and test_b):
            self.faults = Faults()
            axis = self.SKS_SETTINGS['faults']['axis']
            self.faults.set_data(data=faults, grid=self.grid, axis=axis)
            costs = self.SKS_SETTINGS['faults']['costs']
            self.faults._set_costs(costs)
            self.faults._compute_statistics(self.grid)
        return None

    @wp._debug_level(3.5)
    @wp._parameters_validation('outlets', 'required')
    @wp._memoize('outlets')
    @wp._logging()
    def _build_model_outlets(self) -> None:
        """Builds the outlets."""
        self._set_rng('outlets')
        self.outlets = self._construct_feature_points('outlets')
        if self.SKS_SETTINGS['outlets']['shuffle']:
            pts = self.outlets.sample(frac=1, random_state=self.rng['sks'])
            self.outlets = pts.reset_index(drop=True)
        return None

    @wp._debug_level(3.6)
    @wp._parameters_validation('inlets', 'required')
    @wp._memoize('inlets')
    @wp._logging()
    def _build_model_inlets(self) -> None:
        """Builds the inlets."""
        self._set_rng('inlets')
        self.inlets = self._construct_feature_points('inlets')
        if self.SKS_SETTINGS['inlets']['shuffle']:
            pts = self.inlets.sample(frac=1, random_state=self.rng['sks'])
            self.inlets = pts.reset_index(drop=True)
        return None
    
    def _construct_feature_points(self, kind):
        """Constructs the inlets / outlets.

        Four cases possible:
         1. No points declared
         2. More points required than provided : Generates additional random
         points
         3. Less points required than provided : Pick random points among
         provided ones
         4. Points required equals points declared
        """
        ### Creates a point generator instance
        point_manager = PointGenerator(
            rng=self.rng[kind],
            mode=self.SKS_SETTINGS[kind]['mode'],
            domain=self.domain,
            geology=self.geology,
            geologic_ids=self.SKS_SETTINGS[kind]['geology']
        )

        ### Gets existing points

        # Loads points if needed
        logical_test_01 = isinstance(self.SKS_SETTINGS[kind]['data'], (str))
        logical_test_02 = not (self.SKS_SETTINGS[kind]['data'] == '')
        if logical_test_01 and logical_test_02:
            path = self.SKS_SETTINGS[kind]['data']
            self.SKS_SETTINGS[kind]['data'] = np.genfromtxt(path)
        
        ### Inspects validity of points
        points = self.SKS_SETTINGS[kind]['data']
        
        # 2D points # TODO - logging ?
        points_2D = [point for point in points if len(point) == 2]
        validated_points_2D = [point for point in points_2D
                               if point_manager._is_coordinate_2D_valid(point)]
        validated_points_3D = [
            point_manager._generate_3D_coord_from_2D_coord(point)
            for point in validated_points_2D
        ]
        
        # 3D points # TODO - logging ?
        points_3D = [point for point in points if len(point) == 3]
        validated_points_3D += [point for point in points_3D if point_manager._is_coordinate_3D_valid(point)]

        diff = len(points) - len(validated_points_3D)
        if diff > 0:
            # TODO - LOG - VERBOSITY
            # TODO - log name is not correct
            msg = '{}/{} {} have been discarded because out of domain.'.format(diff, len(self.SKS_SETTINGS[kind]['data']), kind)
            this.logger.warning(msg)
        self.SKS_SETTINGS[kind]['data'] = validated_points_3D

        ### Get new points according to the right case
        # Case 1 - No points declared
        if (self.SKS_SETTINGS[kind]['data'] == '') or (self.SKS_SETTINGS[kind]['data'] == []):
            points = point_manager._generate_coordinates(size=self.SKS_SETTINGS[kind]['number'])
        # Case 2 - More points required than provided
        elif (self.SKS_SETTINGS[kind]['number'] > len(self.SKS_SETTINGS[kind]['data'])):
            n_points = self.SKS_SETTINGS[kind]['number'] - len(self.SKS_SETTINGS[kind]['data'])
            points = np.append(np.array(self.SKS_SETTINGS[kind]['data']), point_manager._generate_coordinates(n_points), axis=0)
        # Case 3 - Less points required than provided
        elif (self.SKS_SETTINGS[kind]['number'] < len(self.SKS_SETTINGS[kind]['data'])):
            points = self.rng['sks'].choice(self.SKS_SETTINGS[kind]['data'], self.SKS_SETTINGS[kind]['number'], replace=False)
        # Case 4 - Points required equals points declared
        else:
            points = self.SKS_SETTINGS[kind]['data']

        ### Populates the DataFrame
        x, y, z = zip(*points)
        data = {
            'x': x,
            'y': y,
            'z': z,
        }
        
        # Deletes the point manager
        del point_manager

        return pd.DataFrame(data=data)
    
    @wp._debug_level(3.8)
    @wp._parameters_validation('fractures', 'optional')
    @wp._memoize('fractures')
    @wp._logging()
    def _build_model_fractures(self) -> None:
        """Builds the fractures."""
        self._set_rng('fractures')
        fractures = self.SKS_SETTINGS['fractures']['data']
        test_a = isinstance(fractures, (str))
        test_b = (fractures == '')
        test_c = ('settings' in self.SKS_SETTINGS['fractures'])
        if (not (test_a and test_b)) or test_c:
            self.fractures = Fractures(rng=self.rng['fractures'])
            axis = self.SKS_SETTINGS['fractures']['axis']
            
            # Generates fractures families
            if 'settings' in self.SKS_SETTINGS['fractures']:
                frac_settings = self.SKS_SETTINGS['fractures']['settings']
                for frac_name, frac_settings in frac_settings.items():
                    self.fractures.generate_fracture_family(frac_name,
                                                            self.grid,
                                                            **frac_settings)
                self.fractures.generate_model()
            # Loads data
            else:
                self.fractures.set_data(fractures, self.grid, axis)
                costs = self.SKS_SETTINGS['fractures']['costs']
                self.fractures._set_costs(costs)
            
            self.fractures._compute_statistics(self.grid)
        return None
    
    @wp._debug_level(3.9)
    @wp._logging()
    def _build_conceptual_model(self) -> None:
        """Builds the conceptual model.
        
        Data range attributions :
         - 100 - 199 : Geology
         - 200 - 299 : Fractures
         - 300 - 399 : Faults
         - 0 - Out
        """
        # Initialization
        conceptual_model = np.zeros_like(self.grid.data_volume)
        conceptual_model_table = []
        
        # 100 - Geology
        geology_items = [(100 + i, 'Geology', id, cost) for (i, (id, cost)) in enumerate(self.geology.costs.items())]
        for (id, feature, id_, cost) in geology_items:
            conceptual_model = np.where(self.geology.data_volume == id_,
                                        id,
                                        conceptual_model)
        conceptual_model_table += geology_items
        
        # 200 - Fractures
        if self.fractures is not None:
            fractures_items = [(200 + i, 'Fractures', id, cost) for (i, (id, cost)) in enumerate(self.fractures.costs.items())]
            for (id, feature, id_, cost) in fractures_items:
                conceptual_model = np.where(self.fractures.data_volume == id_,
                                            id,
                                            conceptual_model)
            conceptual_model_table += fractures_items
            
        # 300 - Faults
        if self.faults is not None:
            faults_items = [(300 + i, 'Faults', id, cost) for (i, (id, cost)) in enumerate(self.faults.costs.items())]
            for (id, feature, id_, cost) in faults_items:
                conceptual_model = np.where(self.faults.data_volume == id_,
                                            id,
                                            conceptual_model)
            conceptual_model_table += faults_items
        
        # 0 - Out
        out_item = [(0, 'Out', np.nan, self.SKS_SETTINGS['sks']['costs']['out'])]
        conceptual_model = np.where(self.domain.data_volume == 0,
                                    0,
                                    conceptual_model)
        conceptual_model_table += out_item
        conceptual_model_table = pd.DataFrame(conceptual_model_table, columns=['id', 'feature', 'id-feature', 'cost']).set_index('id').sort_values('id')

        self.conceptual_model = conceptual_model
        self.conceptual_model_table = conceptual_model_table
        return None
    
    ##### 4 - FMM #########################
    @wp._debug_level(4)
    @wp._logging('fmm', 'construction')
    def _construct_fmm_variables(self) -> None:
        """Constructs the fast-marching method variables."""
        # Distributes inlets and outlets among karstic generations
        self._initialize_fmm_iterations()
        # Distributes inlets and outlets among computing iterations
        self._construct_fmm_iterations()
        # Initializes useful variables for the farst-marching method
        self._initialize_fmm_variables()   
        return None
        
    def _initialize_fmm_iterations(self):
        # Defining some variables
        outlets_nbr = len(self.outlets)
        inlets_nbr = len(self.inlets)
        self.outlets_importance = self.SKS_SETTINGS['outlets']['importance']
        self.inlets_importance = self.SKS_SETTINGS['inlets']['importance']
        inlets_per_outlet = self.SKS_SETTINGS['inlets']['per_outlet']

        # Calculating inlets and outlets repartitions
        self.nbr_iteration = len(self.outlets_importance) * len(self.inlets_importance)        # total number of iterations that will occur
        outlets_repartition = self._repartition_points(outlets_nbr, self.outlets_importance)    # correct for outlets_importance not summing to correct number of actual outlets
        inlets_repartition = self._repartition_points(inlets_nbr, inlets_per_outlet)          # correct for inlets_per_outlet not summing to correct number of actual inlets

        # Distributing inlets and outlets iterations
        outlets_distribution = pd.Series([k for (k, n) in enumerate(outlets_repartition) for j in range(n)], name='outlet_iteration')
        inlets_distribution = pd.Series([k for (k, n) in enumerate(inlets_repartition)  for j in range(n)], name='outlet_key')
        self._outlets = pd.concat([self.outlets, outlets_distribution], axis=1)  # store as a semi-private variable for internal use only
        self._inlets = pd.concat([self.inlets, inlets_distribution], axis=1)  # store as a semi-private variable for internal use only

        # Distributing iterations for each inlet
        for (outlet_key, row) in self._outlets.iterrows():
            inlets_test = self._inlets['outlet_key']==outlet_key
            inlets_current = self._inlets[inlets_test]
            inlets_nbr = len(inlets_current)
            repartition = self._repartition_points(inlets_nbr, self.inlets_importance)
            distribution = pd.Series([k for (k, n) in enumerate(repartition) for j in range(n)], name='inlet_iteration', index=inlets_current.index, dtype='object')
            self._inlets.loc[inlets_test, 'inlet_iteration'] = distribution

        return None
    
    def _construct_fmm_iterations(self):
        # Set up iteration structure:
        iteration = 0
        inlets = []
        outlets = []

        # Loops over outlet iterations
        for outlet_iteration in range(len(self.outlets_importance)):

            # Loops over inlet iterations
            for inlet_iteration in range(len(self.inlets_importance)):
                
                # Gets the outlets assigned to the current outlet iteration
                outlets_current = self._outlets[self._outlets['outlet_iteration'] == outlet_iteration]

                # Gets the inlets assigned to the current inlet iteration
                inlets_current = self._inlets[self._inlets['outlet_key'] == outlet_iteration]
                inlets_current_iteration = inlets_current[inlets_current['inlet_iteration'] == inlet_iteration]
                
                # Appends the list
                outlets.append(outlets_current.index.to_list())
                inlets.append(inlets_current_iteration.index.to_list())

                # Increments total iteration number by 1
                iteration = iteration + 1
        
        # Creates iterations dataframe
        self.iterations = pd.DataFrame({
            'outlets': outlets,
            'inlets': inlets,
        })
        self.iterations.index.name = 'iteration'
        
        return None
    
    def _repartition_points(self, nbr_points, importance):
        """
        Corrects for integers in importance factors list not summing correctly
        to total number of points.
        # TODO ?
        # """
        total_importance = float(sum(importance))                                               # total number of points as assigned (this may not be the actual total)
        proportion = [float(i) / total_importance for i in importance]                # percent of points to use this iteration
        repartition = [round(proportion[i] * nbr_points) for i in range(len(proportion))]    # number of points to use this iteration (need to round bc percentage will not result in whole number)
        repartition[-1] = nbr_points - sum(repartition[0:-1])                                    # leftover points are assignd to last iteration
        return repartition

    def _initialize_fmm_variables(self):
        ### Raster maps
        self.maps = {
            'outlets': np.full(self.grid.shape, np.nan),  # map of null values where each cell with an outlet will have the index of that outlet
            'nodes': np.full(self.grid.shape, np.nan),  # map of null values where each cell that has a node will be updated with that node index
            'cost': [],   # cost of travel through each cell
            'alpha': [],  # cost of travel along gradient through each cell
            'beta': [],   # cost of travel perpendicular to gradient through each cell
            'time': [],   # travel time to outlet from each cell
            'karst': [],  # presence/absence of karst conduit in each cell
        }
        
        ### Defines outlets map according to outlets emplacements
        for (i, outlet) in self._outlets.iterrows():  # assign outlet indices. Compute the outlets map (array indicating location of outlets as their index and everywhere else as nan).
            X = self.grid.get_i(outlet.x)
            Y = self.grid.get_j(outlet.y)
            Z = self.grid.get_k(outlet.z)
            self.maps['outlets'][X][Y][Z] = i

        ### Vector maps
        self.vectors = {
            'nodes': {},      # empty dict to store nodes (key: nodeID, val: [x, y, z, type])
            'edges': {},      # empty dict to store edges (key: edgeID, val: [inNode, outNode])
            'n': 0,           # start node counter at zero
            'e': 0,           # start edge counter at zero
            'geodesics': [],  # empty list to store raw fast-marching path output
        }

        ### Set up fast-marching
        self.fmm = {
            'algorithm': self.SKS_SETTINGS['sks']['algorithm'],
            'riemannMetric': [],                                            # this changes at every iteration, but cannot be stored?
            'fastMarching': agd.Eikonal.dictIn({
                'model': self.SKS_SETTINGS['sks']['algorithm'], # set algorithm from settings file ('Isotropic2', 'Isotropic3', 'Riemann2', 'Riemann3')
                'order': 2,                                     # recommended setting: 2 (replace by variable)
                'exportValues': 1,                                     # export the travel time field
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

###############################################################################
### KARST NETWORK SIMULATION ###
################################

    def compute(self) -> None:
        """
        Computes the karst network according to the parameters. Must be called
        after ``build()`` method. This method will :
        1. Computes conduits for each generation and stores nodes & edges for
        network.
        2. Stores all the relevant data for this network in specific
        dictionaries.
        3. Exports the state of the project, this file could be read by the
        'analysis' and 'visualization' sub-packages.

        Examples
        --------
        >>> import pykasso as pk
        >>> simulation = pk.SKS()
        >>> simulation.build()
        >>> simulation.compute()
        """
        self._compute_karst_network()
        self._export_results()
        self._export_project_state()
        return None
    
    def _compute_karst_network(self):
        this.logger = logging.getLogger("fmm.modelisation")
        this.logger.info("Computing karst network")
        
        self.iteration = 0
        
        for self.iteration in range(self.nbr_iteration):
        
            # Compute travel time maps and conduit network
            if self.fmm['algorithm'] == 'Isotropic3':
                self._compute_cost_map()        # 2.1.1
                self._compute_time_map()        # 2.1.2

            elif self.fmm['algorithm'] == 'Riemann3':
                self._compute_cost_map()        # 2.1.1
                self._compute_alpha_map()       # 2.1.4
                self._compute_beta_map()        # 2.1.5
                self._compute_riemann_metric()  # 2.1.6
                self._compute_time_map()        # 2.1.7
            
            self._compute_karst_map()
            self._voxelize_karst_network()

            msg = "iteration : {}/{}".format(self.iteration + 1,
                                             self.nbr_iteration)
            this.logger.info(msg)
        return None

    ### 2.1.1 Iso- and anisotropic case
    @wp._debug_level(1, True)
    @wp._logging()
    def _compute_cost_map(self):
        """
        Computes the cost map (how difficult it is to traverse each cell).

        TODO
        """
        # During the first iteration, iniatializes the cost map according to
        # the conceptual model.
        if self.iteration == 0:
            zeros_array = np.zeros((self.grid.nx, self.grid.ny, self.grid.nz))
            self.maps['cost'].append(zeros_array)
            for (i, row) in self.conceptual_model_table.iterrows():
                logical_test = self.conceptual_model == row.name
                self.maps['cost'][0] = np.where(logical_test,
                                                row.cost,
                                                self.maps['cost'][0])
                
        # During the rest of the iterations
        else:
            self.maps['cost'].append(self.maps['cost'][self.iteration - 1])
            
            self.maps['cost'][self.iteration] = (
                np.where(self.maps['karst'][self.iteration - 1] > 0,
                         self.SKS_SETTINGS['sks']['costs']['conduits'],
                         self.maps['cost'][self.iteration])
            )
            
        return None
    
    ### 2.1.4 Anisotropic case
    @wp._debug_level(2, True)
    @wp._logging()
    def _compute_alpha_map(self) -> None:
        """
        Computes the alpha map: travel cost in the same direction as the
        gradient.
        """
        alpha_map_0 = self.maps['cost'][self.iteration].copy()
        alpha_map = self.maps['cost'][self.iteration].copy()
        
        ### Option A
        if self.SKS_SETTINGS['sks']['mode'] == 'A':
            # Vadose zone = Phreatic zone = Bedrock zone
            domain = self.grid.data_volume
            alpha_map = self._set_map_from_elevation(alpha_map, domain)
            # Bedrock
            if self.domain.bedrock is not None:
                alpha_map = self._set_map_from_value(
                    alpha_map,
                    alpha_map.max(),
                    self.domain.bedrock.data_volume)
        
        ### Option B
        if self.SKS_SETTINGS['sks']['mode'] == 'B':
            # Vadose zone = Bedrock zone
            domain = self.domain.phreatic['vadose_zone']
            alpha_map = self._set_map_from_elevation(alpha_map, domain)
            # Bedrock
            if self.domain.bedrock is not None:
                alpha_map = self._set_map_from_value(
                    alpha_map,
                    alpha_map.max(),
                    self.domain.bedrock.data_volume)
            # Phreatic zone
            alpha_map = self._set_map_from_value(
                alpha_map,
                alpha_map_0,
                self.domain.phreatic['phreatic_zone'])
            
        ### Option C
        if self.SKS_SETTINGS['sks']['mode'] == 'C':
            # Vadose zone = Bedrock zone
            factor = 100
            alpha_map = self._set_map_from_value(
                alpha_map,
                alpha_map * factor,
                self.domain.phreatic['vadose_zone'])
            # Bedrock
            if self.domain.bedrock is not None:
                alpha_map = self._set_map_from_value(
                    alpha_map,
                    alpha_map.max() * 2,
                    self.domain.bedrock.data_volume)
            # Phreatic zone
            if self.domain.water_level is not None:
                alpha_map = self._set_map_from_value(
                    alpha_map,
                    alpha_map_0,
                    self.domain.phreatic['phreatic_zone'])
        
        ### Option D
        if self.SKS_SETTINGS['sks']['mode'] == 'D':
            # Vadose zone
            domain = np.logical_and(
                np.logical_not(self.domain.bedrock_domains['bedrock_vadose']),
                self.domain.phreatic['vadose_zone']
            )
            factor_01 = 100
            alpha_map = self._set_map_from_value(
                alpha_map,
                alpha_map * factor_01,
                domain)
            # Bedrock zone
            factor_02 = 50
            alpha_map = self._set_map_from_value(
                alpha_map,
                alpha_map * factor_02,
                self.domain.bedrock_domains['bedrock_vadose'])
            # Bedrock
            if self.domain.bedrock is not None:
                alpha_map = self._set_map_from_value(
                    alpha_map,
                    alpha_map.max() * 2,
                    self.domain.bedrock.data_volume)
            # Phreatic zone
            alpha_map = self._set_map_from_value(
                alpha_map,
                alpha_map_0,
                self.domain.phreatic['phreatic_zone'])
        
        self.maps['alpha'].append(alpha_map)
        return None
    
    def _set_map_from_elevation(self, array_map, domain) -> np.ndarray:
        """"""
        if (self.grid.nz == 1) and (self.domain.bedrock is not None):
            bedrock = self.domain.bedrock.data_surface.reshape(self.grid.shape)
            array_map = array_map * bedrock
            return array_map
        else:
            test = (domain == 1)
            new_array_map = array_map * self.grid.Z
            array_map = np.where(test, new_array_map, array_map)
            return array_map
    
    def _set_map_from_value(self, array_map, value, domain) -> np.ndarray:
        """ """
        test = (domain == 1)
        array_map = np.where(test, value, array_map)
        return array_map

    ### 2.1.5 Anisotropic case
    @wp._debug_level(3, True)
    @wp._logging()
    def _compute_beta_map(self) -> None:
        """
        Computes the beta map: travel cost perpendicular to the gradient. If
        beta is higher than alpha, conduits will follow the steepest gradient.
        If beta is lower than alpha, conduits will follow contours.
        """
        ratio = self.SKS_SETTINGS['sks']['costs']['ratio']
        alpha_map_0 = self.maps['cost'][self.iteration].copy()
        alpha_map = self.maps['alpha'][self.iteration].copy()
        
        ### Option A & default situation
        # Vadose zone = Phreatic zone = Bedrock zone
        beta_map = alpha_map / ratio
            
        ### Options B / C / D
        if self.SKS_SETTINGS['sks']['mode'] in ['B', 'C', 'D']:
            # Phreatic zone
            beta_map = self._set_map_from_value(
                beta_map,
                alpha_map_0,
                self.domain.phreatic['phreatic_zone'])
        
        # beta_map[:, :, :] = 1
        self.maps['beta'].append(beta_map)
        return None

    ### 2.1.6 Anisotropic case
    @wp._debug_level(4, True)
    @wp._logging()
    def _compute_riemann_metric(self) -> None:
        """
        Compute the riemann metric: Define the Riemannian metric needed as
        input for the anisotropic fast marching.

        TODO : à terminer
        """
        ### Option A & default situation
        grad_x = np.zeros_like(self.grid.data_volume)
        grad_y = np.zeros_like(self.grid.data_volume)
        grad_z = np.ones_like(self.grid.data_volume)
        
        ### Options B / C / D
        if self.SKS_SETTINGS['sks']['mode'] in ['B', 'C', 'D']:
            if self.domain.bedrock is not None:
                grad_x, grad_y, grad_z = (
                    self._set_gradient_from_bedrock(grad_x,
                                                    grad_y,
                                                    grad_z)
                )
        
        # Sets the needle
        alpha = self.maps['alpha'][self.iteration]
        beta = self.maps['beta'][self.iteration]
        self.orientation = (grad_x, grad_y, grad_z)  # TODO - to remove ??
        self.fmm['riemannMetric'] = agd.Metrics.Riemann.needle(
            [grad_x, grad_y, grad_z],
            alpha,
            beta
        )
        return None
    
    def _set_gradient_from_bedrock(self, grad_x, grad_y, grad_z) -> None:
        """"""
        test_grid_nx = (self.grid.nx == 1)
        test_grid_ny = (self.grid.ny == 1)
        test_grid_nz = (self.grid.nz == 1)
        test_domain = (self.domain.bedrock_domains['bedrock_vadose'] == 1)
        
        if test_grid_nx or test_grid_ny or test_grid_nz:
            if test_grid_nx:
                bedrock = np.roll(self.domain.bedrock.data_volume, 1, axis=2)
                bedrock[:, :, 0] = 1
                gradient_y, gradient_z = np.gradient(bedrock, self.grid.dy,
                                                     self.grid.dz, axis=(1, 2))
                gradient_x = np.full_like(gradient_y, 0)
            elif test_grid_ny:
                bedrock = np.roll(self.domain.bedrock.data_volume, 1, axis=2)
                bedrock[:, :, 0] = 1
                gradient_x, gradient_z = np.gradient(bedrock, self.grid.dx,
                                                     self.grid.dz, axis=(0, 2))
                gradient_y = np.full_like(gradient_x, 0)
            elif test_grid_nz:
                bedrock = self.domain.bedrock.data_surface
                gradient_x, gradient_y = np.gradient(bedrock, self.grid.dx,
                                                     self.grid.dy, axis=(0, 1))
                gradient_x = gradient_x.reshape(self.grid.shape)
                gradient_y = gradient_y.reshape(self.grid.shape)
                gradient_z = np.full_like(gradient_x, 0)
                bedrock = self.domain.bedrock.data_volume
            
            grad_x = np.where(test_domain, gradient_x, grad_x)
            grad_y = np.where(test_domain, gradient_y, grad_y)
            grad_z = np.where(test_domain, gradient_z, grad_z)
            
        else:
            # On commence par calculer le gradient en faisant gaffe à l'intervertion x y
            bedrock = self.domain.bedrock.data_surface
            gradient_x, gradient_y = np.gradient(bedrock,
                                                 self.grid.dx,
                                                 self.grid.dy)
            
            # Calcule du vecteur avec ses trois composantes (vx, vy, vz)

            # On créé les matrices de composantes vx, vy, vz vides
            vx = np.zeros(gradient_x.shape)
            vy = np.zeros(gradient_x.shape)
            vz = np.zeros(gradient_x.shape)
            
            # positions ou le gradient en x est non nul
            idx_grad_x_not0 = (gradient_x != 0)
            
            # vx = +/- 1 dans la direction oposée au gradient
            vx[idx_grad_x_not0] = - np.sign(gradient_x[idx_grad_x_not0])
            
            # Calcul de vy (pour respecter la direction horizontale)
            vy[idx_grad_x_not0] = (vx[idx_grad_x_not0]
                                   * gradient_y[idx_grad_x_not0]
                                   / gradient_x[idx_grad_x_not0])

            # On traite le cas particulier pour lequel gx = 0 et gy != 0
            idx_gx_is0 = ((gradient_x == 0) & (gradient_y != 0))

            # Dans ce cas la on normalise vy
            vy[idx_gx_is0] = - np.sign(gradient_y[idx_gx_is0])

            # Calcul de vz partout
            vz = gradient_x * vx + gradient_y * vy

            ### Dernier cas particulier problématique: la surface horizontale
            # Chercher les occurences
            idx_gxgy_are0 = ((gradient_x == 0) & (gradient_y == 0))

            # Par convention
            vx[idx_gxgy_are0] = 1
            vy[idx_gxgy_are0] = 0
            vz[idx_gxgy_are0] = 0

            ### Finalement on normalise le vecteur
            norm = np.sqrt(vx**2 + vy**2 + vz**2)
            vx /= norm
            vy /= norm
            vz /= norm
            
            args = np.argwhere(test_domain)
            i, j, k = zip(*args)

            grad_x[i, j, k] = -vx[i, j]
            grad_y[i, j, k] = -vy[i, j]
            grad_z[i, j, k] = -vz[i, j]
        
        out = (grad_x, grad_y, grad_z)
        return out

    ### 2.1.2
    @wp._debug_level(5, True)
    @wp._logging()
    def _compute_time_map(self):
        """
        Compute the travel time map (how long it takes to get to the outlet
        from each cell), using the ani- or isotropic agd-hfm fast-marching
        algorithm, and store travel time map.
        
        Note: the AGD-HFM library uses different indexing, so x and y indices
        are reversed for inlets and outlets. TODO - not so sure???
        TODO
        """
        logical_test = self.iterations.index == self.iteration
        # Set the outlets for this iteration
        outlets_ids = self.iterations[logical_test].outlets.values[0]
        iteration_outlets = self.outlets[self.outlets.index.isin(outlets_ids)]
        seeds_x = iteration_outlets.x
        seeds_y = iteration_outlets.y
        seeds_z = iteration_outlets.z
        seeds = list(zip(seeds_x, seeds_y, seeds_z))
        self.fmm['fastMarching']['seeds'] = seeds

        # Select inlets for current iteration
        inlets_ids = self.iterations[logical_test].inlets.values[0]
        iteration_inlets = self.inlets[self.inlets.index.isin(inlets_ids)]
        tips_x = iteration_inlets.x
        tips_y = iteration_inlets.y
        tips_z = iteration_inlets.z
        tips = list(zip(tips_x, tips_y, tips_z))
        self.fmm['fastMarching']['tips'] = tips

        # Set the travel cost through each cell
        if self.fmm['algorithm'] == 'Isotropic3':
            self.fmm['fastMarching']['cost'] = (
                self.maps['cost'][self.iteration]
            )
        elif self.fmm['algorithm'] == 'Riemann3':
            self.fmm['fastMarching']['metric'] = self.fmm['riemannMetric']

        # Set verbosity of hfm run
        verbosity = self.SKS_SETTINGS['verbosity']['agd']
        self.fmm['fastMarching']['verbosity'] = verbosity

        # Run the fast marching algorithm and store the outputs
        self.fmm['fastMarchingOutput'] = self.fmm['fastMarching'].Run()

        # Store travel time maps
        self.maps['time'].append(self.fmm['fastMarchingOutput']['values'])

        # Store fastest travel paths
        self.vectors['geodesics'].append(
            self.fmm['fastMarchingOutput']['geodesics']
        )
        return None

    ### 2.1.3
    @wp._debug_level(6, True)
    @wp._logging()
    def _compute_karst_map(self):
        """
        Compute the karst map based on the paths from agd-hfm.
        Array of all zeros, with ones in cells containing a karst conduit.
        """
        ### Loop over conduit paths generated by fast marching:
        for path in self.fmm['fastMarchingOutput']['geodesics']: #loop over conduit paths in this iteration (there is one path from each inlet)
            merge = False                                        #reset indicator for whether this conduit has merged with an existing conduit
            for p in range(path.shape[1]):                       #loop over points making up this conduit path
                point = path[:,p]                                #get coordinates of current point
                [[ix, iy, iz], error] = self.fmm['fastMarching'].IndexFromPoint(point) #convert to coordinates to indices, /!\ returning iy first then ix TODO ???
                #Place nodes and links:
                if np.isnan(self.maps['nodes'][ix, iy, iz]):                                    #if there is no existing conduit node here
                    if ~np.isnan(self.maps['outlets'][ix, iy, iz]):                              #if there is an outlet here (cell value is not nan)
                        outlet = self._outlets.iloc[int(self.maps['outlets'][ix, iy, iz])]         #get the outlet coordinates using the ID in the outlets map
                        self.vectors['nodes'][self.vectors['n']] = [outlet.x, outlet.y, outlet.z, 'outfall']           #add a node at the outlet coordinates (with the node type for SWMM)
                        self.maps['nodes'][ix, iy, iz] = self.vectors['n']                                   #update node map with node index
                        if p > 0:                                                           #if this is not the first point (i.e. the inlet) in the current path
                            if merge == False:                                               #if this conduit has not merged with an existing conduit
                                self.vectors['edges'][self.vectors['e']] = [self.vectors['n']-1, self.vectors['n']]                       #add an edge connecting the previous node to the current node
                                self.vectors['e'] = self.vectors['e'] + 1                                             #increment edge counter up by one
                            else:                                                          #if this conduit HAS merged with an existing conduit
                                [[fromix, fromiy, fromiz], error] = self.fmm['fastMarching'].IndexFromPoint(path[:, p-1]) #get xyz indices of previous point in current conduit path
                                n_from = self.maps['nodes'][fromix, fromiy, fromiz]           #get node index of the node already in the cell where the previous point was
                                self.vectors['edges'][self.vectors['e']] = [n_from, self.vectors['n']]                         #add an edge connecting existing conduit node to current node
                                self.vectors['e'] = self.vectors['e'] + 1                                             #increment edge counter up by one
                        self.vectors['n'] = self.vectors['n'] + 1                                                   #increment node counter up by one
                    else:                                                                  #if there is NOT an outlet here
                        if p > 0:                                                           #if this is not the first point in the current path
                            #possible improvement: if the next point on the path is on an existing point, skip the current point.
                            self.vectors['nodes'][self.vectors['n']] = [point[0], point[1], point[2], 'junction']            #add a junction node here (with the node type for SWMM)
                            self.maps['nodes'][ix, iy, iz] = self.vectors['n']                               #update node map with node index
                            if merge == False:                                              #if this conduit has not merged with an existing conduit
                                self.vectors['edges'][self.vectors['e']] = [self.vectors['n']-1, self.vectors['n']]                      #add and edge connecting the previous node to the current node
                                self.vectors['e'] = self.vectors['e'] + 1                                            #increment edge counter up by one
                            else:                                                           #if this conduit HAS merged with an existing conduit
                                [[fromix, fromiy, fromiz], error] = self.fmm['fastMarching'].IndexFromPoint(path[:, p-1]) #get xy indices of previous point in current conduit path
                                n_from = self.maps['nodes'][fromix, fromiy, fromiz]                   #get node index of the node already in the cell where the previous point was
                                self.vectors['edges'][self.vectors['e']] = [n_from, self.vectors['n']]                        #add an edge connecting existing conduit node to current node
                                self.vectors['e'] = self.vectors['e'] + 1                                            #increment edge counter up by one
                                merge = False                                                #reset merge indicator to show that current conduit has left                                                              #if this is the first point in current path
                        else:                                                                #if this is the first point in the current path (counter <= 0, therefore it is an inlet)
                            self.vectors['nodes'][self.vectors['n']] = [point[0], point[1], point[2], 'inlet']               #add an inlet node here (with the node type for SWMM)
                            self.maps['nodes'][ix, iy, iz] = self.vectors['n']                               #update node map with node index
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
                        else:                                                                #if this is the first point in the current path (i.e. the inlet is on an exising conduit)
                            self.vectors['nodes'][self.vectors['n']] = [point[0], point[1], point[2], 'inlet']                #add a node here (with the node type for SWMM)- this will cause there to be two nodes in the same cell
                            self.maps['nodes'][ix, iy, iz] = self.vectors['n']                                #update node map with node index
                            self.vectors['n'] = self.vectors['n'] + 1                                                 #increment node counter by 1
                # self.maps['karst'][self.iteration][ix, iy, iz] = 1                               #update karst map to put a conduit in current cell
        return None
    
    def _voxelize_karst_network(self) -> None:
        
        # Gets karst map from previous iteration
        # except for the very first iteration
        if self.iteration == 0:
            karst_map_0 = np.zeros(self.grid.shape)
            self.maps['karst'].append(karst_map_0)
        else:
            karst_map_copy = self.maps['karst'][self.iteration - 1].copy()
            self.maps['karst'].append(karst_map_copy)

        # Cleans edges and nodes
        edges = list(self.vectors['edges'].values())
        nodes = {}
        for i_node in self.vectors['nodes']:
            nodes[i_node] = self.vectors['nodes'][i_node][:3]

        # Retrieves the points
        points = []
        for (i1, i2) in edges:
            point1 = nodes[i1]
            point2 = nodes[i2]
            x1, y1, z1 = point1
            x2, y2, z2 = point2
            x = x2 - x1
            y = y2 - y1
            z = z2 - z1
            x3 = x1 + x / 2
            y3 = y1 + y / 2
            z3 = z1 + z / 2
            points.extend([(x1, y1, z1)])
            points.extend([(x2, y2, z2)])
            points.extend([(x3, y3, z3)])

        # Gets the indices
        try:
            x, y, z = zip(*points)
            i, j, k = self.grid.get_indices(x, y, z)

            # Calculates the new karst array
            new_karst = np.zeros_like(self.grid.data_volume)
            new_karst[i, j, k] = 1
            self.maps['karst'][self.iteration] = new_karst
        except:
            pass
        
        return None
    
    def _export_results(self):
        """
        TODO
        """
        path = self.CORE_SETTINGS['simulation_directory'] + 'results'
        with open(path + '.pickle', 'wb') as handle:
            results = {
                'maps': self.maps.copy(),
                'vectors': self.vectors.copy(),
                'inlets': self.inlets.copy(),
                'outlets': self.outlets.copy(),
                'grid': self.grid,
                'domain': self.domain,
                'geology': self.geology,
                'faults': self.faults,
                'fractures': self.fractures,
                # TODO
                'orientation': self.orientation,
            }
            pickle.dump(results, handle, protocol=pickle.HIGHEST_PROTOCOL)
        return None
    
    def _export_project_state(self):
        """
        TODO
        """
        path = self.CORE_SETTINGS['outputs_directory'] + 'project_state'
        with open(path + '.yaml', 'w') as handle:
            STATE_PROJECT = {
                'grid': this.ACTIVE_PROJECT['settings']['grid'],
                'project_directory': this.ACTIVE_PROJECT['project_directory'],
                'n_simulation': this.ACTIVE_PROJECT['n_simulation'],
                'simulation_locations': this.ACTIVE_PROJECT['simulation_locations'],
            }
            yaml.dump(STATE_PROJECT, handle)
        return None
