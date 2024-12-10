"""
This module contains a class modeling the karstic network generator tool.
"""

### Internal dependencies
import os
import copy
import pickle
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
import pykasso.model._wrappers as wp
from .domain import Domain, Delimitation, Topography, Bedrock, WaterTable
from .geologic_features import Geology, Faults
from .fracturation import Fractures
from pykasso._utils.array import normalize_array
from pykasso.core._namespaces import DEFAULT_FMM_COSTS

### Typing
from pykasso._typing import Project


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
    
    def __init__(self,
                 project: Project,
                 ) -> None:
        """_summary_

        Returns
        -------
        _type_
            _description_
        """
        ### Initialization
        self.project = project
        self.grid = project.grid

    def generate(self,
                 model_parameters: dict = {},
                 export_settings: dict = {},
                 ) -> None:
        """_summary_

        Parameters
        ----------
        model_parameters : dict, optional
            _description_, by default {}
        export_settings : dict, optional
            _description_, by default {}

        """
        self.model_parameters = model_parameters
        self.export_settings = export_settings
        
        ### Initialize the model
        self._initialize()
        
        ### Build the geological model
        self._build()
        
        ### Compute the karst network
        self._compute()
        
        return None
    
    def _initialize(self) -> None:
        """
        Initialize and configure the basic settings.
        """
        # Initialize class attributes
        self.domain = None
        self.geology = None
        self.faults = None
        self.fractures = None
        self.inlets = None
        self.outlets = None
        self.conceptual_model = None
        self.conceptual_model_table = None
        self.rng = {}
        
        # Load model parameters
        self._load_model_parameters()
        
        # Create simulation directory
        self._create_sim_dir()
        
        # Check verbosity levels
        self._check_verbosity()
        
        # Update the log
        self._update_log()
    
        return None
    
    def _load_model_parameters(self) -> None:
        """
        Load the model parameters.
        """
        model_parameters = copy.deepcopy(self.model_parameters)
        model_parameters_f = self.project.core['filenames']['parameters']
        
        # If model parameters dictionary is empty then load default conf file
        if model_parameters == {}:
            project_dir = self.project.core['paths']['inputs_dir']
            model_parameters = project_dir + model_parameters_f
        if isinstance(model_parameters, str):
            with open(model_parameters, 'r') as f:
                self.model_parameters = yaml.safe_load(f)
        if isinstance(model_parameters, dict):
            self.model_parameters = model_parameters
        
        # Inject the inputs directory within the data paths
        inputs_dir = self.project.core['paths']['inputs_dir']
        for key, dict_parameters in self.model_parameters.items():
            if key == 'domain':
                for key_, value in dict_parameters.items():
                    if not isinstance(value, str):
                        continue
                    new_path = inputs_dir + value
                    dict_parameters[key_] = new_path
                self.model_parameters['domain'].update(dict_parameters)
            else:
                for key_, value in dict_parameters.items():
                    if key_ == 'data':
                        if not isinstance(value, str):
                            continue
                        new_path = inputs_dir + value
                        dict_parameters[key_] = new_path
        return None
    
    def _create_sim_dir(self) -> None:
        """
        Create simulation directory
        """
        self.project._increment_n_simulations()
        outputs_dir = self.project.core['paths']['outputs_dir']
        sim_dir = 'simulation_{}/'.format(self.project.n_simulations)
        sim_path = outputs_dir + sim_dir
        os.makedirs(sim_path, exist_ok=True)
        self.project.simulations.append(sim_path)
        self.simulation_path = sim_path
    
    def _check_verbosity(self) -> None:
        """
        Check if verbosity levels dictionary has been declared
        """
        if 'verbosity' not in self.model_parameters:
            self.model_parameters['verbosity'] = {
                'logging': 0,
                'agd': 0,
            }
        else:
            if 'logging' not in self.model_parameters['verbosity']:
                self.model_parameters['logging'] = 0
            if 'agd' not in self.model_parameters['verbosity']:
                self.model_parameters['agd'] = 0
        return None
        
    def _update_log(self) -> None:
        """
        Print information about the current simulation in the log.
        """
        # Set logging level
        level = self.model_parameters['verbosity']['logging']
        if level > 4:
            level = 4
        if level < 0:
            level = 0
        logging_level = self.project._logging_levels[level]
        self.logger = logging.getLogger('♠').setLevel(logging_level)
        
        # Print current simulation number
        n_sim = self.project.n_simulations
        l_sim = len(str(n_sim))
        date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.logger = logging.getLogger("►")
        self.logger.info('***********************{}****'.format('*' * l_sim))
        self.logger.info('*** pyKasso simulation {} ***'.format(n_sim))
        self.logger.info('***********************{}****'.format('*' * l_sim))
        self.logger.info('Path          : ' + self.simulation_path)
        self.logger.info('Date          : ' + date)
        self.logger.info("Logging level : " + str(logging_level))
        self.logger.info('---')
        return None
    
##########################
### BUILDING THE MODEL ###
##########################

    def _build(self) -> None:
        """
        Build the geological model.
        
        Model building procedure:
         1. Define the main model parameters: ``self.rng``;
         2. Build the geological features: ``self.geology``, ``self.faults``,
         ``self.fractures``;
         3. Build the model extension: ``self.domain``;
         4. Build the point features: ``self.outlets``, ``self.inlets``;
         5. Construct the conceptual model: ``self.conceptual_model``,
         ``self.conceptual_model_table``;
         6. Prepare the fast-marching method variables.
         
        If a geological feature has not been defined in the settings, calling
        its attribute will return ``None``.
        """
        self._build_model_parameters()     # Set the main parameters
        self._build_geological_features()  # Construct the geological features
        self._build_domain()               # Construct the domain
        self._build_point_features()       # Construct the points features
        self._build_conceptual_model()     # Construct the conceptual model
        self._construct_fmm_variables()    # Construct fmm features
        return None
    
    ##### MODEL PARAMETERS #########################
    @wp._parameters_validation('sks', 'optional')
    def _build_model_parameters(self) -> None:
        """
        Build main model parameters:
         - Random Number Generator
        """
        # Define the main seed
        if self.model_parameters['sks']['seed'] == 0:
            main_seed = np.random.default_rng().integers(low=0, high=10**9)
            self.model_parameters['sks']['seed'] = main_seed
        else:
            main_seed = self.model_parameters['sks']['seed']
        self.rng['sks'] = np.random.default_rng(main_seed)
        return None
    
    def _set_rng(self, attribute: str) -> None:
        """
        Assign the seed for the specified attribute. If the seed has
        been set by the user, then assign that value, otherwise either
        generate a random seed triggered by the global 'sks' seed if the
        seed is None, or generate a completely independent random seed
        if the seed is 0.
        """
        
        # If the actual seed is None, derive a seed from 'sks' seed.
        if self.model_parameters[attribute]['seed'] is None:
            seed = self.rng['sks'].integers(low=0, high=10**9)

        # If the seed is 0, generate an independent random seed
        elif self.model_parameters[attribute]['seed'] == 0:
            seed = np.random.default_rng().integers(low=0, high=10**9)

        # Otherwise, get the seed from the dictionary if it exists
        else:
            seed = self.model_parameters[attribute]['seed']
        
        self.rng[attribute] = np.random.default_rng(seed)
        return None
    
    ##### GEOLOGICAL FEATURES #########################
    def _build_geological_features(self) -> None:
        """
        Build the geological features:
         - Geology model
         - Faults model
         - Fracturation model
        """
        self._build_model_geology()
        self._build_model_faults()
        self._build_model_fractures()
        return None
    
    @wp._parameters_validation('geology', 'optional')
    @wp._memoize('geology')
    @wp._logging()
    def _build_model_geology(self) -> None:
        """
        Build the geology.
        """
        # Retrieve geology parameters
        geology_settings = self.model_parameters['geology']

        # Test value of 'data'
        if not isinstance(geology_settings['data'], np.ndarray):
            if geology_settings['data'] == '':
                geology_settings['data'] = None
        
        # Create the geology
        self.geology = Geology(grid=self.grid, **geology_settings)
        return None
    
    @wp._parameters_validation('faults', 'optional')
    @wp._memoize('faults')
    @wp._logging()
    def _build_model_faults(self) -> None:
        """
        Build the faults.
        """
        # Retrieve geology parameters
        faults_settings = self.model_parameters['faults']
        
        # Test value of 'data'
        if not isinstance(faults_settings['data'], np.ndarray):
            if faults_settings['data'] == '':
                faults_settings['data'] = None
                
        # Create the faults
        self.faults = Faults(grid=self.grid, **faults_settings)
        return None
    
    @wp._parameters_validation('fractures', 'optional')
    @wp._memoize('fractures')
    @wp._logging()
    def _build_model_fractures(self) -> None:
        """
        Build the fractures.
        """
        # Retrieve fractures parameters
        fractures_settings = self.model_parameters['fractures']
        
        # Test value of 'data'
        if not isinstance(fractures_settings['data'], np.ndarray):
            if fractures_settings['data'] == '':
                fractures_settings['data'] = None
        
        # Set the random seed for the fracturation
        self._set_rng('fractures')
        
        # Create the fractures
        self.fractures = Fractures(grid=self.grid,
                                   rng=self.rng['fractures'])
            
        # Generate fractures families
        if 'generate' in fractures_settings:
            frac_settings = fractures_settings['generate']
            for frac_name, frac_settings in frac_settings.items():
                default_cost = DEFAULT_FMM_COSTS['fractures']
                cost = frac_settings.pop('cost', default_cost)
                self.fractures.generate_fracture_family(frac_name,
                                                        frac_settings,
                                                        cost)
            self.fractures.compute_model()
            self.fractures.compute_statistics()
        return None
    
    ##### DOMAIN #########################
    @wp._parameters_validation('domain', 'optional')
    @wp._memoize('domain')
    @wp._logging()
    def _build_domain(self) -> None:
        """
        Build the domain.
        """
        # Define each limit and build the final domain
        delimitation = self._build_domain_delimitation()
        topography = self._build_domain_topography()
        bedrock = self._build_domain_bedrock()
        water_table = self._build_domain_water_table()
        self.domain = Domain(self.grid,
                             delimitation=delimitation,
                             topography=topography,
                             bedrock=bedrock,
                             water_table=water_table,
                             geology=self.geology.get_data_model())
        return None
    
    @wp._logging()
    def _build_domain_delimitation(self) -> None:
        """
        Build the delimitation.
        """
        delimitation = self.model_parameters['domain']['delimitation']
        if delimitation is not None:
            if isinstance(delimitation, (str)):
                delimitation = np.genfromtxt(delimitation).tolist()
                self.model_parameters['domain']['delimitation'] = delimitation
            instance = Delimitation(vertices=delimitation, grid=self.grid)
            return instance
        else:
            return None

    @wp._logging()
    def _build_domain_topography(self) -> None:
        """
        Build the topography.
        """
        topography = self.model_parameters['domain']['topography']
        if topography is not None:
            instance = Topography(grid=self.grid, data=topography)
            instance.data_volume = instance._surface_to_volume('<=', self.grid)
            return instance
        else:
            return None
        
    @wp._logging()
    def _build_domain_bedrock(self) -> None:
        """
        Build the bedrock elevation.
        """
        bedrock = self.model_parameters['domain']['bedrock']
        if bedrock is not None:
            instance = Bedrock(grid=self.grid, data=bedrock)
            instance.data_volume = instance._surface_to_volume('<=', self.grid)
            return instance
        else:
            return None
        
    @wp._logging()
    def _build_domain_water_table(self) -> None:
        """
        Build the water table.
        """
        water_table = self.model_parameters['domain']['water_table']
        if water_table is not None:
            instance = WaterTable(grid=self.grid, data=water_table)
            instance.data_volume = instance._surface_to_volume('<=', self.grid)
            return instance
        else:
            return None

    ##### POINT FEATURES #########################
    def _build_point_features(self) -> None:
        """
        Build the point features:
        - Outlet locations
        - Inlet locations
        """
        self._points = {
            'outlets': {},
            'inlets': {},
        }
        self._build_model_outlets()
        self._build_model_inlets()
        # del self._points
        return None
    
    @wp._parameters_validation('outlets', 'required')
    @wp._memoize('outlets')
    @wp._logging()
    def _build_model_outlets(self) -> None:
        """
        Build the outlets.
        """
        self._set_rng('outlets')
        self.outlets = self._construct_feature_points('outlets')
        if self.model_parameters['outlets']['shuffle']:
            pts = self.outlets.sample(frac=1, random_state=self.rng['sks'])
            self.outlets = pts.reset_index(drop=True)
        return None

    @wp._parameters_validation('inlets', 'required')
    @wp._memoize('inlets')
    @wp._logging()
    def _build_model_inlets(self) -> None:
        """
        Build the inlets.
        """
        self._set_rng('inlets')
        self.inlets = self._construct_feature_points('inlets')
        if self.model_parameters['inlets']['shuffle']:
            pts = self.inlets.sample(frac=1, random_state=self.rng['sks'])
            self.inlets = pts.reset_index(drop=True)
        return None
    
    def _construct_feature_points(self, kind: str) -> pd.DataFrame:
        """
        Construct a set of points.

        Four cases possible:
         1. No points declared
         2. More points required than provided : Generates additional random
         points
         3. Less points required than provided : Pick random points among
         provided ones
         4. Points required equals points declared
         
         Parameters
         ----------
         kind : str
            Kind of points : "outlets" or "inlets".
         
         Returns
         -------
         pd.DataFrame
        """
        logger = logging.getLogger("construction.{}".format(kind))
        
        ### Prepare data for point generation
        self._define_point_generation_vars(kind)
        
        ### Inspect validity of declared points
        points = self.model_parameters[kind]['data']
        subdomain = self.model_parameters[kind]['subdomain']
        validated_points = []
        if len(points) != 0:
            for point in points:
                if len(point) == 2:
                    if self._is_point_valid(point, kind):
                        point = self._2D_to_3D_point(kind, point)
                        validated_points.append(point)
                    else:
                        msg = ("Point with coordinates {} is invalid in"
                               " subdomain '{}'.").format(point, subdomain)
                        logger.error(msg)
                        raise ValueError(msg)
                elif len(point) == 3:
                    if self._is_point_valid(point, kind):
                        validated_points.append(point)
                    else:
                        msg = ("Point with coordinates {} is invalid in"
                               " subdomain '{}'.").format(point, subdomain)
                        logger.error(msg)
                        raise ValueError(msg)
                else:
                    msg = "Point with coordinates {} is invalid.".format(point)
                    logger.error(msg)
                    raise ValueError(msg)
            self.model_parameters[kind]['data'] = validated_points

        ### Get new points according to the right case
        n = self.model_parameters[kind]['number']
        data = self.model_parameters[kind]['data']
        
        # Case 1 - No points declared
        if (data == '') or (data == []):
            points = self._generate_points(kind, size=n)
        # Case 2 - More points required than provided
        elif (n > len(data)):
            n_points = n - len(data)
            points = np.append(
                np.array(data),
                self._generate_points(kind, n_points),
                axis=0
            )
        # Case 3 - Less points required than provided
        elif (n < len(data)):
            # points = self.rng['sks'].choice(data, n, replace=False) # ?
            points = data[:n]
        # Case 4 - Points required equals points declared
        else:
            points = data

        ### Populate the DataFrame
        x, y, z = zip(*points)
        data = {
            'x': x,
            'y': y,
            'z': z,
        }

        return pd.DataFrame(data=data)
    
    def _define_point_generation_vars(self, kind: str) -> None:
        """
        TODO
        """
        # Retrieve subdomain model
        subdomain = self.model_parameters[kind]['subdomain']
        self._points[kind]['subdomain'] = self.domain.get_subdomain(subdomain)
        
        # Retrieve valid geological model
        geologic_ids = self.model_parameters[kind]['geology']
        if geologic_ids is not None:
            valid_geologic_ids = self._controls_geologic_ids(geologic_ids)
        else:
            df = self.geology.overview()
            valid_geologic_ids = df[df['model']==True].index.to_list()
        self._points[kind]['geology'] = self.geology.get_data_units(valid_geologic_ids)
        
        # Cross subdomain with geological domain
        self._points[kind]['data'] = np.logical_and(
            self._points[kind]['subdomain'] > 0,
            self._points[kind]['geology'] > 0,
        ).astype(int)
        
        # Define valid cells
        indices = np.indices(self.grid.shape)
        filter = self._points[kind]['data'].astype(bool)
        valid_indices = indices[:, filter]
        valid_cells = pd.DataFrame(valid_indices.T, columns=['i', 'j', 'k'])
        self._points[kind]['valid_cells'] = valid_cells
        return None
    
    def _controls_geologic_ids(self, geologic_ids: list) -> list:
        """
        Control if the geologic IDs specified in the ``geologic_ids`` list
        match with the geologic IDs from the geologic model.
        
        Returns
        -------
        list
            List of matching geologic IDs
        """
        values = self.geology.stats.index.to_list()
        validated_geology_ids = []
        
        for geologic_id in geologic_ids:
            if geologic_id in values:
                validated_geology_ids.append(geologic_id)
            else:
                msg = ("Declared geologic ID #{} is not present in "
                       "geology model.".format(geologic_id))
                logging.warning(msg)
            
        if len(validated_geology_ids) == 0:
            msg = ("None of the geologic IDs declared are present in the "
                   "geologic model, geologic constraints are ignored.")
            logging.warning(msg)

        return validated_geology_ids
    
    def _is_point_valid(self, point: tuple, kind: str) -> bool:
        """
        Check if a 2D or 3D point is in domain.
        
        Returns
        -------
        bool
        """
        if self.grid.is_inbox(point):
            if len(point) == 2:
                i, j = self.grid.get_indices(point)
                # out = self.domain.data_surfaces['z'][i, j] > 0
                out = self._points[kind]['data'].max(axis=2)[i, j] > 0
            elif len(point) == 3:
                i, j, k = self.grid.get_indices(point)
                # out = self.domain.data_volume[i, j, k] > 0
                out = self._points[kind]['data'][i, j, k] > 0
            return out
        else:
            return False
        
    def _2D_to_3D_point(self, kind: str, point: tuple) -> tuple:
        """
        TODO
        """
        i, j = self.grid.get_indices(point)
        i = i[0]
        j = j[0]
        valid_cells = self._points[kind]['valid_cells']
        column_valid_cells = valid_cells[(valid_cells['i'] == i)
                                         & (valid_cells['j'] == j)]
        i_, j_, k_ = self.rng[kind].choice(column_valid_cells)
        z = (self.grid.zmin + (k_ + self.rng[kind].random()) * self.grid.dz)
        x, y = point
        return (x, y, z)
    
    def _generate_points(self, kind: str, size: int = 1) -> np.ndarray:
        """
        TODO
        """
        RNG = self.rng[kind]
        indices = RNG.choice(self._points[kind]['valid_cells'], size=size)
        i, j, k = zip(*indices)
        i, j, k = np.array(i), np.array(j), np.array(k)
        x = (self.grid.xmin + (i + RNG.random()) * self.grid.dx)
        y = (self.grid.ymin + (j + RNG.random()) * self.grid.dy)
        z = (self.grid.zmin + (k + RNG.random()) * self.grid.dz)
        return np.dstack((x, y, z))[0]
    
    @wp._logging()
    def _build_conceptual_model(self) -> None:
        """
        Build the conceptual model.
        
        Data range attributions:
         - 100 - 199: Geology
         - 200 - 299: Fractures
         - 300 - 399: Faults
         - 0 - Out of domain
        """
        ### Initialization
        conceptual_model = np.zeros_like(self.grid.data_volume)
        items = []
        
        ### 100-199: Geology
        df_geology = self.geology.overview()
        df_geology = df_geology[df_geology['model'].isin([True, 1])]
        df_geology = df_geology[['names', 'costs']]
        dt_geology = df_geology.to_dict(orient='index')
        
        geology_items = [(100 + i, 'Geology', id,
                          dict_['names'], dict_['costs'])
                         for (i, (id, dict_))
                         in enumerate(dt_geology.items())]
        
        for (model_id, feature, data_id, name, cost) in geology_items:
            conceptual_model = np.where(self.geology.data_volume == data_id,
                                        model_id,
                                        conceptual_model)
        
        items += geology_items
        
        ### 200-299: Fractures
        df_fractures = self.fractures.overview()
        df_fractures = df_fractures[df_fractures['model'].isin([True, 1])]
        df_fractures = df_fractures[['names', 'costs']]
        dt_fractures = df_fractures.to_dict(orient='index')
        
        fractures_items = [(200 + i, 'Fractures', id,
                            dict_['names'], dict_['costs'])
                           for (i, (id, dict_))
                           in enumerate(dt_fractures.items())]
        
        for (model_id, feature, data_id, name, cost) in fractures_items:
            conceptual_model = np.where(self.fractures.data_volume == data_id,
                                        model_id,
                                        conceptual_model)
        
        items += fractures_items
            
        ### 300-399: Faults
        df_faults = self.faults.overview()
        df_faults = df_faults[df_faults['model'].isin([True, 1])]
        df_faults = df_faults[['names', 'costs']]
        dt_faults = df_faults.to_dict(orient='index')
        
        faults_items = [(300 + i, 'Faults', id,
                         dict_['names'], dict_['costs'])
                        for (i, (id, dict_))
                        in enumerate(dt_faults.items())]
        
        for (model_id, feature, data_id, name, cost) in faults_items:
            conceptual_model = np.where(self.faults.data_volume == data_id,
                                        model_id,
                                        conceptual_model)
        
        items += faults_items
        
        ### 0 - Out of domain
        cost_out = self.model_parameters['sks']['costs']['out']
        out_item = [(0, 'Out', np.nan, np.nan, cost_out)]
        conceptual_model = np.where(self.domain.data_volume == 0,
                                    0,
                                    conceptual_model)
        items += out_item
        
        ### Create data
        columns = ['model_id', 'feature', 'data_id', 'name', 'cost']
        table = pd.DataFrame(items, columns=columns)
        table = table.set_index('model_id').sort_values('model_id')

        self.conceptual_model = conceptual_model
        self.conceptual_model_table = table
        return None
    
    ##### 4 - FMM #########################
    @wp._logging('fmm', 'construction')
    def _construct_fmm_variables(self) -> None:
        """
        Construct the fast-marching method variables.
        """
        # Distribute inlets and outlets among karstic generations
        self._initialize_fmm_iterations()
        # Distribute inlets and outlets among computing iterations
        self._construct_fmm_iterations()
        # Initialize useful variables for the farst-marching method
        self._initialize_fmm_variables()
        return None
        
    def _initialize_fmm_iterations(self):
        # Defining some variables
        outlets_nbr = len(self.outlets)
        inlets_nbr = len(self.inlets)
        self.outlets_importance = self.model_parameters['outlets']['importance']
        self.inlets_importance = self.model_parameters['inlets']['importance']
        
        # Calculate the total number of iterations that will occur
        outlets_importance_len = len(self.outlets_importance)
        inlets_importance_len = len(self.inlets_importance)
        self.nbr_iteration = outlets_importance_len * inlets_importance_len
        
        # Compute outlets partition
        outlets_partition = self._partition_points(outlets_nbr,
                                                   self.outlets_importance)
        
        # Compute inlets partition
        inlets_partition = self._partition_points(inlets_nbr,
                                                  outlets_partition)
        
        # TODO - per_outlet optional
        # per_outlets = [1] * outlets_nbr
        # self.model_parameters['inlets'].setdefault('per_outlet', per_outlets)
        # inlets_per_outlet = self.model_parameters['inlets']['per_outlet']
        # inlets_repartition = self._partition_points(inlets_nbr,
        #                                             inlets_per_outlet)
        
        # Distribute outlets iterations
        outlets_distribution = pd.Series([k for (k, n) in enumerate(outlets_partition) for j in range(n)], name='outlet_iteration')
        self._outlets = pd.concat([self.outlets, outlets_distribution], axis=1)
        
        # Distribute inlets iterations
        inlets_distribution = pd.Series([k for (k, n) in enumerate(inlets_partition) for j in range(n)], name='outlet_key')
        self._inlets = pd.concat([self.inlets, inlets_distribution], axis=1)

        # Distribute iterations in inlets
        for (outlet_key, row) in self._outlets.iterrows():
            inlets_test = self._inlets['outlet_key'] == outlet_key
            inlets_current = self._inlets[inlets_test]
            inlets_nbr = len(inlets_current)
            partition = self._partition_points(inlets_nbr, self.inlets_importance)
            distribution = pd.Series([k for (k, n) in enumerate(partition) for j in range(n)], name='inlet_iteration', index=inlets_current.index, dtype='object')
            self._inlets.loc[inlets_test, 'inlet_iteration'] = distribution

        return None
    
    def _construct_fmm_iterations(self):
        # Set up iteration structure:
        iteration = 0
        inlets = []
        outlets = []

        # Loop over outlet iterations
        for outlet_iteration in range(len(self.outlets_importance)):

            # Loop over inlet iterations
            for inlet_iteration in range(len(self.inlets_importance)):
                
                # Retrieve the outlets assigned to the current outlet iteration
                test = self._outlets['outlet_iteration'] == outlet_iteration
                outlets_current = self._outlets[test]

                # test = self._inlets['inlet_iteration'] == inlet_iteration
                # inlets_current = self._inlets[test]
                # TODO - per_outlet
                # Retrieve the inlets assigned to the current inlet iteration
                test = self._inlets['outlet_key'] == outlet_iteration
                inlets_current = self._inlets[test]
                test = inlets_current['inlet_iteration'] == inlet_iteration
                inlets_current_iteration = inlets_current[test]
                
                # Append the list
                outlets.append(outlets_current.index.to_list())
                # inlets.append(inlets_current.index.to_list())
                inlets.append(inlets_current_iteration.index.to_list())

                # Increment iteration counter
                iteration = iteration + 1
        
        # Create iterations dataframe
        self.iterations = pd.DataFrame({
            'outlets': outlets,
            'inlets': inlets,
        })
        self.iterations.index.name = 'iteration'
        
        # Shuffle
        # self.iterations = self.iterations.sample(frac=1, ignore_index=True)
        
        #, random_state=self.rng['sks'])
        
        return None
    
    def _partition_points(self, nbr_points, importance):
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
            'cost': [],         # cost of travel through each cell
            'alpha': [],        # cost of travel along gradient through each cell
            'beta': [],         # cost of travel perpendicular to gradient through each cell
            'gradient': [],     # bebdrock gradient
            'time': [],         # travel time to outlet from each cell
            'karst': [],        # presence/absence of karst conduit in each cell
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
            'algorithm': self.model_parameters['sks']['algorithm'],
            'riemannMetric': [],                                            # this changes at every iteration, but cannot be stored?
            'fastMarching': agd.Eikonal.dictIn({
                'model': self.model_parameters['sks']['algorithm'], # set algorithm from settings file ('Isotropic2', 'Isotropic3', 'Riemann2', 'Riemann3')
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

    def _compute(self) -> None:
        """
        Compute the karst conduit network according to the parameters. It will:
        1. Compute conduits for each generation and store the nodes & edges.
        2. Store all the relevant data for this network in specific
        dictionaries.
        3. Export the state of the project, this file could be read by the
        'analysis' and 'visualization' sub-packages.
        """
        self._compute_karst_network()
        self._export_results()
        self.project._export_project_file()
        return None
    
    def _compute_karst_network(self) -> None:
        """
        TODO
        """
        # Logging operations
        self.logger = logging.getLogger("modelisation.karst")
        self.logger.info("Computing karst conduit network...")
        
        # Initialize the iteration counter
        self.iteration = 0
        
        # Start the conduit generation loop
        for self.iteration in range(self.nbr_iteration):
        
            # Conduits generation with isotropic fast marching
            if self.fmm['algorithm'] == 'Isotropic3':
                self._compute_cost_map()
                self._compute_time_map()

            # Conduits generation with anisotropic fast marching
            elif self.fmm['algorithm'] == 'Riemann3':
                self._compute_cost_map()
                self._compute_alpha_map()
                self._compute_beta_map()
                self._compute_gradient()
                self._compute_riemann_metric()
                self._compute_time_map()
            
            # Transform results from fast marching into karst conduit network
            self._compute_karst_map()
            
            # Voxelize the karst conduit network
            self._voxelize_karst_network()

            # Log the current iteration
            msg = "iteration : {}/{}".format(self.iteration + 1,
                                             self.nbr_iteration)
            self.logger.info(msg)
        return None

    @wp._logging()
    def _compute_cost_map(self) -> None:
        """
        Compute the cost map (how difficult it is to traverse each cell).
        
        TODO
        """
        # If first iteration, iniatialize the cost map according to the
        # conceptual model.
        if self.iteration == 0:
            zeros_array = np.zeros(self.grid.shape)
            self.maps['cost'].append(zeros_array)
            for (i, row) in self.conceptual_model_table.iterrows():
                logical_test = self.conceptual_model == row.name
                self.maps['cost'][0] = np.where(logical_test,
                                                row.cost,
                                                self.maps['cost'][0])
                
        # During the rest of the iterations
        else:
            # Append a copy of the cost map from the previous iteration
            self.maps['cost'].append(self.maps['cost'][self.iteration - 1])
            
            # Update the array with the cost of the karstic conduits computed
            # during the previous iteration
            self.maps['cost'][self.iteration] = (
                np.where(self.maps['karst'][self.iteration - 1] > 0,
                         self.model_parameters['sks']['costs']['conduits'],
                         self.maps['cost'][self.iteration])
            )
        return None
    
    @wp._logging()
    def _compute_alpha_map(self) -> None:
        """
        Compute the alpha map: travel cost in the same direction as the
        gradient.
        """
        # Retrieve the cost map from the current iteration
        cost_map = self.maps['cost'][self.iteration].copy()
        
        # Select the appropriate method to set the alpha map
        if self.project.dimension == '2D':
            alpha_map = self._compute_alpha_2D_map(cost_map)
            
        elif self.project.dimension == '3D':
            alpha_map = self._compute_alpha_3D_map(cost_map)
        
        # Increase cost from cells outside the domain. Retrieve the maximum
        # cost value from the domain and set it outside the domain.
        logical_test = self.domain.data_volume.astype('bool')
        max_value = np.where(logical_test, alpha_map, 0).max()
        logical_test = np.invert(logical_test)
        alpha_map = np.where(logical_test, max_value * 10, alpha_map)

        # Append the alpha map from the current iteration to the list
        self.maps['alpha'].append(alpha_map)
        return None

    def _compute_alpha_3D_map(self, cost_map: np.ndarray) -> np.ndarray:
        """
        Compute the alpha map: travel cost in the same direction as the
        gradient.
        """
        ### Option A
        # Vadose, bedrock, and phreatic zones are equivalent
        if self.model_parameters['sks']['mode'] == 'A':
            
            # The cost is depending from the elevation
            alpha_map = self._set_alpha_3D_from_elevation(cost_map)
        
        ### Option B
        # Vadose, and bedrock zones are equivalent. Phreatic zone is isotropic.
        if self.model_parameters['sks']['mode'] == 'B':

            # The cost is depending from the elevation
            alpha_map = self._set_alpha_3D_from_elevation(cost_map)
            
            # Phreatic zone: isotropic cost
            if (self.grid.nz > 1) and (self.domain.water_table is not None):
                logical_test = self.domain.get_subdomain('phreatic_zone')
                alpha_map = np.where(logical_test, cost_map, alpha_map)
            
        ### Option C
        # Vadose, and bedrock zones are multiplied by a F cost factor.
        # Phreatic zone is isotropic.
        if self.model_parameters['sks']['mode'] == 'C':
            
            # Vadose and bedrock zones: multiplied by a cost factor
            F = self.model_parameters['sks']['factors']['F']
            logical_test = self.domain.get_subdomain('vadose_zone')
            alpha_map = np.where(logical_test, cost_map * F, cost_map)
            
            # Phreatic zone: isotropic cost
            if self.domain.water_table is not None:
                logical_test = self.domain.get_subdomain('phreatic_zone')
                alpha_map = np.where(logical_test, cost_map, alpha_map)
        
        ### Option D
        # Vadose zone is multiplied by a F1 cost factor.
        # Bedrock vadose zone is multiplied by a F2 cost factor.
        # Phreatic zone is isotropic.
        if self.model_parameters['sks']['mode'] == 'D':
            
            # Retrieve factor values
            F1 = self.model_parameters['sks']['factors']['F1']
            F2 = self.model_parameters['sks']['factors']['F2']
            
            # Vadose zone without bedrock vadose zone:
            # multiplied by a f1 cost factor
            bedrock_vadose = self.domain.get_subdomain('bedrock_vadose')
            logical_test = np.logical_not(bedrock_vadose.astype(bool))
            alpha_map = np.where(logical_test, cost_map * F1, cost_map)
            
            # Bedrock vadose zone: multiplied by a f2 cost factor
            if self.domain.bedrock is not None:
                logical_test = self.domain.get_subdomain('bedrock_vadose')
                alpha_map = np.where(logical_test, cost_map * F2, alpha_map)
            
            # Phreatic zone: isotropic cost
            if self.domain.water_table is not None:
                logical_test = self.domain.get_subdomain('phreatic_zone')
                alpha_map = np.where(logical_test, cost_map, alpha_map)
                
        return alpha_map
    
    def _set_alpha_3D_from_elevation(self, cost_map: np.ndarray) -> np.ndarray:
        """"""
        X, Y, Z = self.grid.get_meshgrids()
        normalized_z = normalize_array(Z) + 1
        cost_map = cost_map * normalized_z
        return cost_map
    
    def _compute_alpha_2D_map(self, cost_map: np.ndarray) -> np.ndarray:
        """
        Compute the alpha map: travel cost in the same direction as the
        gradient.
        """
        ### Option X  # TODO
        # Conduits follow the bedrock gradient
        if self.domain._is_defined['bedrock']:
            alpha_map = self._set_alpha_2D_from_bedrock(cost_map)
         
        ### Option Y  # TODO
        # TODO
        
        ### Default option
        else:
            alpha_map = cost_map
        
        return alpha_map
    
    def _set_alpha_2D_from_bedrock(self, cost_map: np.ndarray) -> np.ndarray:
        """"""
        if (self.grid.nz == 1):
            bedrock = self.domain.bedrock.data_surface
            bedrock = bedrock.reshape(self.grid.shape)
            normalized_bedrock = normalize_array(bedrock) + 1
            alpha_map = cost_map * normalized_bedrock
        else:
            X, Y, Z = self.grid.get_meshgrids()
            normalized_z = normalize_array(Z) + 1
            alpha_map = cost_map * normalized_z
        return alpha_map
    
    @wp._logging()
    def _compute_beta_map(self) -> None:
        """
        Computes the beta map: travel cost perpendicular to the gradient. If
        beta is higher than alpha, conduits will follow the steepest gradient.
        If beta is lower than alpha, conduits will follow contours.
        """
        # Retrieve the costs maps from the current iteration
        ratio = self.model_parameters['sks']['costs']['ratio']
        cost_map = self.maps['cost'][self.iteration].copy()
        alpha_map = self.maps['alpha'][self.iteration].copy()
        
        # Calculate the default beta map
        beta_map = alpha_map / ratio
        
        # Select the appropriate method to alter the beta map
        if self.project.dimension == '2D':
            beta_map = self._compute_beta_2D_map(cost_map,
                                                 alpha_map,
                                                 beta_map)
            
        elif self.project.dimension == '3D':
            beta_map = self._compute_beta_3D_map(cost_map,
                                                 alpha_map,
                                                 beta_map)
        
        # Append the beta map from the current iteration to the list
        self.maps['beta'].append(beta_map)
        return None

    def _compute_beta_2D_map(self,
                             cost_map: np.ndarray,
                             alpha_map: np.ndarray,
                             beta_map: np.ndarray,
                             ) -> np.ndarray:
        """"""
        return beta_map

    def _compute_beta_3D_map(self,
                             cost_map: np.ndarray,
                             alpha_map: np.ndarray,
                             beta_map: np.ndarray,
                             ) -> np.ndarray:
        """"""
        ### Option A
        # Vadose zone = Bedrock zone = Phreatic zone
        pass
            
        ### Options B / C / D
        # Vadose zone = Bedrock zone ≠ Phreatic zone
        # Isotropic fast marching in phreatic zone
        if self.model_parameters['sks']['mode'] in ['B', 'C', 'D']:
            if self.domain._is_defined['water_table']:
                test = self.domain.get_subdomain('phreatic_zone')
                beta_map = np.where(test, cost_map, beta_map)
        return beta_map

    #########################
    ##### GRADIENT MAPS #####
    #########################
    
    @wp._logging()
    def _compute_gradient(self) -> None:
        """Compute the gradient in the x, y and z-axis."""
        
        ### Default situation
        grad_x = np.full(self.grid.shape, 0, dtype=np.float32)
        grad_y = np.full(self.grid.shape, 0, dtype=np.float32)
        grad_z = np.full(self.grid.shape, 1, dtype=np.float32)
        grad = [grad_x, grad_y, grad_z]
        
        # Select the appropriate method to alter the gradient maps
        if self.project.dimension == '2D':
            grad = self._compute_gradient_2D(*grad)
            
        elif self.project.dimension == '3D':
            grad = self._compute_gradient_3D(*grad)
            
        # Append the gradient maps from the current iteration to the list
        self.maps['gradient'].append(grad)
        return None
    
    def _compute_gradient_2D(self,
                             grad_x: np.ndarray,
                             grad_y: np.ndarray,
                             grad_z: np.ndarray
                             ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Compute gradient maps in the case of a 2D project."""
        
        grad = (grad_x, grad_y, grad_z)
        
        ### Option ?
        if self.domain._is_defined['bedrock']:
            grad = self._set_gradient_2D_from_bedrock(*grad)
        
        return grad
        
    def _compute_gradient_3D(self,
                             grad_x: np.ndarray,
                             grad_y: np.ndarray,
                             grad_z: np.ndarray
                             ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Compute gradient maps in the case of a 3D project."""
        ### Option A
        # Vadose zone = Bedrock zone = Phreatic zone
        if self.model_parameters['sks']['mode'] in ['A']:
            pass
    
        ### Options B / C / D
        elif self.model_parameters['sks']['mode'] in ['B', 'C', 'D']:
            
            if self.domain._is_defined['water_table']:
                subdomain = self.domain.get_subdomain('phreatic_zone')
                test_subdomain = (subdomain == 1)
                grad_x = np.where(test_subdomain, 1, grad_x)
                grad_y = np.where(test_subdomain, 1, grad_y)
                grad_z = np.where(test_subdomain, 1, grad_z)
        
            if self.domain._is_defined['bedrock']:
                grad_x, grad_y, grad_z = (
                    self._set_gradient_3D_from_bedrock(grad_x,
                                                       grad_y,
                                                       grad_z)
                )
        grad = (grad_x, grad_y, grad_z)
        return grad
    
    def _set_gradient_2D_from_bedrock(self,
                                      grad_x: np.ndarray,
                                      grad_y: np.ndarray,
                                      grad_z: np.ndarray,
                                      ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Calculate the gradient from the bedrock for a 2D project."""
        
        # Test which dimension equals 1
        test_grid_nx = (self.grid.nx == 1)
        test_grid_ny = (self.grid.ny == 1)
        test_grid_nz = (self.grid.nz == 1)
        
        # X-axis
        if test_grid_nx:
            bedrock = np.roll(self.domain.bedrock.data_volume, 1, axis=2)
            bedrock[:, :, 0] = 1
            gradient_y, gradient_z = np.gradient(bedrock,
                                                 self.grid.dy,
                                                 self.grid.dz,
                                                 axis=(1, 2))
            gradient_x = np.full_like(gradient_y, 0)
        
        # Y-axis
        elif test_grid_ny:
            bedrock = np.roll(self.domain.bedrock.data_volume, 1, axis=2)
            bedrock[:, :, 0] = 1
            gradient_x, gradient_z = np.gradient(bedrock,
                                                 self.grid.dx,
                                                 self.grid.dz,
                                                 axis=(0, 2))
            gradient_y = np.full_like(gradient_x, 0)
        
        # Z-axis
        elif test_grid_nz:
            bedrock = self.domain.bedrock.data_surface
            gradient_x, gradient_y = np.gradient(bedrock,
                                                 self.grid.dx,
                                                 self.grid.dy,
                                                 axis=(0, 1))
            # Correct shape from N2 to N3
            gradient_x = gradient_x.reshape(self.grid.shape)
            gradient_y = gradient_y.reshape(self.grid.shape)
            gradient_z = np.full_like(gradient_x, 0)
        
        # Apply the new gradient only in domain
        test_domain = (self.domain.get_subdomain('domain') == 1)
        grad_x = np.where(test_domain, gradient_x, grad_x)
        grad_y = np.where(test_domain, gradient_y, grad_y)
        grad_z = np.where(test_domain, gradient_z, grad_z)
        grad = (grad_x, grad_y, grad_z)
        return grad
    
    def _set_gradient_3D_from_bedrock(self,
                                      grad_x: np.ndarray,
                                      grad_y: np.ndarray,
                                      grad_z: np.ndarray
                                      ) -> np.ndarray:
        """Calculate the gradient from the bedrock for a 3D project."""
        # On commence par calculer le gradient en faisant gaffe
        # à l'intervertion x y
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
        
        test_domain = self.domain.get_subdomain('bedrock_vadose')
        test_domain = test_domain.astype('bool')
        args = np.argwhere(test_domain)
        
        if len(args) != 0:
            
            i, j, k = zip(*args)

            grad_x[i, j, k] = -vx[i, j]
            grad_y[i, j, k] = -vy[i, j]
            grad_z[i, j, k] = -vz[i, j]
        
        out = (grad_x, grad_y, grad_z)
        return out
            
    def _compute_riemann_metric(self) -> None:
        """
        Compute the riemann metric: Define the Riemannian metric needed as
        input for the anisotropic fast marching.

        TODO : à terminer
        """
        # Give the parameters for the fast marching algorithm
        alpha = self.maps['alpha'][self.iteration]
        beta = self.maps['beta'][self.iteration]
        grad_x, grad_y, grad_z = self.maps['gradient'][self.iteration]

        self.fmm['riemannMetric'] = agd.Metrics.Riemann.needle(
            [grad_x, grad_y, grad_z],
            alpha,
            beta
        )
        return None

    ### 2.1.2
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
        verbosity = self.model_parameters['verbosity']['agd']
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
    @wp._logging()
    def _compute_karst_map(self):
        """
        Compute the karst map based on the paths from agd-hfm.
        Array of all zeros, with ones in cells containing a karst conduit.
        """
        ### Loop over conduit paths generated by fast marching:
        for path in self.fmm['fastMarchingOutput']['geodesics']: # loop over conduit paths in this iteration (there is one path from each inlet)
            merge = False                                        # reset indicator for whether this conduit has merged with an existing conduit
            for p in range(path.shape[1]):                       # loop over points making up this conduit path
                point = path[:, p]                                # get coordinates of current point
                [[ix, iy, iz], error] = self.fmm['fastMarching'].IndexFromPoint(point) #convert to coordinates to indices, /!\ returning iy first then ix TODO ???
                # Place nodes and links:
                if np.isnan(self.maps['nodes'][ix, iy, iz]):                                    #if there is no existing conduit node here
                    if ~np.isnan(self.maps['outlets'][ix, iy, iz]):                              #if there is an outlet here (cell value is not nan)
                        outlet = self._outlets.iloc[int(self.maps['outlets'][ix, iy, iz])]         #get the outlet coordinates using the ID in the outlets map
                        self.vectors['nodes'][self.vectors['n']] = [outlet.x, outlet.y, outlet.z, 'outlet']           #add a node at the outlet coordinates (with the node type for SWMM)
                        self.maps['nodes'][ix, iy, iz] = self.vectors['n']                                   #update node map with node index
                        if p > 0:                                                           #if this is not the first point (i.e. the inlet) in the current path
                            if merge is False:                                               #if this conduit has not merged with an existing conduit
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
                            # possible improvement: if the next point on the path is on an existing point, skip the current point.
                            self.vectors['nodes'][self.vectors['n']] = [point[0], point[1], point[2], 'junction']            #add a junction node here (with the node type for SWMM)
                            self.maps['nodes'][ix, iy, iz] = self.vectors['n']                               #update node map with node index
                            if merge is False:                                              #if this conduit has not merged with an existing conduit
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
                    if merge is True:                                                       #if this conduit has already merged into an existing conduit
                        pass                                                                 #skip this node (there is already a node here)
                    elif n_existing == self.vectors['n'] - 1:                                            #if existing index is only one less than next node to be added index, this is a duplicate node and can be skipped
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
        
        # Store the edges in a DataFrame
        self.vectors['edges_'] = pd.DataFrame.from_dict(
            self.vectors['edges'],
            orient='index',
            columns=['n1', 'n2'],
            dtype=int,
        )

        # Store the nodes in a DataFrame
        self.vectors['nodes_'] = pd.DataFrame.from_dict(
            self.vectors['nodes'],
            orient='index',
            columns=['x', 'y', 'z', 'type'],
        )
        
        # Indicate if nodes are in the vadose zone
        XYZ = self.vectors['nodes_'][['x', 'y', 'z']].to_numpy()
        I, J, K = self.grid.get_indices(XYZ)
        vadose_zone = self.domain.get_subdomain('vadose_zone')
        self.vectors['nodes_']['vadose'] = vadose_zone[I, J, K]
        
        return None
    
    def _voxelize_karst_network(self) -> None:
        """
        Using the calculated nodes and edges, voxelize the segments in a 3D
        array.
        """
        
        # Gets karst map from previous iteration
        # except for the very first iteration
        # TODO ????
        if self.iteration == 0:
            karst_map_0 = np.zeros(self.grid.shape)
            self.maps['karst'].append(karst_map_0)
        else:
            karst_map_copy = self.maps['karst'][self.iteration - 1].copy()
            self.maps['karst'].append(karst_map_copy)
        ############### ????

        ### Retrieve edges and nodes
        edges = self.vectors['edges_'].copy()
        nodes = self.vectors['nodes_'][['x', 'y', 'z']].copy()
        edges = edges.join(nodes, on='n1').join(nodes, on='n2', rsuffix='_')
        
        ### For each segment, create midpoints
        ndx = np.ceil(
            np.abs(edges.x.to_numpy() - edges.x_.to_numpy())
            / self.grid.dx + 1
        ).astype('int')
        ndy = np.ceil(
            np.abs(edges.y.to_numpy() - edges.y_.to_numpy())
            / self.grid.dy + 1
        ).astype('int')
        ndz = np.ceil(
            np.abs(edges.z.to_numpy() - edges.z_.to_numpy())
            / self.grid.dz + 1
        ).astype('int')
        edges['n'] = np.vstack([ndx, ndy, ndz]).T.max(axis=1)

        edges['X'] = edges.apply(
            lambda x: np.linspace(x.x, x.x_, int(x.n)),
            axis=1,
        )
        edges['Y'] = edges.apply(
            lambda x: np.linspace(x.y, x.y_, x.n),
            axis=1,
        )
        edges['Z'] = edges.apply(
            lambda x: np.linspace(x.z, x.z_, x.n),
            axis=1,
        )
        
        ### Retrieve all the coordinates
        X = [x for x_ in edges.X.tolist() for x in x_]
        Y = [y for y_ in edges.Y.tolist() for y in y_]
        Z = [z for z_ in edges.Z.tolist() for z in z_]
        
        ### Retrieve indices and voxelize
        coords = np.vstack([X, Y, Z])
        i, j, k = self.grid.get_indices(coords.T)
        new_karst = np.zeros_like(self.grid.data_volume)
        new_karst[i, j, k] = 1
        self.maps['karst'][self.iteration] = new_karst
        return None
    
    def _export_results(self):
        """
        TODO
        """
        if self.export_settings == {}:
            # TODO
            pass
        
        path = self.simulation_path + 'results.pickle'
        
        results = {
            'maps': self.maps.copy(),
            'vectors': self.vectors.copy(),
            'inlets': self.inlets.copy(),
            'outlets': self.outlets.copy(),
            'domain': self.domain,
            'geology': self.geology,
            'faults': self.faults,
            'fractures': self.fractures,
            'parameters': self.model_parameters.copy()
        }
        
        with open(path, 'wb') as handle:
            pickle.dump(results, handle, protocol=pickle.HIGHEST_PROTOCOL)
        return None
