"""
This module contains a class modeling the karstic network generator tool.
"""

### Internal dependencies
import os
import sys
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
from .domain import Domain, Delimitation, Topography, Bedrock, WaterLevel
from .geologic_features import Geology, Faults
from .fracturation import Fractures
from .points import PointGenerator
from pykasso._utils.array import normalize_array, transform_array_where

### Typing
from pykasso._typing import Project

### Access 'Application' instance memory
# app = sys.modules['pykasso.core.application']

##### Module variables ########################################
this = sys.modules[__name__]

# Define default fast-marching costs
this.default_fmm_costs = {
    # 'out': 0.999,
    'out': 10,  # TODO
    'aquifer': 0.4,
    'aquiclude': 0.8,
    'beddings': 0.35,
    'faults': 0.2,
    'fractures': 0.2,
    'karst': 0.1,
    'conduits': 0.1,
    'ratio': 0.5,
}

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
    
    def __init__(self, project: Project) -> None:
        """_summary_

        Returns
        -------
        _type_
            _description_
        """
        ### Initialization
        self.project = project
        self.grid = project.grid
        self.domain = None
        self.geology = None
        self.faults = None
        self.fractures = None
        self.inlets = None
        self.outlets = None
        self.conceptual_model = None
        self.conceptual_model_table = None
        self.rng = {}

    def generate(self, model_parameters: dict = {},
                 export_settings: dict = {}) -> None:
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
        Initialize and configure the basic settings
        """
        ### Load model parameters
        self._load_model_parameters()
        
        ### Create simulation directory
        self._create_sim_dir()
        
        ### Check verbosity levels
        self._check_verbosity()
        
        ### Update log
        self._update_log()
    
        return None
    
    def _load_model_parameters(self) -> None:
        """
        Load the model parameters
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
        
        return None
    
    def _create_sim_dir(self) -> None:
        """
        Create simulation directory
        """
        self.project.n_simulations += 1
        outputs_dir = self.project.core['paths']['outputs_dir']
        sim_dir = 'simulation_{}\\'.format(self.project.n_simulations)
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
        Update log
        """
        ### LOGGING ###
        
        # Sets logging level
        level = self.model_parameters['verbosity']['logging']
        if level > 4:
            level = 4
        if level < 0:
            level = 0
        logging_level = self.project._logging_levels[level]
        self.logger = logging.getLogger('♠').setLevel(logging_level)
        
        # Prints current simulation number
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
         1. Define the model extension : ``self.domain`` ;
         2. Build the geological features : ``self.geology``, ``self.faults``,
         ``self.inlets``, ``self.outlets``, ``self.fractures``,
         ``self.conceptual_model``, ``self.conceptual_model_table``;
         3. Prepare the fast-marching method variables.
         
        If a geological feature has not been defined in the settings, calling
        its attribute will return ``None``.
        """
        self._build_domain()             # 1 - Constructs the domain
        self._build_model()              # 2 - Constructs the model
        self._construct_fmm_variables()  # 3 - Constructs fmm features
        return None
    
    ##### 1 - DOMAIN #########################
    @wp._parameters_validation('domain', 'optional')
    @wp._memoize('domain')
    @wp._logging()
    def _build_domain(self) -> None:
        """
        Build the domain.
        """
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
    
    @wp._logging()
    def _build_domain_delimitation(self):
        """
        Build the delimitation.
        """
        delimitation = self.model_parameters['domain']['delimitation']
        test_a = isinstance(delimitation, (str))
        if not (test_a and (delimitation == '')):
            if isinstance(delimitation, (str)):
                delimitation = np.genfromtxt(delimitation).tolist()
                self.model_parameters['domain']['delimitation'] = delimitation
            instance = Delimitation(vertices=delimitation, grid=self.grid)
            return instance
        else:
            return None

    @wp._logging()
    def _build_domain_topography(self):
        """
        Build the topography.
        """
        topography = self.model_parameters['domain']['topography']
        test_a = isinstance(topography, (str))
        if not (test_a and (topography == '')):
            instance = Topography()
            instance.set_data(data=topography, grid=self.grid)
            instance._surface_to_volume('<=', self.grid)
            return instance
        else:
            return None
        
    @wp._logging()
    def _build_domain_bedrock(self):
        """
        Build the bedrock elevation.
        """
        bedrock = self.model_parameters['domain']['bedrock']
        test_a = isinstance(bedrock, (str))
        if not (test_a and (bedrock == '')):
            instance = Bedrock()
            instance.set_data(data=bedrock, grid=self.grid)
            instance._surface_to_volume('<=', self.grid)
            return instance
        else:
            return None
        
    @wp._logging()
    def _build_domain_water_level(self):
        """
        Build the water level elevation.
        """
        water_level = self.model_parameters['domain']['water_level']
        test_a = isinstance(water_level, (str))
        if not (test_a and (water_level == '')):
            instance = WaterLevel()
            instance.set_data(data=water_level, grid=self.grid)
            instance._surface_to_volume('<=', self.grid)
            return instance
        else:
            return None

    ##### 2 - MODEL #########################
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
        self._build_model_fractures()
        
        ### Constructs the conceptual model
        self._build_conceptual_model()
        return None
    
    def _set_rng(self, attribute):
        """
        Populates the rng seeds dictionary with the defined seed attribute.
        """
        
        # Generates a random seed when the seed equals 0 (default setting)
        if self.model_parameters[attribute]['seed'] == 0:
            seed = np.random.default_rng().integers(low=0, high=10**6)
            self.model_parameters[attribute]['seed'] = seed
            
        # Sets the seed
        seed = self.model_parameters[attribute]['seed']
        self.rng[attribute] = np.random.default_rng(seed)
        return None
    
    @wp._parameters_validation('sks', 'optional')
    def _build_model_parameters(self) -> None:
        """
        Build characteristic model parameters.
        """
        # Define the main seed
        self._set_rng('sks')
        return None
    
    @wp._parameters_validation('geology', 'optional')
    @wp._memoize('geology')
    @wp._logging()
    def _build_model_geology(self) -> None:
        """
        Build the geology.
        """
        geology = self.model_parameters['geology']['data']
        test_a = isinstance(geology, (str))
        self.geology = Geology()
        if not (test_a and (geology == '')):
            axis = self.model_parameters['geology']['axis']
            self.geology.set_data(data=geology, grid=self.grid, axis=axis)
        else:
            self.geology.data_volume = (
                self.geology._set_data_full_3D(grid=self.grid, value=1)
            )
        costs = self.model_parameters['geology']['costs']
        self.geology._set_costs(costs)
        self.geology._compute_statistics(self.grid)
        return None
    
    @wp._parameters_validation('faults', 'optional')
    @wp._memoize('faults')
    @wp._logging()
    def _build_model_faults(self) -> None:
        """
        Build the faults.
        """
        faults = self.model_parameters['faults']['data']
        test_a = isinstance(faults, (str))
        if not (test_a and (faults == '')):
            self.faults = Faults()
            axis = self.model_parameters['faults']['axis']
            self.faults.set_data(data=faults, grid=self.grid, axis=axis)
            costs = self.model_parameters['faults']['costs']
            self.faults._set_costs(costs)
            self.faults._compute_statistics(self.grid)
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
            subdomain=self.model_parameters[kind]['subdomain'],
            domain=self.domain,
            geology=self.geology,
            geologic_ids=self.model_parameters[kind]['geology']
        )

        ### Gets existing points

        # Loads points if needed
        logical_test_1 = isinstance(self.model_parameters[kind]['data'], (str))
        logical_test_2 = not (self.model_parameters[kind]['data'] == '')
        if logical_test_1 and logical_test_2:
            path = self.model_parameters[kind]['data']
            data = np.genfromtxt(path)
            if len(data.shape) == 1:
                data = np.array([data])
            self.model_parameters[kind]['data'] = data
        
        ### Inspects validity of points
        points = self.model_parameters[kind]['data']
        
        # 2D points # TODO - logging ?
        points_2D = [point for point in points if len(point) == 2]
        validated_points_2D = [point for point in points_2D
                               if point_manager._is_point_valid(point)]
        validated_points_3D = [
            point_manager._generate_3D_point_from_2D_point(point)
            for point in validated_points_2D
        ]
        
        # 3D points # TODO - logging ?
        points_3D = [point for point in points if len(point) == 3]
        validated_points_3D += (
            [point for point in points_3D
             if point_manager._is_point_valid(point)]
        )

        diff = len(points) - len(validated_points_3D)
        if diff > 0:
            # TODO - LOG - VERBOSITY
            # TODO - log name is not correct
            msg = ('{}/{} {} have been discarded because out of domain.'
                   .format(diff,
                           len(self.model_parameters[kind]['data']), kind))
            self.logger.warning(msg)
        self.model_parameters[kind]['data'] = validated_points_3D

        ### Get new points according to the right case
        n = self.model_parameters[kind]['number']
        data = self.model_parameters[kind]['data']
        # Case 1 - No points declared
        if (data == '') or (data == []):
            points = point_manager._generate_points(size=n)
        # Case 2 - More points required than provided
        elif (n > len(data)):
            n_points = n - len(data)
            points = np.append(np.array(data), point_manager
                               ._generate_points(n_points), axis=0)
        # Case 3 - Less points required than provided
        elif (n < len(data)):
            points = self.rng['sks'].choice(data, n, replace=False)
        # Case 4 - Points required equals points declared
        else:
            points = data

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
    
    @wp._parameters_validation('fractures', 'optional')
    @wp._memoize('fractures')
    @wp._logging()
    def _build_model_fractures(self) -> None:
        """
        Build the fractures.
        """
        self._set_rng('fractures')
        fractures = self.model_parameters['fractures']['data']
        test_a = isinstance(fractures, (str))
        test_b = ('settings' in self.model_parameters['fractures'])
        if (not (test_a and (fractures == ''))) or test_b:
            self.fractures = Fractures(rng=self.rng['fractures'])
            
            # Generate fractures families
            if 'settings' in self.model_parameters['fractures']:
                
                frac_settings = self.model_parameters['fractures']['settings']
                for frac_name, frac_settings in frac_settings.items():
                    self.fractures.generate_fracture_family(frac_name,
                                                            self.grid,
                                                            **frac_settings)
                # TODO - 'superposition' -> into a parameter
                self.fractures.generate_model('superposition', self.grid)
            # Load data
            else:
                axis = self.model_parameters['fractures']['axis']
                self.fractures.set_data(fractures, self.grid, axis)
                costs = self.model_parameters['fractures']['costs']
                self.fractures._set_costs(costs)
            
            self.fractures._compute_statistics(self.grid)
        return None
    
    @wp._logging()
    def _build_conceptual_model(self) -> None:
        """Build the conceptual model.
        
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
        geology_items = [(100 + i, 'Geology', id, cost) for (i, (id, cost))
                         in enumerate(self.geology.costs.items())]
        for (id, feature, id_, cost) in geology_items:
            conceptual_model = np.where(self.geology.data_volume == id_,
                                        id,
                                        conceptual_model)
        conceptual_model_table += geology_items
        
        # 200 - Fractures
        if self.fractures is not None:
            fractures_items = [(200 + i, 'Fractures', id, cost)
                               for (i, (id, cost))
                               in enumerate(self.fractures.costs.items())]
            for (id, feature, id_, cost) in fractures_items:
                conceptual_model = np.where(self.fractures.data_volume == id_,
                                            id,
                                            conceptual_model)
            conceptual_model_table += fractures_items
            
        # 300 - Faults
        if self.faults is not None:
            faults_items = [(300 + i, 'Faults', id, cost) for (i, (id, cost))
                            in enumerate(self.faults.costs.items())]
            for (id, feature, id_, cost) in faults_items:
                conceptual_model = np.where(self.faults.data_volume == id_,
                                            id,
                                            conceptual_model)
            conceptual_model_table += faults_items
        
        # 0 - Out
        cost_out = self.model_parameters['sks']['costs']['out']
        out_item = [(0, 'Out', np.nan, cost_out)]
        conceptual_model = np.where(self.domain.data_volume == 0,
                                    0,
                                    conceptual_model)
        conceptual_model_table += out_item
        conceptual_model_table = pd.DataFrame(conceptual_model_table, columns=['id', 'feature', 'id-feature', 'cost']).set_index('id').sort_values('id')

        self.conceptual_model = conceptual_model
        self.conceptual_model_table = conceptual_model_table
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
        inlets_per_outlet = self.model_parameters['inlets']['per_outlet']

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
                test = self._outlets['outlet_iteration'] == outlet_iteration
                outlets_current = self._outlets[test]

                # Gets the inlets assigned to the current inlet iteration
                test = self._inlets['outlet_key'] == outlet_iteration
                inlets_current = self._inlets[test]
                test = inlets_current['inlet_iteration'] == inlet_iteration
                inlets_current_iteration = inlets_current[test]
                
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
        
        # print(self.iterations)
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
            'gradient': [], # bebdrock gradient
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
        Computes the karst network according to the parameters. Must be called
        after ``build()`` method. This method will :
        1. Computes conduits for each generation and stores nodes & edges for
        network.
        2. Stores all the relevant data for this network in specific
        dictionaries.
        3. Exports the state of the project, this file could be read by the
        'analysis' and 'visualization' sub-packages.
        """
        self._compute_karst_network()
        self._export_results()
        self.project._export_project_file()
        return None
    
    def _compute_karst_network(self):
        self.logger = logging.getLogger("fmm.modelisation")
        self.logger.info("Computing karst network")
        
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
            self.logger.info(msg)
        return None

    ### 2.1.1 Iso- and anisotropic case
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
                         self.model_parameters['sks']['costs']['conduits'],
                         self.maps['cost'][self.iteration])
            )
            
        return None
    
    ### 2.1.4 Anisotropic case
    @wp._logging()
    def _compute_alpha_map(self) -> None:
        """
        Computes the alpha map: travel cost in the same direction as the
        gradient.
        """
       
        cost_map = self.maps['cost'][self.iteration].copy()
        
        ### Option A
        if self.model_parameters['sks']['mode'] == 'A':
            # Vadose zone = Phreatic zone = Bedrock zone
            alpha_map = self._set_alpha_from_elevation(cost_map)
            
            # Bedrock
            if (self.grid.nz > 1) and (self.domain.bedrock is not None):
                alpha_map = transform_array_where(
                    alpha_map,
                    alpha_map.max(),
                    self.domain.bedrock.data_volume)
        
        ### Option B
        if self.model_parameters['sks']['mode'] == 'B':
            # Vadose zone = Bedrock zone
            alpha_map = self._set_alpha_from_elevation(cost_map)
            # Bedrock
            if (self.grid.nz > 1) and (self.domain.bedrock is not None):
                alpha_map = transform_array_where(
                    alpha_map,
                    alpha_map.max(),
                    self.domain.bedrock.data_volume)
            # Phreatic zone
            if (self.grid.nz > 1) and (self.domain.water_level is not None):
                alpha_map = transform_array_where(
                    alpha_map,
                    cost_map,
                    self.domain.get_subdomain('phreatic_zone'))
            
        ### Option C
        if self.model_parameters['sks']['mode'] == 'C':
            # Vadose zone = Bedrock zone
            factor = self.model_parameters['sks']['factors']['F']
            alpha_map = transform_array_where(
                cost_map,
                cost_map * factor,
                self.domain.get_subdomain('vadose_zone'))
            # Bedrock
            if self.domain.bedrock is not None:
                alpha_map = transform_array_where(
                    alpha_map,
                    alpha_map.max() * 10,
                    self.domain.bedrock.data_volume)
            # Phreatic zone
            if self.domain.water_level is not None:
                alpha_map = transform_array_where(
                    alpha_map,
                    cost_map,
                    self.domain.get_subdomain('phreatic_zone'))
        
        ### Option D
        if self.model_parameters['sks']['mode'] == 'D':
            # Vadose zone
            domain = np.logical_and(
                np.logical_not(self.domain.get_subdomain('bedrock_vadose')),
                self.domain.get_subdomain('vadose_zone')
            )
            factor_01 = self.model_parameters['sks']['factors']['F1']
            alpha_map = transform_array_where(
                cost_map,
                cost_map * factor_01,
                domain)
            
            if self.domain.bedrock is not None:
                # Bedrock vadose zone
                factor_02 = self.model_parameters['sks']['factors']['F2']
                alpha_map = transform_array_where(
                    alpha_map,
                    cost_map * factor_02,
                    self.domain.get_subdomain('bedrock_vadose'))
                
                # Bedrock
                alpha_map = transform_array_where(
                    alpha_map,
                    alpha_map.max() * 10,
                    self.domain.bedrock.data_volume)
                
            # Phreatic zone
            if self.domain.water_level is not None:
                alpha_map = transform_array_where(
                    alpha_map,
                    cost_map,
                    self.domain.get_subdomain('phreatic_zone'))
            
        # Out of domain
        if self.domain.topography is not None:
            domain = np.logical_not(self.domain.topography.data_volume)
            alpha_map = transform_array_where(
                alpha_map,
                alpha_map.max() * 10,
                domain)
        
        self.maps['alpha'].append(alpha_map)
        return None
    
    def _set_alpha_from_elevation(self, cost_map: np.ndarray) -> np.ndarray:
        if (self.grid.nz == 1) and (self.domain.bedrock is not None):
            # print('flag')
            bedrock = self.domain.bedrock.data_surface.reshape(self.grid.shape)
            normalized_bedrock = normalize_array(bedrock) + 1
            cost_map = cost_map * normalized_bedrock
        else:
            X, Y, Z = self.grid.get_meshgrids()
            normalized_z = normalize_array(Z) + 1
            cost_map = cost_map * normalized_z
        return cost_map
    
    ### 2.1.5 Anisotropic case
    @wp._logging()
    def _compute_beta_map(self) -> None:
        """
        Computes the beta map: travel cost perpendicular to the gradient. If
        beta is higher than alpha, conduits will follow the steepest gradient.
        If beta is lower than alpha, conduits will follow contours.
        """
        ratio = self.model_parameters['sks']['costs']['ratio']
        cost_map = self.maps['cost'][self.iteration].copy()
        alpha_map = self.maps['alpha'][self.iteration].copy()
        
        ### Option A & default situation
        # Vadose zone = Phreatic zone = Bedrock zone
        beta_map = alpha_map / ratio
            
        ### Options B / C / D
        if self.model_parameters['sks']['mode'] in ['B', 'C', 'D']:
            # Phreatic zone
            if self.domain.water_level is not None:
                beta_map = transform_array_where(
                    beta_map,
                    cost_map,
                    self.domain.get_subdomain('phreatic_zone'))
        
        self.maps['beta'].append(beta_map)
        return None

    ### 2.1.6 Anisotropic case
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
        if self.model_parameters['sks']['mode'] in ['B', 'C', 'D']:
            
            if self.domain.water_level is not None:
                subdomain = self.domain.get_subdomain('phreatic_zone')
                test_subdomain = (subdomain == 1)
                grad_x = np.where(test_subdomain, 1, grad_x)
                grad_y = np.where(test_subdomain, 1, grad_y)
                grad_z = np.where(test_subdomain, 1, grad_z)
        
            if self.domain.bedrock is not None:
                grad_x, grad_y, grad_z = (
                    self._set_gradient_from_bedrock(grad_x,
                                                    grad_y,
                                                    grad_z)
                )
        
        # Saves the gradient
        self.maps['gradient'].append((grad_x, grad_y, grad_z))
        
        # Sets the needle
        alpha = self.maps['alpha'][self.iteration]
        beta = self.maps['beta'][self.iteration]
        
        print(self.model_parameters['sks']['mode'])
        # from pykasso.visualisation.visualiser import Visualiser
        # Visualiser.show_array(grad_x)
        # Visualiser.show_array(grad_y)
        # Visualiser.show_array(grad_z)
        # Visualiser.show_array(alpha)
        # Visualiser.show_array(beta)
        
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
        test_domain = (self.domain.get_subdomain('bedrock_vadose') == 1)
        
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
                # Corrects shape from N2 to N3
                gradient_x = gradient_x.reshape(self.grid.shape)
                gradient_y = gradient_y.reshape(self.grid.shape)
                gradient_z = np.full_like(gradient_x, 0)
                
            test_domain = (self.domain.get_subdomain('domain') == 1)
            grad_x = np.where(test_domain, gradient_x, grad_x)
            grad_y = np.where(test_domain, gradient_y, grad_y)
            grad_z = np.where(test_domain, gradient_z, grad_z)
            
        else:
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
            
            args = np.argwhere(test_domain)
            
            if len(args) != 0:
                
                i, j, k = zip(*args)

                grad_x[i, j, k] = -vx[i, j]
                grad_y[i, j, k] = -vy[i, j]
                grad_z[i, j, k] = -vz[i, j]
        
        out = (grad_x, grad_y, grad_z)
        return out

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
        
        # dx = self.grid.dx
        # dy = self.grid.dy
        # dz = self.grid.dz
        # n_grid = np.sqrt(dx**2 + dy**2 + dz**2)
        
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

        ### Retrieves the points
        points = []
        for (i1, i2) in edges:
            
            # Retrieves point coordinates
            point1 = nodes[i1]
            point2 = nodes[i2]
            x1, y1, z1 = point1
            x2, y2, z2 = point2
            
            # Calculates the number of segments to create
            ndx = np.ceil(np.abs(x2 - x1) / self.grid.dx).astype('int')
            ndy = np.ceil(np.abs(y2 - y1) / self.grid.dy).astype('int')
            ndz = np.ceil(np.abs(z2 - z1) / self.grid.dz).astype('int')
            n = max(ndx, ndy, ndz)
            if n == 1:
                n = 2

            # Calculates coordinates of points
            x = np.linspace(x1, x2, n)
            y = np.linspace(y1, y2, n)
            z = np.linspace(z1, z2, n)
            
            # Saves points
            points.extend(list(zip(x, y, z)))

        # Gets the indices
        i, j, k = self.grid.get_indices(points)

        # Calculates the new karst array
        new_karst = np.zeros_like(self.grid.data_volume)
        new_karst[i, j, k] = 1
        
        # import pykasso.visualization as pkv
        # pkv.show_array(new_karst, ghost=True)
        self.maps['karst'][self.iteration] = new_karst
        
        return None
    
    def _export_results(self):
        """
        TODO
        """
        path = self.simulation_path + 'results.pickle'
        with open(path, 'wb') as handle:
            results = {
                'maps': self.maps.copy(),
                'vectors': self.vectors.copy(),
                'inlets': self.inlets.copy(),
                'outlets': self.outlets.copy(),
                'domain': self.domain,
                'geology': self.geology,
                'faults': self.faults,
                'fractures': self.fractures,
            }
            pickle.dump(results, handle, protocol=pickle.HIGHEST_PROTOCOL)
        return None
