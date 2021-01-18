"""
TODO :
- retirer geology via image ???
- _set_points
- import methods : faire le point
"""

"""
pyKasso
=======

pyKasso is an object-oriented library intented for karst network simulation.

License
-------
Released under the MIT license:
   Copyright (C) 2020 Univertsity of Neuchâtel - CHYN
   François Miville <francois.miville@unine.ch>
   Chloé Fandel     <cfandel@email.arizona.edu>
   Philippe Renard  <philippe.renard@unine.ch>
"""

from .grid           import Grid
from .polygon        import Polygon
from .geologymanager import GeologyManager
from .pointmanager   import PointManager

import os
import sys

import yaml
import numpy    as np
import numpy.ma as ma
import pandas   as pd

class SKS():
    """
    A super-class to manage all the data class and to simulate karst networks.
    """

    def __init__(self, yaml_settings_file = None, rand_seed = None):
        """
        Construct a SKS class according to the specified settings datafile.

        Parameters
        ----------
        yaml_settings_file : string
            YAML settings file location.

        Examples
        --------
            >>> catchment = pk.SKS()
            >>> catchment = pk.SKS('exemple.yaml')
        """
        if yaml_settings_file is None:
            yaml_settings_file = os.path.dirname(os.path.abspath(__file__)) + '/' + 'default_files/settings.yaml'

        try:
            with open(yaml_settings_file, 'r') as stream:
                try:
                    settings = yaml.safe_load(stream)
                except yaml.YAMLError as exc:
                    print(exc)
                    raise
        except:
            print("/!\\ Error : unable to read the YAML settings file!")
            raise

        self.settings = settings
        
        #Set travel cost through each geologic formation
        geology_cost = []
        for elem in self.settings['geology_cost']:
            geology_cost.append(self.settings[elem])
        self.settings['geology_cost'] = geology_cost
        
        if self.settings['rand_seed'] > 0:
            np.random.seed(self.settings['rand_seed'])
            
        if rand_seed != None:
            np.random.seed(rand_seed)

        self._get_reference_statistics()
        self._load_data()
        self.karst_simulations = []
        
    ################
    # get_ methods #
    ################

    def get_x0(self):
        return self.settings['x0']

    def get_y0(self):
        return self.settings['y0']

    def get_z0(self):
        return self.settings['z0']

    def get_nx(self):
        return self.settings['nx']

    def get_ny(self):
        return self.settings['ny']

    def get_nz(self):
        return self.settings['nz']

    def get_dx(self):
        return self.settings['dx']

    def get_dy(self):
        return self.settings['dy']

    def get_dz(self):
        return self.settings['dz']
        
    def get_data_has_polygon(self):
        return self.settings['data_has_polygon']

    def get_polygon_data(self):
        return self.settings['polygon_data']

    def get_inlets_mode(self):
        return self.settings['inlets_mode']

    def get_inlets_data(self):
        return self.settings['inlets_data']

    def get_inlets_number(self):
        return self.settings['inlets_number']

    def get_inlets_shuffle(self):
        return self.settings['inlets_shuffle']

    def get_outlets_mode(self):
        return self.settings['outlets_mode']

    def get_outlets_data(self):
        return self.settings['outlets_data']

    def get_outlets_number(self):
        return self.settings['outlets_number']

    def get_outlets_shuffle(self):
        return self.settings['outlets_shuffle']

    def get_geological_mode(self):
        return self.settings['geological_mode']

    def get_geological_datafile(self):
        return self.settings['geological_datafile']

    def get_topography_mode(self):
        return self.settings['topography_mode']

    def get_topography_datafile(self):
        return self.settings['topography_datafile']

    def get_orientation_mode(self):
        return self.settings['orientation_mode']

    def get_orientation_datafile(self):
        return self.settings['orientation_datafile']

    def get_faults_mode(self):
        return self.settings['faults_mode']

    def get_faults_datafile(self):
        return self.settings['faults_datafile']

    def get_fractures_mode(self):
        return self.settings['fractures_mode']

    def get_fractures_datafile(self):
        return self.settings['fractures_datafile']

    def get_fractures_densities(self):
        return self.settings['fractures_densities']

    def get_fractures_min_orientation(self):
        return self.settings['fractures_min_orientation']

    def get_fractures_max_orientation(self):
        return self.settings['fractures_max_orientation']
        
    def get_fractures_min_dip(self):
        return self.settings['fractures_min_dip']

    def get_fractures_max_dip(self):
        return self.settings['fractures_max_dip']

    def get_fractures_alpha(self):
        return self.settings['fractures_alpha']

    def get_fractures_min_length(self):
        return self.settings['fractures_min_length']

    def get_fractures_max_length(self):
        return self.settings['fractures_max_length']

    def get_algorithm(self):
        return self.settings['algorithm']

    def get_cost_out(self):
        return self.settings['cost_out']

    def get_cost_aquifer(self):
        return self.settings['cost_aquifer']

    def get_cost_aquiclude(self):
        return self.settings['cost_aquiclude']

    def get_cost_faults(self):
        return self.settings['cost_faults']

    def get_cost_fractures(self):
        return self.settings['cost_fractures']

    def get_cost_conduits(self):
        return self.settings['cost_conduits']

    def get_cost_ratio(self):
        return self.settings['cost_ratio']

    def get_geology_id(self):
        return self.settings['geology_id']
    
    """ ???
    def get_geology_velocity(self):
        return self.settings['geology_velocity']
    """
 
    def get_geology_cost(self):
        return self.settings['geology_cost']

    def get_importance_factor(self):
        return self.settings['importance_factor']

    def get_rand_seed(self):
        return self.settings['rand_seed']

    def get_verbosity(self):
        return self.settings['verbosity']
        
    ##################################
    # get_ methods for uploaded data #
    ##################################

    def get_polygon_vertices(self):
        """
        Get the polygon vertices as a list.
        """
        if self.polygon.polygon is None:
            print("No polygon set.")
            return None
        else:
            return self.polygon.polygon

    def get_inlets(self):
        """
        Get the inlets as an array.
        """
        return self.inlets

    def get_outlets(self):
        """
        Get the outlets as an array.
        """
        return self.outlets

    def get_geology(self):
        """
        Get the geological data as a numpy-array.
        """
        return self.geology.data['geology']['data']

    def get_topography(self):
        """
        Get the topography data as a numpy-array.
        """
        return self.geology.data['topography']['data']

    def get_contact_surface(self):
        """
        Get the lower contact surface of the karst unit as a numpy-array.
        """
        return self.geology.data['contact']['data']

    def get_orientation(self):
        """
        Get the orientation data as two numpy arrays (for x and y components).
        """
        return [self.geology.data['orientationx']['data'], self.geology.data['orientationy']['data']]

    def get_faults(self):
        """
        Get the faults data as a numpy-array.
        """
        return self.geology.data['faults']['data']

    def get_fractures(self, fractures_family=None):
        """
        Get the fractures as a numpy-array.
        
        Parameter
        ---------
        fractures_family : integer
            First family is 0
        """
        if fractures_family is None:
            return self.geology.data['fractures']['data']
        else:
            return self.geology.fractures[fractures_family]['frac_map']

    def get_fractures_numbers(self):
        """
        Get the number of fractures.
        
        Return
        ------
        dictionnary
            A dictionnary with the number of fractures given by user, and with the number of fractures calculated.
        """
        user = self.fractures_numbers
        model = []
        for frac_family in self.geology.fractures:
            model.append(self.geology.fractures[frac_family]['frac_nbr'])
        fracs = {'user' : user, 'model' : model}
        return fracs
    
    def get_mask(self):
        """
        Get the numpy-array mask calculated if a polygon is given.
        """
        if self.mask is None:
            return 'No mask to return.'
        else:
            return self.mask
            
    ################
    # set_ methods #
    ################

    def set_data_has_polygon(self, data_has_polygon):
        """
        Define if study area will be focused on a polygon delimitations.
        
        Parameter
        ---------
        data_has_polygon : bool
            If true, a polygon will be required. 
        """
        self.settings['data_has_polygon'] = data_has_polygon
        return None

    def set_polygon_data(self, polygon_data):
        """
        Define the polygon datafile path.
        Useful only when data_has_polygon is true.
        
        Parameter
        ---------
        polygon_data : string or list
            Polygon datafile path or list of vertices coordinates. 
        """
        self.settings['polygon_data'] = polygon_data
        return None

    def set_inlets_mode(self, inlets_mode):
        """
        Define the inlets mode.
        
        Parameter
        ---------
        inlets_mode : string
            'random'    - Full random points 
            'import'    - Import points
            'composite' - Add n random points to imported points 
        """
        self.settings['inlets_mode'] = inlets_mode
        return None

    def set_inlets_data(self, inlets_data):
        """
        Define the inlets datafile path.
        Useful only when inlets mode is on 'import' or 'composite'.
        
        Parameter
        ---------
        inlets_data : string or list
            Inlets datafile path or list of inlets coordinates. 
        """
        self.settings['inlets_data'] = inlets_data
        return None

    def set_inlets_number(self, inlets_number):
        """
        Define the number of inlets to generate.
        Useful only when inlets mode is on 'random' or 'composite'.
        
        Parameter
        ---------
        inlets_number : string
            Number of inlets to generate. 
        """
        self.settings['inlets_number'] = inlets_number
        return None

    def set_inlets_shuffle(self, inlets_shuffle):
        """
        Define whether to shuffle the order of the inlets randomly.
        Useful only when iterating over inlets.
        
        Parameter
        ---------
        inlets_shuffle : string
            0: don't shuffle, 1: shuffle randomly.
        """
        self.settings['inlets_shuffle'] = inlets_shuffle
        return None

    def set_inlets_per_outlet(self, inlets_per_outlet):
        """
        Define the proportion of inlets to be distributed across each outlet.
        Length of array indicates number of outlets, 
        each integer indicates number of inlets to assign to that outlet, 
        sum of integers = total number of inlets
        Useful only when iterating over inlets and outlets.
        
        Parameter
        ---------
        inlets_per_outlet : [1]: a single iteration with all inlets to one outlet, [1,1,1]: three outlets with one inlet in each,
            [1,2,3]: three outlets with one inlet in the first, 2 inlets in the second, and 3 inlets in the third. 
        """
        self.settings['inlets_per_outlet'] = inlets_per_outlet
        return None
    
    def set_inlets_importance(self, inlets_importance):
        """
        Define the proportion of inlets to be distributed across each iteration.
        Length of array indicates number of inlet iterations, 
        each integer indicates number of inlets to run in that iteration, 
        sum of integers = total number of inlets
        Useful only when iterating over inlets.
        
        Parameter
        ---------
        inlets_importance : [1]: a single iteration with all inlets, [1,1,1]: three iterations with one inlet in each,
            [1,2,3]: three iterations with one inlet in the first, 2 inlets in the second, and 3 inlets in the third. 
        """
        self.settings['inlets_importance'] = inlets_importance
        return None

    def set_outlets_mode(self, outlets_mode):
        """
        Define the outlets mode.
        
        Parameter
        ---------
        outlets_mode : string
            'random'    - Full random points 
            'import'    - Import points
            'composite' - Add n random points to imported points 
        """
        self.settings['outlets_mode'] = outlets_mode
        return None

    def set_outlets_data(self, outlets_data):
        """
        Define the outlets datafile path.
        Useful only when outlets mode is on 'import' or 'composite'.
        
        Parameter
        ---------
        outlets_data : string or list
            Outlets datafile path or list of outlets coordinates. 
        """
        self.settings['outlets_data'] = outlets_data
        return None

    def set_outlets_number(self, outlets_number):
        """
        Define the number of outlets to generate.
        Useful only when outlets mode is on 'random' or 'composite'.
        
        Parameter
        ---------
        outlets_number : string
            Number of outlets to generate. 
        """
        self.settings['outlets_number'] = outlets_number
        return None

    def set_outlets_shuffle(self, outlets_shuffle):
        """
        Define whether to shuffle the order of the outlets randomly.
        Useful only when iterating over outlets.
        
        Parameter
        ---------
        outlets_shuffle : string
            0: don't shuffle, 1: shuffle randomly. 
        """
        self.settings['outlets_shuffle'] = outlets_shuffle
        return None

    def set_outlets_importance(self, outlets_importance):
        """
        Define the proportion of outlets to be distributed across each iteration.
        Length of array indicates number of outlet iterations, 
        each integer indicates number of outlets to run in that iteration, 
        sum of integers = total number of outlets
        Useful only when iterating over outlets.
        
        Parameter
        ---------
        outlets_importance : [1]: a single iteration with all outlets, [1,1,1]: three iterations with one outlet in each,
            [1,2,3]: three iterations with one outlet in the first, 2 outlets in the second, and 3 outlets in the third. 
        """
        self.settings['outlets_importance'] = outlets_importance
        return None
    
    def set_geological_mode(self, geological_mode):
        """
        Define the geological mode.
        
        Parameter
        ---------
        geological_mode : string
            'null'   - No geology
            'import' - Import from old gslib file format
            'gslib'  - Import from new gslib format
            'csv'    - Import from csv
            'image'  - Import via image
        """
        self.settings['geological_mode'] = geological_mode
        return None

    def set_geological_datafile(self, geological_datafile):
        """
        Define the geological datafile path.
        Useful only when geological mode is not 'null'.
        
        Parameter
        ---------
        geological_datafile : string
            Geological datafile path. 
        """
        self.settings['geological_datafile'] = geological_datafile
        return None

    def set_topography_mode(self, topography_mode):
        """
        Define the topography mode.
        
        Parameter
        ---------
        topography_mode : string
            'null'   - No topography
            'import' - Import from old gslib file format
            'gslib'  - Import from new gslib format
            'csv'    - Import from csv
            'image'  - Import via image
        """
        self.settings['topography_mode'] = topography_mode
        return None

    def set_topography_datafile(self, topography_datafile):
        """
        Define the topography datafile path.
        Useful only when topography mode is not 'null'.
        
        Parameter
        ---------
        topography_datafile : string
            topography datafile path. 
        """
        self.settings['topography_datafile'] = topography_datafile
        return None

    def set_orientation_mode(self, orientation_mode):
        """
        Define the orientation mode.
        
        Parameter
        ---------
        orientation_mode : string
            'null'    - No orientation
            'topo'    - Calculate from topography
            'surface' - Calculate from csv file of a surface (useful if using lower surface of karst unit)
        """
        self.settings['orientation_mode'] = orientation_mode
        return None

    def set_orientation_datafile(self, orientation_datafile):
        """
        Define the orientation datafile path.
        Usefull only when orientation mode is not 'null'.
        
        Parameter
        ---------
        orientation_datafile : string 
            orientation datafile path. 
        """
        self.settings['orientation_datafile'] = orientation_datafile
        return None

    def set_faults_mode(self, faults_mode):
        """
        Define the mode for the faults.
        
        Parameter
        ---------
        faults_mode : string
            'null'   - No faults
            'import' - Import from old gslib file format
            'gslib'  - Import from new gslib format
            'csv'    - Import from csv
            'image'  - Import via image
        """
        self.settings['faults_mode'] = faults_mode
        return None

    def set_faults_datafile(self, faults_datafile):
        """
        Define the faults datafile path.
        Useful only when the mode for faults is on 'import' or 'image'.
        
        Parameter
        ---------
        faults_datafile : string
            Faults datafile path. 
        """
        self.settings['faults_datafile'] = faults_datafile
        return None

    def set_fractures_mode(self, fractures_mode):
        """
        Define the mode for the fractures.
        
        Parameter
        ---------
        fractures_mode : string
            'null'   - No fractures
            'import' - Import from old gslib file format
            'gslib'  - Import from new gslib format
            'csv'    - Import from csv
            'image'  - Import via image
            'random' - Generate random fractures
        """
        self.settings['fractures_mode'] = fractures_mode
        return None

    def set_fractures_datafile(self, fracture_datafile):
        """
        Define the fractures datafile path.
        Useful only when the mode for fractures is on 'import' or 'image'.
        
        Parameter
        ---------
        fracture_datafile : string
            Fractures datafile path. 
        """
        self.fracture_datafile = fracture_datafile
        return None

    def set_fractures_densities(self, fractures_densities):
        """
        Define the fractures densitiy for each fractures family.
        Useful only when the mode for fractures is on 'random'.
        
        Parameter
        ---------
        fractures_densities : list
            One density for each family. 
        """
        self.settings['fractures_densities'] = fractures_densities
        return None

    def set_fractures_min_orientation(self, fractures_min_orientation):
        """
        Define the minimum orientation of the fracturation for each fractures family.
        Useful only when the mode for fractures is on 'random'.
        
        Parameter
        ---------
        fractures_min_orientation : list
            One orientation for each family. 
        """
        self.settings['fractures_min_orientation'] = fractures_min_orientation
        return None

    def set_fractures_max_orientation(self, fractures_max_orientation):
        """
        Define the maximum orientation of the fracturation for each fractures family.
        Useful only when the mode for fractures is on 'random'.
        
        Parameter
        ---------
        fractures_max_orientation : list
            One orientation for each family. 
        """
        self.settings['fractures_max_orientation'] = fractures_max_orientation
        return None
        
    def set_fractures_min_dip(self, fractures_min_dip):
        """
        Define the minimum dip of the fracturation for each fractures family.
        Useful only when the mode for fractures is on 'random'.
        
        Parameter
        ---------
        fractures_min_dip : list
            One dip for each family. 
        """
        self.settings['fractures_min_dip'] = fractures_min_dip
        return None

    def set_fractures_max_dip(self, fractures_max_dip):
        """
        Define the maximum dip of the fracturation for each fractures family.
        Useful only when the mode for fractures is on 'random'.
        
        Parameter
        ---------
        fractures_max_dip : list
            One dip for each family. 
        """
        self.settings['fractures_max_dip'] = fractures_max_dip
        return None

    def set_fractures_alpha(self, fractures_alpha):
        """
        Define alpha, a parameter in the fracturation law.
        Useful only when the mode for fractures is on 'random'.
        
        Parameter
        ---------
        fractures_alpha : list
            One alpha for each family.
        """
        self.settings['fractures_alpha'] = fractures_alpha
        return None

    def set_fractures_min_length(self, fractures_min_length):
        """
        Define the minimum lenght for all the fractures.
        Useful only when the mode for fractures is on 'random'.
        
        Parameter
        ---------
        fractures_min_length : float 
        """
        self.settings['fractures_min_length'] = fractures_min_length
        return None

    def set_fractures_max_length(self, fractures_max_length):
        """
        Define the maximum lenght for all the fractures.
        Useful only when the mode for fractures is on 'random'.
        
        Parameter
        ---------
        fractures_max_length : float 
        """
        self.settings['fractures_max_length'] = fractures_max_length
        return None

    def set_algorithm(self, algorithm):
        """
        Define the algorithm to use when calculating travel time to spring.
        
        Parameter
        ---------
        algorithm : string
            'Isotropic2' - 
            'Isotropic3' - 
            'Riemann2'   - 
            'Riemann3'   - 
            See AGD-HFM documentation for full list of options.
        """
        self.settings['algorithm'] = algorithm
        return None
    
    def set_cost_out(self, cost_out):
        """
        Define the fast-marching value for the outside of the study area.
        The value must be between 0 and 1 and should be high to avoid unrealistic conduits.
        
        Parameter
        ---------
        cost_out : float  (default: 0.999)
        """
        self.settings['cost_out'] = cost_out
        return None

    def set_cost_aquifer(self, cost_aquifer):
        """
        Define the fast-marching value for the aquifer cells.
        Should be between 0 and 1 and lower than aquiclude but higher than conduits.
        
        Parameter
        ---------
        cost_aquifer : float  (default: 0.3)
        """
        self.settings['cost_aquifer'] = cost_aquifer
        return None

    def set_cost_aquiclude(self, cost_aquiclude):
        """
        Define the fast-marching value for the aquiclude cells.
        Should be between 0 and 1 and higher than aquiclude but lower than cost_out
        
        Parameter
        ---------
        cost_aquiclude : float  (default: 0.8)
        """
        self.settings['cost_aquiclude'] = cost_aquiclude
        return None

    def set_cost_faults(self, cost_faults):
        """
        Define the fast-marching value for the faults cells.
        Should be between 0 and 1 and between conduits and cost_out. Higher = conduit will avoid faults, lower = conduits will follow faults
        
        Parameter
        ---------
        cost_faults : float  (default: 0.2)
        """
        self.settings['cost_faults'] = cost_faults
        return None

    def set_cost_fractures(self, cost_fractures):
        """
        Define the fast-marching value for the fractures cells.
        Should be between 0 and 1 and between conduits and cost_out. Higher = conduit will avoid fractures, lower = conduits will follow fractures
        
        Parameter
        ---------
        cost_fractures : float (default: 0.2)
        """
        self.settings['cost_fractures'] = cost_fractures
        return None

    def set_cost_conduits(self, cost_conduits):
        """
        Define the fast-marching value for the conduits cells.
        Should be between 0 and 1 but lower than aquifer (for conduits to preferentially follow each other)
        
        Parameter
        ---------
        cost_conduits : float  (default: 0.01)
        """
        self.settings['cost_conduits'] = cost_conduits
        return None

    def set_cost_ratio(self, cost_ratio):
        """
        Define the fast-marching ratio of travel cost parallel to gradient / travel cost prependicular to gradient.
        Should be between 0 and 0.5. 
        
        Parameter
        ---------
        cost_ratio : float  (default: 0.25)
        """
        self.settings['cost_ratio'] = cost_ratio
        return None

    def set_geology_id(self, geology_id):
        """
        Define the geology id (from geology datafile) to consider in the simulation.
        Only needed in 'import' mode.
        
        Parameter
        ---------
        geology_id : list
        """
        self.settings['geology_id'] = geology_id
        return None
    
    def set_geology_cost(self, geology_cost):
        """
        Define the travel cost to apply to each geology id.
        
        Parameter
        ---------
        geology_cost : list
        """
        self.settings['geology_cost'] = geology_cost
        return None

    def set_rand_seed(self, rand_seed):
        """
        Define the random seed.
        May help for reproduicity.
        
        Parameter
        ---------
        rand_seed : integer
        """
        self.settings['rand_seed'] = rand_seed
        if self.settings['rand_seed'] == 0:
            np.random.seed()
        else:
            np.random.seed(self.settings['rand_seed'])
        return None

    def increment_rand_seed(self):
        """
        Increment by one the value of the random seed.
        """
        self.settings['rand_seed'] += 1
        np.random.seed(self.settings['rand_seed'])
        return None

    def set_verbosity(self, verbosity):
        """
        Define the verbosity (how much output to print during runs).
        
        Parameter
        ---------
        verbosity: integer 
            0: print minimal output,  1: print medium output,  2: print max output
        """
        self.settings['verbosity'] = verbosity
        return None
        
    ###################
    # update_ methods #
    ###################
    
    def update_polygon(self):
        """
        Update the polygon settings.
        """
        if self.settings['data_has_polygon']:
            self.polygon.set_polygon(self.settings['polygon_data'])
            self.mask = self.polygon.mask
            self.polygon.inspect_polygon()
        else:
            self.polygon.polygon = None
            self.polygon.mask = None
            self.mask = self.polygon.mask
        return None

    def update_inlets(self):
        """
        Update the inlets settings.
        """
        if self.settings['inlets_mode'] == 'random':
            self.points.generate_points('inlets', self.settings['inlets_number'])
        elif self.settings['inlets_mode'] == 'import':
            self.points.set_points('inlets', self.settings['inlets_data'])
        elif self.settings['inlets_mode'] == 'composite':
            self.points.composite_points('inlets', self.settings['inlets_data'], self.settings['inlets_number'])
        else:
            print('/!\\ Error - unrecognized inlets setting : ', self.settings['inlets_mode'])
            sys.exit()
        self.points.inspect_points()
        self.inlets = self.points.points['inlets'][:]
        return None

    def update_outlets(self):
        """
        Update the outlets settings.
        """
        if self.settings['outlets_mode'] == 'random':
            self.points.generate_points('outlets', self.settings['outlets_number'])
        elif self.settings['outlets_mode'] == 'import':
            self.points.set_points('outlets', self.settings['outlets_data'])
        elif self.settings['outlets_mode'] == 'composite':
            self.points.composite_points('outlets', self.settings['outlets_data'], self.settings['outlets_number'])
        else:
            print('/!\\ Error - unrecognized outlets setting : ', self.settings['outlets_mode'])
            sys.exit()
        self.points.inspect_points()
        self.outlets = self.points.points['outlets'][:]
        return None

    def update_geology(self):
        """
        Update the geology settings.
        """
        if self.settings['geological_mode'] == 'null':
            self.geology.set_data_null('geology')
        elif self.settings['geological_mode'] == 'import':
            self.geology.set_data('geology', self.settings['geological_datafile'])
        elif self.settings['geological_mode'] == 'image':
            self.geology.set_data_from_image('geology', self.settings['geological_datafile'])
        elif self.settings['geological_mode'] == 'gslib':
            self.geology.set_data_from_gslib('geology', self.settings['geological_datafile'])
        elif self.settings['geological_mode'] == 'csv':
            self.geology.set_data_from_csv('geology', self.settings['geological_datafile'])
        else:
            print('/!\\ Error - unrecognized geological mode : ', self.settings['geological_mode'])
            sys.exit() 
        self.geology.compute_stats_on_data()
        if self.settings['polygon_data']:
            self.geology_masked['geology'] = ma.MaskedArray(self.geology.data['geology']['data'], mask=self.mask)
        return None

    def update_topography(self):
        """
        Update the topography settings.
        """
        if self.settings['topography_mode'] == 'null':
            self.geology.set_data_null('topography')
        elif self.settings['topography_mode'] == 'import':
            self.geology.set_data('topography', self.settings['topography_datafile'])
        elif self.settings['topography_mode'] == 'image':
            self.geology.set_data_from_image('topography', self.settings['topography_datafile'])
        elif self.settings['topography_mode'] == 'gslib':
            self.geology.set_data_from_gslib('topography', self.settings['topography_datafile'])
        elif self.settings['topography_mode'] == 'csv':
            self.geology.set_data_from_csv('topography', self.settings['topography_datafile'])
        else:
            print('/!\\ Error - unrecognized topography mode : ', self.settings['topography_mode'])
            sys.exit() 
        self.geology.compute_stats_on_data()
        if self.settings['polygon_data']:
            self.geology_masked['topography'] = ma.MaskedArray(self.geology.data['topography']['data'], mask=self.mask)
        return None

    def update_orientation(self):
        """
        Update the orientation settings.
        """
        if self.settings['orientation_mode'] == 'null':
            self.geology.set_data_null('orientationx')
            self.geology.set_data_null('orientationy')
        elif self.settings['orientation_mode'] == 'topo':
            self.geology.generate_orientations(self.geology.data['topography']['data'])
        elif self.settings['orientation_mode'] == 'surface':
            self.geology.set_data_from_csv('surface', self.settings['orientation_datafile'])
            self.geology.generate_orientations(self.geology.data['surface']['data'])
        else:
            print('/!\\ Error - unrecognized orientation mode : ', self.settings['orientation_mode'])
            sys.exit()
        self.geology.compute_stats_on_data()
        if self.settings['polygon_data']:
            self.geology_masked['orientationx'] = ma.MaskedArray(self.geology.data['orientationx']['data'], mask=self.mask)
            self.geology_masked['orientationy'] = ma.MaskedArray(self.geology.data['orientationy']['data'], mask=self.mask)
        return None

    def update_faults(self):
        """
        Update the faults settings.
        """
        if self.settings['faults_mode'] == 'null':
            self.geology.set_data_null('faults')
        elif self.settings['faults_mode'] == 'import':
            self.geology.set_data('faults', self.settings['faults_datafile'])
        elif self.settings['faults_mode'] == 'image':
            self.geology.set_data_from_image('faults', self.settings['faults_datafile'])
        else:
            print('/!\\ Error - unrecognized faults mode : ', self.settings['faults_mode'])
            sys.exit()
        self.geology.compute_stats_on_data()
        if self.settings['polygon_data']:
            self.geology_masked['faults'] = ma.MaskedArray(self.geology.data['faults']['data'], mask=self.mask)
        return None

    def update_fractures(self):
        """
        Update the fractures settings.
        """
        if self.settings['fractures_mode'] == 'null':
            self.geology.set_data_null('fractures')
        elif self.settings['fractures_mode'] == 'import':
            self.geology.set_data('fractures', self.settings['fractures_datafile'])
        elif self.settings['fractures_mode'] == 'image':
            self.geology.set_data_from_image('fractures', self.settings['fractures_datafile'])
        elif self.settings['fractures_mode'] == 'random':
            self.geology.generate_fractures(self.settings['fractures_numbers'], self.settings['fractures_min_orientation'], self.settings['fractures_max_orientation'], self.settings['fractures_alpha'], self.settings['fractures_min_length'], self.settings['fractures_max_length'])
        else:                                                                                  
            print('/!\\ Error - unrecognized fractures mode : ', self.settings['fractures_mode'])
            sys.exit()
        self.geology.compute_stats_on_data()
        if self.settings['polygon_data']:
            self.geology_masked['fractures'] = ma.MaskedArray(self.geology.data['fractures']['data'], mask=self.mask)
        return None
    
    def update_all(self):
        """
        Update all the parameters.
        """
        self.update_polygon()
        self.update_inlets()
        self.update_outlets()
        self.update_geology()
        self.update_topography()
        self.update_orientation()
        self.update_faults()
        self.update_fractures()
        return None

    def shuffle_inlets(self):
        """
        Shuffle the inlets order.
        """
        np.random.shuffle(self.inlets)
        return None

    def shuffle_outlets(self):
        """
        Shuffle the outlets order.
        """
        np.random.shuffle(self.outlets)
        return None
        
    ############################
    # initialization functions #
    ############################

    # __init__()
    def _get_reference_statistics(self):
        """
        Get the reference statistics for comparing it with karstnet outputs.
        """
        path = os.path.dirname(os.path.abspath(__file__)) + '/' + 'default_files/statistics.xlsx'
        dataframe = pd.read_excel(path, usecols = "B:H,K")
        self.reference_statistics = {}
        for key in dataframe:
            min_val = np.nanmin(dataframe[key])
            max_val = np.nanmax(dataframe[key])
            self.reference_statistics[key] = (min_val,max_val)
        return None

    # __init__()
    def _load_data(self):
        """
        Initialize all the data.
        """
        # Set the grid
        self.grid    = self._set_grid(self.settings['x0'],self.settings['y0'],self.settings['z0'],self.settings['nx'],self.settings['ny'],self.settings['nz'],self.settings['dx'],self.settings['dy'],self.settings['dz'])
        # Set the polygon and its mask
        self.polygon = self._set_polygon(self.grid)
        self.mask    = self.polygon.mask
        # Set the points
        self.points  = self._set_points(self.grid, self.polygon)
        self.inlets  = self.points.points['inlets'][:]
        self.outlets = self.points.points['outlets'][:]
        if self.settings['inlets_shuffle'] == True:
            np.random.shuffle(self.inlets)
        if self.settings['outlets_shuffle'] == True:
            np.random.shuffle(self.outlets)
        # Set the geologic data
        self.geology = self._set_geology(self.grid)
        self.geology.compute_stats_on_data()
        
        self.geology_masked = {}
        if self.mask is not None:
            for key in self.geology.data:
                self.geology_masked[key] = ma.MaskedArray(self.geology.data[key]['data'], mask=self.mask)
        
        return None

    def _set_grid(self, x0, y0, z0, nx, ny, nz, dx, dy, dz):
        """
        Set the grid object.
        """
        grid = Grid(x0, y0, z0, nx, ny, nz, dx, dy, dz)
        return grid

    def _set_polygon(self, grid):
        """
        Set the polygon object.
        """
        polygon = Polygon(grid)
        if self.settings['data_has_polygon']:
            polygon.set_polygon(self.settings['polygon_data'])
            polygon.inspect_polygon()
        return polygon

    def _set_points(self, grid, polygon):
        """
        Set the point manager object.
        """
        points = PointManager(grid, polygon)
        
        # Inlets
        if self.settings['inlets_mode'] == 'random':
            points.generate_points('inlets', self.settings['inlets_number'])
        elif self.settings['inlets_mode'] == 'import':
            points.set_points('inlets', self.settings['inlets_data'])
        elif self.settings['inlets_mode'] == 'composite':
            points.composite_points('inlets', self.settings['inlets_data'], self.settings['inlets_number'])
        else:
            print('/!\\ Error - unrecognized inlets setting : ', self.settings['inlets_mode'])
            sys.exit()
            
        # Outlets
        if self.settings['outlets_mode'] == 'random':
            points.generate_points('outlets', self.settings['outlets_number'])
        elif self.settings['outlets_mode'] == 'import':
            points.set_points('outlets', self.settings['outlets_data'])
        elif self.settings['outlets_mode'] == 'composite':
            points.composite_points('outlets', self.settings['outlets_data'], self.settings['outlets_number'])
        else:
            print('/!\\ Error - unrecognized outlets setting : ', self.settings['outlets_mode'])
            sys.exit()
        points.inspect_points()
        return points

    def _set_geology(self, grid):
        """
        Set the geology manager object.
        """
        geology = GeologyManager(grid)
        
        # Geology
        if self.settings['geological_mode'] == 'null':
            geology.set_data_null('geology')
        elif self.settings['geological_mode'] == 'import':
            geology.set_data('geology', self.settings['geological_datafile'])
        elif self.settings['geological_mode'] == 'image':
            geology.set_data_from_image('geology', self.settings['geological_datafile'])
        elif self.settings['geological_mode'] == 'csv':
            geology.set_data_from_csv('geology', self.settings['geological_datafile'])
        elif self.settings['geological_mode'] == 'gslib':
            geology.set_data_from_gslib('geology', self.settings['geological_datafile'])
        else:
            print('/!\\ Error - unrecognized geological mode : ', self.settings['geological_mode'])
            sys.exit()
            
        # Topography
        if self.settings['topography_mode'] == 'null':
            geology.set_data_null('topography')
        elif self.settings['topography_mode'] == 'csv':
            geology.set_data_from_csv('topography', self.settings['topography_datafile'])
        else:
            print('/!\\ Error - unrecognized topography mode : ', self.settings['topography_mode'])
            sys.exit()
            
        # Orientation
        if self.settings['orientation_mode'] == 'null':
            geology.set_data_null('orientationx')
            geology.set_data_null('orientationy')
        elif self.settings['orientation_mode'] == 'topo':
            geology.generate_orientations(geology.data['topography']['data'])
        elif self.settings['orientation_mode'] == 'surface':
            geology.set_data_from_csv('surface', self.settings['orientation_datafile'])
            geology.generate_orientations(geology.data['surface']['data'])          
        else:
            print('/!\\ Error - unrecognized orientation mode : ', self.settings['orientation_mode'])
            sys.exit()
            
        # Faults
        if self.settings['faults_mode'] == 'null':
            geology.set_data_null('faults')
        elif self.settings['faults_mode'] == 'import':
            geology.set_data('faults', self.settings['faults_datafile'])
        elif self.settings['faults_mode'] == 'image':
            geology.set_data_from_image('faults', self.settings['faults_datafile'])
        else:
            print('/!\\ Error - unrecognized faults mode : ', self.settings['faults_mode'])
            sys.exit()
            
        # Fractures
        if self.settings['fractures_mode'] == 'null':
            geology.set_data_null('fractures')
        elif self.settings['fractures_mode'] == 'import':
            geology.set_data('fractures', self.settings['fractures_datafile'])
        elif self.settings['fractures_mode'] == 'image':
            geology.set_data_from_image('fractures', self.settings['fractures_datafile'])
        elif self.settings['fractures_mode'] == 'random':
            geology.generate_fractures(self.settings['fractures_densities'], self.settings['fractures_alpha'], self.settings['fractures_min_orientation'], self.settings['fractures_max_orientation'], self.settings['fractures_min_dip'], self.settings['fractures_max_dip'], self.settings['fractures_min_length'], self.settings['fractures_max_length'])
        else:
            print('/!\\ Error - unrecognized fractures mode : ', self.settings['fractures_mode'])
            sys.exit()
        return geology
        
""" WORK IN PROGRESS """