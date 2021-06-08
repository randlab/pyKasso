"""
TODO :
- Line 1561 (support 3D)
- Line 1113 - add multithreading ?

- faire du ménage dans les fonctions de visu
- function : get_fractures_numbers ? marche pas ??
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
from .karstnetwork   import KarstNetwork
from .functions      import get_settings

import os
import sys
import copy
import math

import yaml
import numpy    as np
import numpy.ma as ma
import pandas   as pd

import matplotlib
import matplotlib.pyplot as plt

import karstnet as kn
import agd
from agd import Eikonal
from agd.Metrics import Riemann

import concurrent.futures

class SKS():
    """
    Super-class managing all the data classes and simulating stochastic karst networks.
    """

    def __init__(self, yaml_settings_file = None, rand_seed = None):
        """
        Construct a SKS class according to the specified settings datafile.

        Parameters
        ----------
        yaml_settings_file : str
            YAML settings file location.

        Examples
        --------
        >>> catchment = pk.SKS()
        >>> catchment = pk.SKS('example.yaml')
        """

        # print('CAUTION: You are using the development version of this package.') #uncomment this to tell apart the development version and the main version

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

        self.settings['fractures_numbers'] = [int(self.settings['nx']*self.settings['ny']*self.settings['dx']*self.settings['dy']*float(i)) for i in self.settings['fractures_densities']]

        #Set travel cost through each geologic formation
        geology_cost = []
        for elem in self.settings['geology_cost']:
            geology_cost.append(self.settings[elem])
        self.settings['geology_cost'] = geology_cost

        if self.settings['rand_seed'] > 0:
            np.random.seed(self.settings['rand_seed'])

        if rand_seed != None:
            np.random.seed(rand_seed)

        # Avoid issues with text and numbers
        parameters = ['fractures_densities', 'fractures_min_orientation', 'fractures_max_orientation', 'fractures_min_dip', 'fractures_max_dip', 'fractures_alpha', 'fractures_min_length', 'fractures_max_length']
        for parameter in parameters:
            self.settings[parameter] = [float(elem) for elem in self.settings[parameter]]

        self._get_reference_statistics()
        self._load_data()
        self.karst_simulations = []

    ################
    # get_ methods #
    ################

    def get_x0(self):
        """
        Get the origin x-coordinate used in the model.

        Examples
        --------
        >>> x0 = catchment.get_x0()
        """
        return self.settings['x0']

    def get_y0(self):
        """
        Get the origin y-coordinate used in the model.

        Examples
        --------
        >>> y0 = catchment.get_y0()
        """
        return self.settings['y0']

    def get_z0(self):
        """
        Get the origin z-coordinate used in the model.

        Examples
        --------
        >>> z0 = catchment.get_z0()
        """
        return self.settings['z0']

    def get_nx(self):
        """
        Get the number of nodes in the x direction used in the model.

        Examples
        --------
        >>> nx = catchment.get_nx()
        """
        return self.settings['nx']

    def get_ny(self):
        """
        Get the number of cells in the y direction used in the model.

        Examples
        --------
        >>> ny = catchment.get_ny()
        """
        return self.settings['ny']

    def get_nz(self):
        """
        Get the number of cells in the z direction used in the model.

        Examples
        --------
        >>> nz = catchment.get_nz()
        """
        return self.settings['nz']

    def get_dx(self):
        """
        Get the width of cells in the x direction used in the model.

        Examples
        --------
        >>> dx = catchment.get_dx()
        """
        return self.settings['dx']

    def get_dy(self):
        """
        Get the width of cells in the y direction used in the model.

        Examples
        --------
        >>> dy = catchment.get_dy()
        """
        return self.settings['dy']

    def get_dz(self):
        """
        Get the width of cells in the z direction used in the model.

        Examples
        --------
        >>> dz = catchment.get_dz()
        """
        return self.settings['dz']

    def get_data_has_polygon(self):
        """
        Get the boolean value which indicates if a polygon is used or not.

        Examples
        --------
        >>> has_polygon = catchment.get_data_has_polygon()
        """
        return self.settings['data_has_polygon']

    def get_polygon_data(self):
        """
        Get the vertices of the polygon.

        Examples
        --------
        >>> polygon_vertices = catchment.get_polygon_data()
        """
        return self.settings['polygon_data']

    def get_outlets_mode(self):
        return self.settings['outlets_mode']

    def get_outlets_data(self):
        return self.settings['outlets_data']

    def get_outlets_number(self):
        return self.settings['outlets_number']

    def get_outlets_shuffle(self):
        return self.settings['outlets_shuffle']

    def get_outlets_importance(self):
        return self.settings['outlets_importance']

    def get_inlets_mode(self):
        return self.settings['inlets_mode']

    def get_inlets_data(self):
        return self.settings['inlets_data']

    def get_inlets_number(self):
        return self.settings['inlets_number']

    def get_inlets_shuffle(self):
        return self.settings['inlets_shuffle']

    def get_inlets_per_outlet(self):
        return self.settings['inlets_per_outlet']

    def get_inlets_importance(self):
        return self.settings['inlets_importance']

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
        return self.geology.data['surface']['data']

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

    def set_geological_mode(self, geological_mode):
        """
        Define the geological mode.

        Parameter
        ---------
        geological_mode : string
            'null'   - No geology
            'gslib'  - Import geology via gslib
            'csv'    - Import geology via csv
            'image'  - Import geology via image
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
            'gslib'  - Import topography via gslib
            'csv'    - Import topography via csv
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
            'gslib'  - Import faults via gslib
            'csv'    - Import faults via csv
            'image'  - Import faults via image
        """
        self.settings['faults_mode'] = faults_mode
        return None

    def set_faults_datafile(self, faults_datafile):
        """
        Define the faults datafile path.
        Useful only when the mode for faults is on 'gslib', 'csv' or 'image'.

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
            'gslib'  - Import fractures via gslib
            'csv'    - Import fractures via csv
            'image'  - Import fractures via image
            'random' - Generate fractures randomly
        """
        self.settings['fractures_mode'] = fractures_mode
        return None

    def set_fractures_datafile(self, fracture_datafile):
        """
        Define the fractures datafile path.
        Useful only when the mode for fractures is on 'gslib', 'csv' or 'image'.

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
        Only needed in 'gslib' or 'csv' mode.

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
            print('/!\\ Error : unrecognized inlets setting.')
            sys.exit()
        self.points.inspect_points()
        self.inlets = self.points.points['inlets'][:]
        if self.settings['inlets_shuffle'] == True:
            np.random.shuffle(self.inlets)
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
            print('/!\\ Error : unrecognized outlets setting.')
            sys.exit()
        self.points.inspect_points()
        self.outlets = self.points.points['outlets'][:]
        if self.settings['outlets_shuffle'] == True:
            np.random.shuffle(self.outlets)
        return None

    def update_geology(self):
        """
        Update the geology settings.
        """
        if self.settings['geological_mode'] == 'null':
            self.geology.set_data_null('geology')
        elif self.settings['geological_mode'] == 'gslib':
            self.geology.set_data_from_gslib('geology', self.settings['geological_datafile'])
        elif self.settings['geological_mode'] == 'csv':
            self.geology.set_data_from_csv('geology', self.settings['geological_datafile'])
        elif self.settings['geological_mode'] == 'image':
            self.geology.set_data_from_image('geology', self.settings['geological_datafile'])
        else:
            print('/!\\ Error : unrecognized geological mode', self.settings['geological_mode'])
            sys.exit()
        self.geology.compute_stats_on_data('geology')
        if self.settings['polygon_data']:
            self.geology_masked['geology'] = ma.MaskedArray(self.geology.data['geology']['data'], mask=self.mask)
        return None

    def update_topography(self):
        """
        Update the topography settings.
        """
        if self.settings['topography_mode'] == 'null':
            self.geology.set_data_null('topography')
        elif self.settings['topography_mode'] == 'csv':
            self.geology.set_data_from_csv('topography', self.settings['topography_datafile'])
        else:
            print('/!\\ Error : unrecognized topography mode', self.settings['topography_mode'])
            sys.exit()
        self.geology.compute_stats_on_data('topography')
        if self.settings['polygon_data']:
            self.geology_masked['topography'] = ma.MaskedArray(self.geology.data['topography']['data'], mask=self.mask[:,:,0])
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
            print('/!\\ Error : unrecognized orientation mode', self.settings['orientation_mode'])
            sys.exit()
        self.geology.compute_stats_on_data('orientationx')
        self.geology.compute_stats_on_data('orientationy')
        if self.settings['polygon_data']:
            self.geology_masked['orientationx'] = ma.MaskedArray(self.geology.data['orientationx']['data'], mask=self.mask[:,:,0])
            self.geology_masked['orientationy'] = ma.MaskedArray(self.geology.data['orientationy']['data'], mask=self.mask[:,:,0])
        return None

    def update_faults(self):
        """
        Update the faults settings.
        """
        if self.settings['faults_mode'] == 'null':
            self.geology.set_data_null('faults')
        elif self.settings['faults_mode'] == 'gslib':
            self.geology.set_data_from_gslib('faults', self.settings['faults_datafile'])
        elif self.settings['faults_mode'] == 'csv':
            self.geology.set_data_from_csv('faults', self.settings['faults_datafile'])
        elif self.settings['faults_mode'] == 'image':
            self.geology.set_data_from_image('faults', self.settings['faults_datafile'])
        else:
            print('/!\\ Error : unrecognized faults mode.')
            sys.exit()
        self.geology.compute_stats_on_data('faults')
        if self.settings['polygon_data']:
            self.geology_masked['faults'] = ma.MaskedArray(self.geology.data['faults']['data'], mask=self.mask)
        return None

    def update_fractures(self):
        """
        Update the fractures settings.
        """
        if self.settings['fractures_mode'] == 'null':
            self.geology.set_data_null('fractures')
        elif self.settings['fractures_mode'] == 'gslib':
            self.geology.set_data_from_gslib('fractures', self.settings['fractures_datafile'])
        elif self.settings['fractures_mode'] == 'csv':
            self.geology.set_data_from_csv('fractures', self.settings['fractures_datafile'])
        elif self.settings['fractures_mode'] == 'image':
            self.geology.set_data_from_image('fractures', self.settings['fractures_datafile'])
        elif self.settings['fractures_mode'] == 'random':
            self.geology.generate_fractures(self.settings['fractures_densities'],
                                            self.settings['fractures_alpha'],
                                            self.settings['fractures_min_orientation'],
                                            self.settings['fractures_max_orientation'],
                                            self.settings['fractures_min_dip'],
                                            self.settings['fractures_max_dip'],
                                            self.settings['fractures_min_length'],
                                            self.settings['fractures_max_length'])
            self.geology.rasterize_fracture_network
        else:
            print('/!\\ Error : unrecognized fractures mode.')
            sys.exit()
        self.geology.compute_stats_on_data('fractures')
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
        # The grid
        self.grid = self._set_grid(self.settings['x0'],self.settings['y0'],self.settings['z0'],
                                   self.settings['nx'],self.settings['ny'],self.settings['nz'],
                                   self.settings['dx'],self.settings['dy'],self.settings['dz'])
        # The polygon and its mask
        self.polygon = self._set_polygon(self.grid)
        self.mask    = self.polygon.mask
        # The geologic data
        self.geology = self._set_geology(self.grid)
        self.geology.compute_stats_on_data('geology')
        self.geology.compute_stats_on_data('faults')
        self.geology.compute_stats_on_data('fractures')
        # The points
        self.points  = self._set_points(self.grid, self.polygon, self.geology)
        self.inlets  = self.points.points['inlets'][:]
        self.outlets = self.points.points['outlets'][:]
        if self.settings['inlets_shuffle'] == True:
            np.random.shuffle(self.inlets)
        if self.settings['outlets_shuffle'] == True:
            np.random.shuffle(self.outlets)

        self.geology_masked = {}
        if self.mask is not None:
            mask = self.mask
            for key in self.geology.data:
                if key in ['topography', 'orientationx', 'orientationy']:
                    mask = self.mask[:,:,0]
                self.geology_masked[key] = ma.MaskedArray(self.geology.data[key]['data'], mask=mask)

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

    def _set_points(self, grid, polygon, geology):
        """
        Set the point manager object.
        """
        points = PointManager(grid, polygon, geology)

        ###Inlets###
        if self.settings['inlets_mode'] == 'random':
            points.generate_points('inlets', self.settings['inlets_number'])
        elif self.settings['inlets_mode'] == 'import':
            points.set_points('inlets', self.settings['inlets_data'])
        elif self.settings['inlets_mode'] == 'composite':
            points.composite_points('inlets', self.settings['inlets_data'], self.settings['inlets_number'])
        else:
            print('/!\\ Error : unrecognized inlets setting.')
            sys.exit()

        ###Outlets###
        if self.settings['outlets_mode'] == 'random':
            points.generate_points('outlets', self.settings['outlets_number'])
        elif self.settings['outlets_mode'] == 'import':
            points.set_points('outlets', self.settings['outlets_data'])
        elif self.settings['outlets_mode'] == 'composite':
            points.composite_points('outlets', self.settings['outlets_data'], self.settings['outlets_number'])
        else:
            print('/!\\ Error : unrecognized outlets setting.')
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
        elif self.settings['geological_mode'] == 'gslib':
            geology.set_data_from_gslib('geology', self.settings['geological_datafile'])
        elif self.settings['geological_mode'] == 'csv':
            geology.set_data_from_csv('geology', self.settings['geological_datafile'])
        elif self.settings['geological_mode'] == 'image':
            geology.set_data_from_image('geology', self.settings['geological_datafile'])
        else:
            print('/!\\ Error : unrecognized geological mode', self.settings['geological_mode'])
            sys.exit()

        # Topography
        if self.settings['topography_mode'] == 'null':
            geology.set_data_null('topography')
        elif self.settings['topography_mode'] == 'csv':
            geology.set_data_from_csv('topography', self.settings['topography_datafile'])
        else:
            print('/!\\ Error : unrecognized topography mode', self.settings['topography_mode'])
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
            print('/!\\ Error : unrecognized orientation mode', self.settings['orientation_mode'])
            sys.exit()

        # Faults
        if self.settings['faults_mode'] == 'null':
            geology.set_data_null('faults')
        elif self.settings['faults_mode'] == 'gslib':
            geology.set_data_from_gslib('faults', self.settings['faults_datafile'])
        elif self.settings['faults_mode'] == 'csv':
            geology.set_data_from_csv('faults', self.settings['faults_datafile'])
        elif self.settings['faults_mode'] == 'image':
            geology.set_data_from_image('faults', self.settings['faults_datafile'])
        else:
            print('/!\\ Error : unrecognized faults mode', self.settings['faults_mode'])
            sys.exit()

        # Fractures
        if self.settings['fractures_mode'] == 'null':
            geology.set_data_null('fractures')
        elif self.settings['fractures_mode'] == 'gslib':
            geology.set_data_from_gslib('fractures', self.settings['fractures_datafile'])
        elif self.settings['fractures_mode'] == 'csv':
            geology.set_data_from_csv('fractures', self.settings['fractures_datafile'])
        elif self.settings['fractures_mode'] == 'image':
            geology.set_data_from_image('fractures', self.settings['fractures_datafile'])
        elif self.settings['fractures_mode'] == 'random':
            geology.generate_fractures(self.settings['fractures_densities'],
                                       self.settings['fractures_alpha'],
                                       self.settings['fractures_min_orientation'],
                                       self.settings['fractures_max_orientation'],
                                       self.settings['fractures_min_dip'],
                                       self.settings['fractures_max_dip'],
                                       self.settings['fractures_min_length'],
                                       self.settings['fractures_max_length'])
            geology.rasterize_fracture_network()
        else:
            print('/!\\ Error : unrecognized fractures mode', self.settings['fractures_mode'])
            sys.exit()

        return geology

    ##############################
    # Karst simulation functions #
    ##############################

    def compute_karst_network(self, verbose = False):
        """
        Compute the karst network according to the parameters.

        Save the results in the `karst_simulations` list attribute.
        """

        # 1 - Initialize the parameters
        self._initialize_karst_network_parameters()

        # 2 - Compute conduits for each generation & store nodes and edges for network
        self._compute_iterations_karst_network()

        # 3 - Calculate the karst network statistics indicators with karstnet and save karst network
        karstnet_edges = list(self.edges.values()) #convert to format needed by karstnet (list)
        karstnet_nodes = copy.deepcopy(self.nodes) #convert to format needed by karstnet (dic with only coordinates) - make sure to only modify a copy!
        for key, value in karstnet_nodes.items():  #drop last item in list (the node type) for each dictionary entry
            value.pop()
        try:
            k = kn.KGraph(karstnet_edges, karstnet_nodes)  #make graph - edges must be a list, and nodes must be a dic of format {nodeindex: [x,y]}
            stats = k.characterize_graph(verbose)
        except:
            print("kn.KGraph() failed, $stats and $k set to None")
            k = None
            stats = None

        # 4 - Store all the relevant data for this network in dictionaries:
        maps = copy.deepcopy(self.maps)                  #store copy of maps
        points = {}
        points['inlets']  = copy.deepcopy(self._inlets)  #store copy of inlets in pandas df format with info about iteration and inlet/outlet assignment (this will be used in plotting later)
        points['outlets'] = copy.deepcopy(self._outlets) #store copy of outlets in pandas df format with info about iteration and inlet/outlet assignment (this will be used in plotting later)
        network = {}
        network['edges'] = copy.deepcopy(self.edges) #store copy of edges list
        network['nodes'] = copy.deepcopy(self.nodes) #store copy of nodes list
        network['karstnet'] = copy.deepcopy(k)       #store copy of karstnet network object (including graph)
        config = copy.deepcopy(self.settings)        #store copy of settings for the run being stored
        self.karst_simulations.append(KarstNetwork(maps, points, network, stats, config))

        # 5 - Return inlets and outlets to their original format:
        #Chloe: this is to convert inlets and outlets back from pandas dataframes
        #if the private self._inlets and self._outlets variables are working correctly, this step may be removed
        self.inlets  = np.asarray(self._inlets)[:,0:2]    #return inlets to array format without info about iteration and inlet/outlet assignment (important for later iterations)
        self.outlets = np.asarray(self._outlets)[:,0:2]   #return outlets to array format without info about iteration and inlet/outlet assignment (important for later iterations)
        return None

    # 1
    def _initialize_karst_network_parameters(self):
        """
        Initialize the karst network parameters.
        """
        #Iterations
        self.nbr_iteration  = len(self.settings['outlets_importance'])*len(self.settings['inlets_importance']) #total number of iterations that will occur
        outlets_repartition = self._repartition_points(self.outlets, self.settings['outlets_importance']) #correct for outlets_importance not summing to correct number of actual outlets
        self._outlets       = self._distribute_outlets(outlets_repartition)  # distribute outlets to iterations and store as a semi-private variable for internal use only
        inlets_repartition  = self._repartition_points(self.inlets, self.settings['inlets_per_outlet']) #correct for inlets_per_outlet not summing to correct number of actual inlets
        self._inlets        = self._distribute_inlets(inlets_repartition)  # distribute inlets to outlets and store as a semi-private variable for internal use only

        inlets  = []
        for o,outlet in enumerate(self._outlets):                 #loop over outlets
            inlets_current = self._inlets[self._inlets[:,2]==o]    #select only inlets assigned to current outlet
            repartition    = self._repartition_points(inlets_current, self.settings['inlets_importance']) #correct for inlets_importance not summing to correct number of actual inlets available for current outlet
            i = 0            #total inlet counter for current outlet
            for k,n in enumerate(repartition):   #get iteration index and number of inlets assigned to that iteration
                for j in range(n):               #loop over each inlet in current iteration
                    inlets.append((inlets_current[i,0], inlets_current[i,1], inlets_current[i,2], inlets_current[i,3], inlets_current[i,4], k))
                    i += 1
        self._inlets  = pd.DataFrame(inlets,        columns = ['x','y', 'outlet', 'outletx', 'outlety', 'inlet_iteration']) #store as pandas df for easier handling
        self._outlets = pd.DataFrame(self._outlets, columns = ['x','y', 'outlet_iteration']) #store as df for easier indexing

        # Raster maps
        self.maps            = {}
        self.maps['outlets'] = np.full((self.grid.nx, self.grid.ny), np.nan) #map of null values where each cell with an outlet will have the index of that outlet
        self.maps['nodes']   = np.full((self.grid.nx, self.grid.ny), np.nan) #map of null values where each cell that has a node will be updated with that node index
        self.maps['cost']    = np.zeros((self.nbr_iteration, self.grid.nx, self.grid.ny)) #cost of travel through each cell
        self.maps['alpha']   = np.zeros((self.nbr_iteration, self.grid.nx, self.grid.ny)) #cost of travel along gradient through each cell
        self.maps['beta']    = np.zeros((self.nbr_iteration, self.grid.nx, self.grid.ny)) #cost of travel perpendicular to gradient through each cell
        self.maps['time']    = np.zeros((self.nbr_iteration, self.grid.nx, self.grid.ny)) #travel time to outlet from each cell
        self.maps['karst']   = np.zeros((self.nbr_iteration, self.grid.nx, self.grid.ny)) #presence/absence of karst conduit in each cell

        # Vector maps
        self.nodes     = {} #empty dic to store nodes (key: nodeID, val: [x, y, type])
        self.edges     = {} #empty dic to store edges (key: edgeID, val: [inNode, outNode])
        self.n         = 0  #start node counter at zero
        self.e         = 0  #start edge counter at zero
        self.geodesics = [] #empty list to store raw fast-marching path output

        # Set up fast-marching:
        # Note: AGD-HFM library has different indexing, so model dimensions must be [ny,nx],
        # and model extent must be [ymin,ymax, xmin,xmax] (NOT x0,y0)
        self.riemannMetric = []                    #this changes at every iteration, but cannot be stored?
        self.fastMarching = agd.Eikonal.dictIn({
            'model'             : self.settings['algorithm'],      #set algorithm from settings file ('Isotropic2', 'Isotropic3', 'Riemann2', 'Riemann3')
            'order'             : 2,               #recommended setting: 2 (replace by variable)
            'exportValues'      : 1,               #export the travel time field
            'exportGeodesicFlow': 1                #export the walker paths (i.e. the conduits)
        })
        self.fastMarching.SetRect(                 #give the fast-marching algorithm the model grid
            sides=[[self.grid.ymin, self.grid.ymax],    #leftmost edge, rightmost edge (NOT centerpoint)
                   [self.grid.xmin, self.grid.xmax]],   #bottom edge,   top edge (NOT centerpoint)
            dims=[self.grid.ny, self.grid.nx])      #number of cells, number of cells
        return None

    # 1
    def _repartition_points(self, points, importance):
        '''Correct for integers in importance factors list not summing correctly to total number of points'''
        repartition = []
        proportion  = []
        total       = float(sum(importance))     #total number of points as assigned (this may not be the actual total)
        for i,n in enumerate(importance):        #get iteration index and number of points being assigned to that iteration
            proportion.append(float(n)/total)    #percent of points to use this iteration
            repartition.append(round(proportion[i]*len(points))) #number of points to use this iteration (need to round bc percentage will not result in whole number)
        repartition[-1] = len(points)-sum(repartition[0:-1])     #leftover points are assignd to last iteration
        return repartition

    #1
    def _distribute_outlets(self, repartition):
        '''Distribute points into iteration groups based on the corrected repartition from the importance settings'''
        outlets = []
        i = 0              #point counter
        for k,n in enumerate(repartition):        #get iteration index and number of points being assigned to that iteration
            for j in range(n):                    #loop over each point in current iteration
                outlets.append((self.outlets[i][0], self.outlets[i][1], k))   #append point with its iteration number
                i += 1                                                                #increment point index by 1
        return np.asarray(outlets)

    #1
    def _distribute_inlets(self, repartition):
        '''Distribute points into iteration groups based on the corrected repartition from the importance settings'''
        inlets = []
        i = 0              #point counter
        for k,n in enumerate(repartition):        #get iteration index and number of points being assigned to that iteration
            for j in range(n):                    #loop over each point in current iteration
                inlets.append((self.inlets[i][0], self.inlets[i][1], k, self.outlets[k,0], self.outlets[k,1]))   #append point with its iteration number
                i += 1                                                                #increment point index by 1
        return np.asarray(inlets)

    # 2
    def _compute_iterations_karst_network(self):
        """
        Compute each generation of karst conduits.
        """
        if self.settings['verbosity'] > 0:
            print('-START-')

        # Define outlets map according to outlets emplacements
        self._compute_outlets_map()  #assign outlet indices

        #Plot for debugging:
        if self.settings['verbosity'] > 1:
            f,ax = plt.subplots(1,1)
            ax.imshow(self.geology.data['geology']['data'], extent=self.grid.extent, origin='lower', cmap='gray_r', alpha=0.5) #for debugging

        ## Set up iteration structure:
        iteration = 0                                   #initialize total iteration counter
        for outlet_iteration in range(len(self.settings['outlets_importance'])):  #loop over outlet iterations
            if self.settings['verbosity'] > 2:
                print('Total Iteration', iteration, 'Outlet iteration:', outlet_iteration)
            outlets_current = self._outlets[self._outlets.outlet_iteration==outlet_iteration]  #get the outlets assigned to the current outlet iteration
            if self.settings['verbosity'] > 1:
                ax.scatter(outlets_current.x,outlets_current.y, c='c')         #debugging
            for o,outlet in outlets_current.iterrows():                     #loop over outlets in current outlet iteration
                if self.settings['verbosity'] > 2:
                    print('\t Current outlet index:', outlet.name)             #debugging
                if self.settings['verbosity'] > 1:
                    ax.annotate(str(outlet_iteration), (outlet.x,outlet.y))
                inlets_outlet = self._inlets[self._inlets.outlet==outlet.name]          #get the inlets assigned to current outlet
                if self.settings['verbosity'] > 1:
                    print('\t Inlets assigned to this outlet:\n', inlets_outlet)
                for inlet_iteration in range(len(self.settings['inlets_importance'])): #loop over inlet iterations
                    if self.settings['verbosity'] > 2:
                        print('\t\t Inlet iteration:', inlet_iteration)
                    inlets_current = inlets_outlet[inlets_outlet.inlet_iteration==inlet_iteration] #get the inlets assigned to the current inlet iteration
                    if self.settings['verbosity'] > 2:
                        print(inlets_current)                                                         #debugging
                    if self.settings['verbosity'] > 1:
                        ax.scatter(inlets_current.x, inlets_current.y)
                    for i,inlet in inlets_current.iterrows():                                 #loop over inlet in current inlet iteration
                        if self.settings['verbosity'] > 1:
                            ax.annotate(str(outlet_iteration)+'-'+str(inlet_iteration), (inlet.x,inlet.y))  #debugging
                        self._outlets.loc[self._outlets.index==outlet.name, 'iteration'] = iteration   #store total iteration counter
                        self._inlets.loc[ self._inlets.index ==inlet.name,  'iteration'] = iteration   #store total iteration counter

                    #Compute travel time maps and conduit network:
                    if self.settings['algorithm'] == 'Isotropic2':
                        self._compute_cost_map(iteration)
                        self._compute_time_map_isotropic(iteration)
                        self._compute_karst_map(iteration)
                    elif self.settings['algorithm'] == 'Riemann2':
                        self._compute_cost_map(iteration)
                        self._compute_alpha_map(iteration)
                        self._compute_beta_map(iteration)
                        self._compute_riemann_metric(iteration)
                        self._compute_time_map_riemann(iteration)
                        self._compute_karst_map(iteration)
                    else:
                        print('Unrecognized algorithm', self.settings['algorithm'])
                    iteration = iteration + 1                                              #increment total iteration number by 1

                    if self.settings['verbosity'] > 0:
                        print('iteration:{}/{}'.format(iteration+1,self.nbr_iteration))
        if self.settings['verbosity'] > 0:
            print('- END -')
        return None

    # 2
    def _compute_outlets_map(self):
        """
        Compute the outlets map (array indicating location of outlets as their index and everywhere else as nan).
        """
        for i,outlet in self._outlets.iterrows():
            X = self.grid.get_i(outlet.x)
            Y = self.grid.get_j(outlet.y)
            self.maps['outlets'][X][Y] = i
        return None

    # 2
    def _compute_cost_map(self, iteration):
        """
        Compute the cost map (how difficult it is to traverse each cell).
        """
        # If it's the first iteration, iniatialize the cost map according to the geological settings.
        if iteration == 0:
            # Geology
            if self.geology.data['geology']['mode'] == 'null':
                self.maps['cost'][0] = np.full((self.grid.nx, self.grid.ny), self.settings['cost_aquifer']) #every cell has the same travel cost and is part of the aquifer
            elif self.geology.data['geology']['mode'] == 'image':
                self.maps['cost'][0] = np.where(self.geology.data['geology']['data'][:,:,0]==1, self.settings['cost_aquiclude'], self.settings['cost_aquifer'])
            elif self.geology.data['geology']['mode'] == 'csv' or self.geology.data['geology']['mode'] == 'gslib':

                tableFMM = {}
                if len(self.settings['geology_id']) != len(self.settings['geology_cost']):
                    print("- _compute_cost_map() - Error : number of lithologies does not match with number of FMM code.")
                    sys.exit()

                for geology, codeFMM in zip(self.settings['geology_id'], self.settings['geology_cost']):
                    if geology in self.geology.data['geology']['stats']['ID']:
                        tableFMM[geology] = codeFMM
                    else:
                        print("- initialize_costMap() - Warning : no geology n {} found.".format(geology))
                        tableFMM[geology] = self.settings['cost_out']

                for z in range(self.grid.nz):
                    for y in range(self.grid.ny):
                        for x in range(self.grid.nx):
                            geology = self.geology.data['geology']['data'][x][y][z]
                            self.maps['cost'][0][x][y] = tableFMM[geology]
            else:
                print('geology mode', self.geology.data['geology']['mode'], 'not supported')
                sys.exit()

            # Faults
            self.maps['cost'][0] = np.where(self.geology.data['faults']['data'][:,:,0] > 0, self.settings['cost_faults'], self.maps['cost'][0])

            # Fractures
            self.maps['cost'][0] = np.where(self.geology.data['fractures']['data'][:,:,0] > 0, self.settings['cost_fractures'], self.maps['cost'][0])

            # If out of polygon
            if self.mask is not None:
                self.maps['cost'][0] = np.where(self.mask[:,:,0]==1, self.settings['cost_out'], self.maps['cost'][0])

        else: #if not the first iteration
            self.maps['cost'][iteration] = self.maps['cost'][iteration-1] #set cost map to be same as previous iteration
            self.maps['cost'][iteration] = np.where(self.maps['karst'][iteration-1] > 0, self.settings['cost_conduits'], self.maps['cost'][iteration]) #where karst conduits are present from previous iteration, set cost to conduit cost, elsewhere, leave unchanged
        return None

    # 2
    def _compute_alpha_map(self, iteration):
        """
        Compute the alpha map: travel cost in the same direction as the gradient.
        Cost map * topography map, so that the cost is higher at higher elevations, encouraging conduits to go downgradient.
        """
        if self.settings['topography_mode'] != 'null':
            self.maps['alpha'][iteration] = self.maps['cost'][iteration] * self.geology.data['topography']['data']
        elif self.settings['orientation_mode'] == 'surface':
            self.maps['alpha'][iteration] = self.maps['cost'][iteration] * self.geology.data['surface']['data']
        else:
            self.maps['alpha'][iteration] = self.maps['cost'][iteration]
        return None

    # 2
    def _compute_beta_map(self, iteration):
        """
        Compute the beta map: travel cost perpendicular to the gradient.
        If beta is higher than alpha, conduits will follow the steepest gradient.
        If beta is lower than alpha, conduits will follow contours.
        """
        self.maps['beta'][iteration] = self.maps['alpha'][iteration] / self.settings['cost_ratio']
        return None

    # 2
    def _compute_riemann_metric(self, iteration):
        """
        Compute the riemann metric: Define the Riemannian metric needed as input for the anisotropic fast marching.
        """
        self.riemannMetric = agd.Metrics.Riemann.needle([self.geology.data['orientationx']['data'],
                                                        self.geology.data['orientationy']['data']],
                                                        self.maps['alpha'][iteration],
                                                        self.maps['beta'][iteration])
        return None

    # 2
    def _compute_time_map_isotropic(self, iteration):
        """
        Compute the travel time map (how long it takes to get to the outlet from each cell),
        using the isotropic agd-hfm fast-marching algorithm, and store travel time map.
        Note: the AGD-HFM library uses different indexing, so x and y indices are reversed for inlets and outlets.
        """
        self.fastMarching['seeds'] = np.rot90([self._outlets[self._outlets.iteration==iteration].x, self._outlets[self._outlets.iteration==iteration].y], k=3)         #set the outlets for this iteration
        self.fastMarching['tips']  = np.rot90([self._inlets[self._inlets.iteration==iteration].x, self._inlets[self._inlets.iteration==iteration].y], k=3) #select inlets for current iteration
        self.fastMarching['cost']  = np.transpose(self.maps['cost'][iteration])
        #self.fastMarching['seeds']  = np.rot90([self._outlets[self._outlets.iteration==iteration].x, self._outlets[self._outlets.iteration==iteration].y], k=3)         #set the outlets for this iteration
        #self.fastMarching['tips']   = np.rot90([self._inlets[self._inlets.iteration==iteration].x, self._inlets[self._inlets.iteration==iteration].y], k=3) #select inlets for current iteration
        #self.fastMarching['cost']   = self.maps['cost'][iteration]                      #set the isotropic travel cost through each cell
        self.fastMarching['verbosity'] = self.settings['verbosity']           #set verbosity of hfm run
        self.fastMarchingOutput      = self.fastMarching.Run()                 #run the fast marching algorithm and store the outputs
        self.maps['time'][iteration] = np.transpose(self.fastMarchingOutput['values'])  #store travel time maps
        return None

    # 2
    def _compute_time_map_riemann(self, iteration):
        """
        Compute the travel time map (how long it takes to get to the outlet from each cell),
        using the anisotropic agd-hfm fast-marching algorithm, and store travel time map.
        Note: the AGD-HFM library uses different indexing, so x and y indices are reversed for inlets and outlets.
        """
        self.fastMarching['seeds']     = np.rot90([self._outlets[self._outlets.iteration==iteration].x, self._outlets[self._outlets.iteration==iteration].y], k=3)   #set the outlets for this iteration
        self.fastMarching['tips']      = np.rot90([self._inlets[self._inlets.iteration==iteration].x, self._inlets[self._inlets.iteration==iteration].y], k=3)       #select inlets for current iteration
        self.fastMarching['metric']    = self.riemannMetric                #set the travel cost through each cell
        self.fastMarching['verbosity'] = self.settings['verbosity']        #set verbosity of hfm run
        self.fastMarchingOutput        = self.fastMarching.Run()           #run the fast marching algorithm and store the outputs
        self.maps['time'][iteration]   = self.fastMarchingOutput['values'] #store travel time maps
        self.geodesics.append(self.fastMarchingOutput['geodesics'])        #store fastest travel paths
        return None

    # 2
    def _compute_karst_map(self, iteration):
        """
        Compute the karst map based on the paths from agd-hfm.
        Array of all zeros, with ones in cells containing a karst conduit.
        """
        if iteration > 0:           # Get karst map from previous iteration (except for the very first iteration)
            self.maps['karst'][iteration] = self.maps['karst'][iteration-1]

        #Debugging:
        #Chloe: this should stay in since it is very useful if there are problems
        #f1,ax1 = plt.subplots(1,1, figsize=(10,10))  #debugging
        #ax1.imshow(self.maps['karst'][iteration], origin='lower', extent=self.grid.extent, cmap='gray_r')
        #ax1.imshow(self.maps['nodes'], origin='lower', extent=self.grid.extent, cmap='gray_r')
        #ax1.scatter(self._outlets[self._outlets.iteration==iteration].x, self._outlets[self._outlets.iteration==iteration].y, c='c', s=100)
        #ax1.scatter(self._inlets[self._inlets.iteration==iteration].x, self._inlets[self._inlets.iteration==iteration].y, c='orange', s=100)

        #Loop over conduit paths generated by fast marching:
        for path in self.fastMarchingOutput['geodesics']:   #loop over conduit paths in this iteration (there is one path from each inlet)
            merge = False                                   #reset indicator for whether this conduit has merged with an existing conduit
            for p in range(path.shape[1]):                  #loop over points making up this conduit path
                point = path[:,p]                           #get coordinates of current point
                [[iy,ix],error]  = self.fastMarching.IndexFromPoint(point) #convert to coordinates to indices, /!\ returning iy first then ix
                #ax1.scatter(point[1],point[0], c='g',s=5)  #debugging

                #Place nodes and links:
                if np.isnan(self.maps['nodes'][ix,iy]):                                    #if there is no existing conduit node here
                    if ~np.isnan(self.maps['outlets'][ix,iy]):                              #if there is an outlet here (cell value is not nan)
                        outlet = self._outlets.iloc[int(self.maps['outlets'][ix,iy])]         #get the outlet coordinates using the ID in the outlets map
                        self.nodes[self.n]             = [outlet.y, outlet.x, 'outfall']     #add a node at the outlet coordinates (with the node type for SWMM)
                        self.maps['nodes'][ix,iy] = self.n                                   #update node map with node index
                        #ax1.scatter(outlet.x,outlet.y, marker='o', c='b')                   #debugging
                        if p > 0:                                                           #if this is not the first point (i.e. the inlet) in the current path
                            if merge == False:                                               #if this conduit has not merged with an existing conduit
                                self.edges[self.e] = [self.n-1, self.n]                       #add an edge connecting the previous node to the current node
                                self.e = self.e+1                                             #increment edge counter up by one
                                #ax1.plot((self.nodes[self.n][0], self.nodes[self.n-1][0]),(self.nodes[self.n][1], self.nodes[self.n-1][1]))
                            else:                                                          #if this conduit HAS merged with an existing conduit
                                [[fromix,fromiy],error]  = self.fastMarching.IndexFromPoint(path[:,p-1]) #get xy indices of previous point in current conduit path
                                n_from = self.maps['nodes'][fromix,fromiy]                    #get node index of the node already in the cell where the previous point was
                                self.edges[self.e] = [n_from, self.n]                         #add an edge connecting existing conduit node to current node
                                self.e = self.e+1                                             #increment edge counter up by one
                                #ax1.plot((self.nodes[self.n].x, self.nodes[n_from].x),(self.nodes[self.n].y, self.nodes[n_from].y))
                        self.n = self.n+1                                                   #increment node counter up by one
                    else:                                                                  #if there is NOT an outlet here
                        if p > 0:                                                           #if this is not the first point in the current path
                            #possible improvement: if the next point on the path is on an existing point, skip the current point.
                            self.nodes[self.n] = [point[0], point[1], 'junction']            #add a junction node here (with the node type for SWMM)
                            self.maps['nodes'][ix,iy] = self.n                               #update node map with node index
                            #ax1.scatter(point[1],point[0], marker='o', edgecolor='g', facecolor='none')   #debugging
                            if merge == False:                                              #if this conduit has not merged with an existing conduit
                                self.edges[self.e] = [self.n-1, self.n]                      #add and edge connecting the previous node to the current node
                                self.e = self.e+1                                            #increment edge counter up by one
                                #ax1.plot((self.nodes[self.n][1], self.nodes[self.n-1][1]),(self.nodes[self.n][0], self.nodes[self.n-1][0]), c='gold', marker=None)
                            else:                                                           #if this conduit HAS merged with an existing conduit
                                [[fromix,fromiy],error]  = self.fastMarching.IndexFromPoint(path[:,p-1]) #get xy indices of previous point in current conduit path
                                n_from = self.maps['nodes'][fromix,fromiy]                   #get node index of the node already in the cell where the previous point was
                                self.edges[self.e] = [n_from, self.n]                        #add an edge connecting existing conduit node to current node
                                self.e = self.e+1                                            #increment edge counter up by one
                                merge = False                                                #reset merge indicator to show that current conduit has left                                                              #if this is the first point in current path
                        else:                                                                #if this is the first point in the current path (counter <= 0, therefore it is an inlet)
                            self.nodes[self.n] = [point[0], point[1], 'inlet']               #add an inlet node here (with the node type for SWMM)
                            self.maps['nodes'][ix,iy] = self.n                               #update node map with node index
                            #ax1.scatter(point[1],point[0], marker='o', edgecolor='sienna', facecolor='none')
                        self.n = self.n+1                                                   #increment node counter up by one
                elif ~np.isnan(self.maps['nodes'][ix,iy]):                                 #if there is already a node in this cell (either because there is a conduit here, or because there are two nodes in the same cell)
                    n_existing = self.maps['nodes'][ix,iy]                                  #get index of node already present in current cell
                    if merge == True:                                                       #if this conduit has already merged into an existing conduit
                        pass                                                                 #skip this node (there is already a node here)
                    elif n_existing == self.n-1:                                            #if existing index is only one less than next node to be added index, this is a duplicate node and can be skipped
                        pass                                                                 #skip this node (duplicate)
                    else:                                                                   #if existing node index is >1 less than next node to be added index
                        if p > 0:                                                           #if this is not the first point in the current path
                            self.edges[self.e] = [self.n-1, n_existing]                      #add an edge connecting most recently added node and existing node in cell
                            self.e = self.e+1                                                #increment edge counter up by one
                            merge = True                                                     #add a flag indicating that this conduit has merged into an existing one
                            #ax1.plot((self.nodes[self.n-1][1], self.nodes[n_existing][1]),(self.nodes[self.n-1][0], self.nodes[n_existing][0]), c='r', marker=None)
                        else:                                                                #if this is the first point in the current path (i.e. the inlet is on an exising conduit)
                            self.nodes[self.n] = [point[0], point[1], 'inlet']                #add a node here (with the node type for SWMM)- this will cause there to be two nodes in the same cell
                            self.maps['nodes'][ix,iy] = self.n                                #update node map with node index
                            self.n = self.n+1                                                 #increment node counter by 1
                            #ax1.scatter(point[1],point[0], marker='o', edgecolor='g', facecolor='none')  #debugging

                self.maps['karst'][iteration][ix, iy] = 1                               #update karst map to put a conduit in current cell
        return None

    ###########################
    # Visualization functions #
    ###########################

    def show_catchment(self, data='geology', title=None, mask=False, cmap='binary'):
        """
        Show the entire study domain.
        """
        fig, ax1 = plt.subplots()
        if title is None:
            title = data
        fig.suptitle(title, fontsize=16)

        # Load data
        d = self.geology.data[data]['data']
        if mask==True:
            if self.mask is not None:
                d = self.geology_masked[data]
            else:
                return "Error : no mask to plot."
        if data in ["topography", "orientationx", "orientationy"]:
            d = np.transpose(d, (1,0)) # imshow read MxN and we have NxM
            im1 = ax1.imshow(d, extent=self.grid.extent, cmap=cmap, origin="lower")
        elif self.geology.data[data]["mode"] in ["image", "csv", "null"]:
            d = np.flipud(np.transpose(d, (1,0,2))) # we need to reverse transformations from geologymanager
            im1 = ax1.imshow(d , extent=self.grid.extent, cmap='gray_r')
        else:
            im1 = ax1.imshow(d, extent=self.grid.extent, cmap=cmap)#, origin="lower")
        fig.colorbar(im1, ax=ax1)
        if self.settings['data_has_polygon']:
            closed_polygon = self.polygon.polygon[:]
            closed_polygon.append(closed_polygon[0])
            x,y = zip(*closed_polygon)
            ax1.plot(x,y, color='red', label='polygon')
        for key in self.points.points:
            x,y = zip(*self.points.points[key])
            ax1.plot(x,y,'o',label=key)
        ax1.set_aspect('equal', 'box')
        plt.legend(loc='upper right')
        plt.show()
        return None

    def _show_maps(self, sim=-1, iteration=-1, cmap='binary'):
        """
        Show the simulated karst network as an image.
        """
        karst_network = self.karst_simulations[sim]

        fig, ([ax1,ax2],[ax3,ax4]) = plt.subplots(2, 2, sharex=True, sharey=True)
        fig.suptitle('Karst Network', fontsize=16)

        ax1.imshow(karst_network.maps['outlets'], extent=self.grid.extent, origin='lower', cmap=cmap)
        ax1.set_title('Outlets')

        ax2.imshow(karst_network.maps['cost'][iteration], extent=self.grid.extent, origin='lower', cmap=cmap)
        ax2.set_title('Cost')

        ax3.imshow(karst_network.maps['time'][iteration], extent=self.grid.extent, origin='lower', cmap=cmap)
        ax3.set_title('Time')

        ax4.imshow(karst_network.maps['karst'][iteration], extent=self.grid.extent, origin='lower', cmap=cmap)
        ax4.set_title('Karst')

        fig.subplots_adjust(hspace=0.5)
        plt.show()
        return None


    def show(self, data=None, title=None, probability=False):
        """
        Show the entire study domain (defaults to showing most recent simulation).
        """
        if data is None:
            data = self.karst_simulations[-1]

        if probability == True:
            data = self._compute_average_paths()

        fig = plt.figure(figsize=(20,10))

        # Cost map
        fig.add_subplot(131, aspect='equal')
        d = data.maps['cost'][-1]
        plt.xlabel('Cost array'+str(d.shape))
        d = np.transpose(d, (1,0)) # imshow read MxN and we have NxM
        plt.imshow(d, extent=self.grid.extent, origin='lower', cmap='gray') #darker=slower
        plt.colorbar(shrink=0.35)

        # Travel time map
        fig.add_subplot(132, aspect='equal')
        d = data.maps['time'][-1]
        plt.xlabel('Travel time array'+str(d.shape))
        d = np.transpose(d, (1,0)) # imshow read MxN and we have NxM
        plt.imshow(d, extent=self.grid.extent, origin='lower', cmap='cividis') #darker=faster
        plt.colorbar(shrink=0.35)

        # Karst map
        fig.add_subplot(133, aspect='equal')
        d = data.maps['karst'][-1]
        plt.xlabel('Karst array'+str(d.shape))
        d = np.transpose(d, (1,0)) # imshow read MxN and we have NxM
        plt.imshow(d, extent=self.grid.extent, origin='lower', cmap='gray_r') #darker=conduits
        plt.colorbar(shrink=0.35)
        i = plt.scatter(data.points['inlets'].x,  data.points['inlets'].y,  c='orange')
        o = plt.scatter(data.points['outlets'].x, data.points['outlets'].y, c='steelblue')
        p = matplotlib.patches.Rectangle((0,0),0,0, ec='r', fc='none')
        if self.settings['data_has_polygon']:
            closed_polygon = self.polygon.polygon[:]
            closed_polygon.append(closed_polygon[0])
            x,y = zip(*closed_polygon)
            plt.plot(x,y, color='red', label='polygon')
        plt.legend([i,o,p], ['inlets', 'outlets', 'catchment'], loc='upper right')

        if title is not None:
            fig.suptitle(title, fontsize=16)
        plt.show()
        return None

    def show_network(self, data=None, simplify=False, ax=None, plot_nodes=True, polygon=True, labels=['inlets', 'outlets'], title=None, cmap=None, color='k', alpha=1, legend=True):
        """
        #Chloe: This is a new function that I use to create all the figures for the paper.
        Show the karst network as a graph with nodes and edges. Defaults to showing latest iteration.
        Inputs:
        data:   karst simulation object containing nodes, edges, points, etc. Can be obtained from self.karst_simulations[i]
        ax:     axis to plot on
        label:  None or list of strings ['nodes','edges','inlets','outlets'], indicating which components to label
        title:  string, title of plot
        cmap:   string, colormap to use when plotting
        color:  string, single color to use when plotting (cannot have both cmap and color)
        alpha:  float, opacity to plot with (1=opaque, 0=transparent)
        legend: True/False, whether to display legend
        plot_nodes:   True/False, whether to display nodes
        polygon: True/False, whether to display the bounding polygon
        """

        if ax == None:
            fig,ax = plt.subplots(figsize=(10,10))
            ax.set_aspect('equal')

        if data == None:
            data = self.karst_simulations[-1]

        if polygon == True:
            if self.settings['data_has_polygon']:
                closed_polygon = self.polygon.polygon[:]
                closed_polygon.append(closed_polygon[0])
                x,y = zip(*closed_polygon)
                ax.plot(x,y, color='maroon')
                p = matplotlib.lines.Line2D([0],[0], color='k')

        if simplify == True:
            nodes = data.network['nodes']   #get all nodes
            nodes_simple = data.network['karstnet'].graph_simpl.nodes  #get indices of only the nodes in the simplified graph
            nodes_simple = {key: nodes[key] for key in nodes_simple}   #make df of only the nodes in the simplified graph, for plotting
            edges = data.network['edges']   #get all edges
            edges_simple = data.network['karstnet'].graph_simpl.edges  #get only the edges in the simplified graph
            edges_simple = {i: edge for i,edge in enumerate(edges_simple)}   #make df of only the edges in the simplified graph, for p
            nodes = pd.DataFrame.from_dict(nodes_simple, orient='index', columns=['x','y','type']) #convert to pandas for easier plotting
            edges = pd.DataFrame.from_dict(edges_simple, orient='index', columns=['inNode','outNode'])
        else:
            nodes = pd.DataFrame.from_dict(data.network['nodes'], orient='index', columns=['x','y','type']) #convert to pandas for easier plotting
            edges = pd.DataFrame.from_dict(data.network['edges'], orient='index', columns=['inNode','outNode'])

        #Set up data for plotting:
        fromX = nodes.x.loc[edges.inNode]      #calculate coordinates for link start and end points
        fromY = nodes.y.loc[edges.inNode]
        toX   = nodes.x.loc[edges.outNode]
        toY   = nodes.y.loc[edges.outNode]

        #Plot nodes and edges:
        if plot_nodes:
            n = ax.scatter(nodes.y,              nodes.x,                  c='k',         alpha=alpha, s=5)  #scatterplot nodes
        i = ax.scatter(data.points['inlets'].x,  data.points['inlets'].y,  c='orange',    s=30) #scatterplot inlets
        o = ax.scatter(data.points['outlets'].x, data.points['outlets'].y, c='steelblue', s=30) #scatterplot outlets
        e = matplotlib.lines.Line2D([0],[0])                                                  #line artist for legend
        for ind in edges.index:                                                               #loop over edge indices
            if cmap is not None:
                ax.plot((fromY.iloc[ind], toY.iloc[ind]), (fromX.iloc[ind], toX.iloc[ind]), c=plt.cm.get_cmap(cmap)(ind/len(edges)), alpha=alpha)  #plot each edge, moving along color gradient to show order
            elif color is not None:
                ax.plot((fromY.iloc[ind], toY.iloc[ind]), (fromX.iloc[ind], toX.iloc[ind]), c=color, alpha=alpha)  #plot each edge in same color

        #Add labels:
        if labels == None:
            pass
        else:
            if 'nodes' in labels:                                         #label node indices
                for ind in nodes.index:                                   #loop over node indices
                    ax.annotate(str(ind), xy=(nodes.y[ind]-10, nodes.x[ind]))  #annotate slightly to left of each node
            if 'edges' in labels:
                for ind in edges.index:
                    ax.annotate(str(ind), xy=(edges.y[ind]-10, edges.x[ind]))  #annotate slightly to left of each edge
            if 'inlets' in labels:
                for index,inlet in data.points['inlets'].iterrows():
                    ax.annotate(str(int(inlet.outlet))+'-'+str(int(inlet.inlet_iteration)),  xy=(inlet.x-(6*self.grid.dx),  inlet.y))
            if 'outlets' in labels:
                for index,outlet in data.points['outlets'].iterrows():
                    ax.annotate(str(int(outlet.name)), xy=(outlet.x-(4*self.grid.dx), outlet.y))

        #Add legend & title:
        if legend:
            if plot_nodes:
                if plot_polygon:
                    ax.legend([i,o,n,e,p],['inlets','outlets','nodes','edges','polygon'])
                else:
                    ax.legend([i,o,n,e],['inlets','outlets','nodes','edges'])
            else:
                if plot_polygon:
                    ax.legend([i,o,e,p],['inlets','outlets','edges','polygon'])
                else:
                    ax.legend([i,o,e],['inlets','outlets','edges','polygon'])
        if title is not None:
            ax.set_title(title, fontsize=16)

        return None


    def _compute_average_paths(self, show=False):
        """
        Compute the mean of all the simulations.
        Chloe: I have not tested this with the new format.
        """
        karst_maps = []
        for karst_simulation in self.karst_simulations:
            karst_maps.append(karst_simulation.maps['karst'][-1])

        karst_prob = sum(karst_maps)/len(karst_maps)

        if self.mask is not None:
            karst_prob = np.ma.MaskedArray(karst_prob, self.mask)

        return karst_prob

    def compare_stats(self, data=None, mean=False):
        """
        Compare statistics between reference indicators and calculated networks.
        """

        cpd  = []
        cv_d = []
        cv_l = []
        o_e  = []
        l_e  = []
        spl  = []
        m_d  = []
        cvd  = []
        vars = [cpd,cv_d,cv_l,o_e,l_e,spl,m_d,cvd]

        if data == None:
            i = -1
        else:
            i = 0

        for karst_network in self.karst_simulations[i:]:
            stats_list = []
            stats_list.append(karst_network.stats['cpd'])
            stats_list.append(karst_network.stats['cv degree'])
            stats_list.append(karst_network.stats['cv length'])
            stats_list.append(karst_network.stats['orientation entropy'])
            stats_list.append(karst_network.stats['length entropy'])
            stats_list.append(karst_network.stats['aspl'])
            stats_list.append(karst_network.stats['mean degree'])
            stats_list.append(karst_network.stats['correlation vertex degree'])

            cpd.append(karst_network.stats['cpd'])
            cv_d.append(karst_network.stats['cv degree'])
            cv_l.append(karst_network.stats['cv length'])
            o_e.append(karst_network.stats['orientation entropy'])
            l_e.append(karst_network.stats['length entropy'])
            spl.append(karst_network.stats['aspl'])
            m_d.append(karst_network.stats['mean degree'])
            cvd.append(karst_network.stats['correlation vertex degree'])

        if mean == False:
            print('\n')
            print('STATS for modelisation')
            print('%-10s%-12s%-12s%-12s%-12s' % ('Var', 'Value', 'Min ref', 'Max ref', 'Result'))
            print(54*'-')
            for stat_calc, key in zip(stats_list, self.reference_statistics):
                if stat_calc > self.reference_statistics[key][1]:
                    result = 'out +'
                elif stat_calc < self.reference_statistics[key][0]:
                    result = 'out -'
                else:
                    result = 'IN'
                print('%-10s%-12s%-12s%-12s%-12s' % (key, round(stat_calc,3), round(self.reference_statistics[key][0],3), round(self.reference_statistics[key][1],3), result))

        else:
            print('\n')
            print('MEAN STATS')
            print('%-10s%-12s%-12s%-12s%-12s' % ('Var', 'Value', 'Min ref', 'Max ref', 'Result'))
            print(54*'-')
            for var, key in zip(vars, self.reference_statistics):
                mean = sum(var)/len(var)
                if mean > self.reference_statistics[key][1]:
                    result = 'out +'
                elif mean < self.reference_statistics[key][0]:
                    result = 'out -'
                else:
                    result = 'IN'
                print('%-10s%-12s%-12s%-12s%-12s' % (key, round(mean,3), round(self.reference_statistics[key][0],3), round(self.reference_statistics[key][1],3), result))
        return None
