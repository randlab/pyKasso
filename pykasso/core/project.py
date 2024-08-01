"""
Module defining the project class.
"""

### Internal dependencies
import os
import pickle
import datetime
import shutil
import logging
import importlib.util
from importlib.metadata import version
from collections import Counter
from collections.abc import Sequence

### External dependencies
import yaml

### Local dependencies
from pykasso.core.grid import Grid
from .._version import __version__
from pykasso.core._namespaces import (MISC_DIR_PATH,
                                      DEFAULT_PARAMETERS_FILENAME,
                                      DEFAULT_PROJECT_FILENAME,
                                      DEFAULT_LOG_FILENAME)


class Project(Sequence):
    """
    Class modeling a project.
    
    This class stores and defines all the variables needed for defining a
    project in pyKasso.
    
    Attributes
    ----------
    name : str
    
    description : str
    
    creation_date : str
    
    pykasso_version : str
        
    grid : Grid
        
    core : dict
        
    n_simulations : int
        
    simulations : list
        
    """
    
    def __init__(self,
                 grid_parameters: dict,
                 project_location: str = None,
                 example: str = None,
                 force: bool = False
                 ) -> None:
        """
        Initialize a pyKasso project.
        """
        ### Initialize values
        self.__name = None
        self.__description = None
        self.__creation_date = None
        self.__pykasso_version = None
        self.__grid = None
        self.__core = None  # TODO
        self.__n_simulations = None
        self.__simulations = None
        self.__dimension = None
        
        ### Prepare the input variables
    
        # Define some useful variables, paths and locations
        package_location = os.path.dirname(os.path.abspath(__file__))
        misc_dir_loc = package_location + MISC_DIR_PATH
        self._pckg_paths = {
            'package_location': package_location,
            'misc_dir_path': MISC_DIR_PATH,
            'misc_dir_location': misc_dir_loc,
        }
        current_location = os.getcwd()
        
        ### Clean and test provided variables
        
        # Test and clean 'project_location' parameter
        is_project_loc_defined = (project_location is not None)
        if is_project_loc_defined:
            project_location = project_location.lower()
            project_location = project_location.strip('/')
        
        # Test and clean 'example' parameter
        is_example_defined = (example is not None)
        if is_example_defined:
            example = example.lower()
            valid_examples = os.listdir(misc_dir_loc + 'cases/')
            if example not in valid_examples:
                msg = ("'example' argument value is not valid. Valid examples"
                       " names : {}".format(valid_examples))
                raise ValueError(msg)
             
        ### Handle input cases
       
        # Case 1 : 'project_location' and 'example' are not defined
        if not is_project_loc_defined and not is_example_defined:
            date = datetime.datetime.now().strftime("%Y%m%d")
            project_dir = 'pyKasso_project_' + date
        
        # Case 2 : only 'project_location' is defined
        elif is_project_loc_defined and not is_example_defined:
            project_dir = project_location
        
        # Case 3 : only 'example' is defined
        elif not is_project_loc_defined and is_example_defined:
            project_dir = example
        
        # Case 4 : 'project_location' and 'example' are both defined
        elif is_project_loc_defined and is_example_defined:
            project_dir = project_location
        
        ### Create pyKasso's project directories
        
        # Define the project subdirectories and filenames
        self.core = {
            'locations': {
                'root': current_location,
                'project': current_location + '/' + project_dir + '/',
            },
            'paths': {
                'project_dir': project_dir + '/',
                'inputs_dir': project_dir + '/inputs/',
                'outputs_dir': project_dir + '/outputs/',
            },
            'filenames': {
                'parameters': DEFAULT_PARAMETERS_FILENAME,
                'project': DEFAULT_PROJECT_FILENAME,
                'log': DEFAULT_LOG_FILENAME,
            }
        }
        
        # Copy the files from queried example
        if is_example_defined:
            source = misc_dir_loc + 'cases/' + example
            destination = current_location + '/' + project_dir
            shutil.copytree(source, destination, dirs_exist_ok=force)
        # Otherwise create the pykasso's project structure
        else:
            os.makedirs(self.core['paths']['project_dir'], exist_ok=force)
            os.makedirs(self.core['paths']['inputs_dir'], exist_ok=force)
            os.makedirs(self.core['paths']['outputs_dir'], exist_ok=force)
            
            # Copy the default parameters.yaml file from the misc directory
            # to the project directory
            source = misc_dir_loc + self.core['filenames']['parameters']
            shutil.copy2(source, self.core['paths']['inputs_dir'])
        
        ### Construct the dictionary used for settings comparison between
        ### actual and previous simulation, used for the memoization operation
        self.__name = project_dir.split('/')[-1]
        self.__description = ''
        date = datetime.datetime.now().strftime("%Y/%m/%d %H:%M")
        self.__creation_date = date
        self.__pykasso_version = __version__
        self.__n_simulations = 0
        self.__simulations = []
            
        ### Construct the memoization dictionary
        self._memoization = {
            'settings': {
                ### static features ###
                'domain': None,
                'geology': None,
                'faults': None,
                ### dynamic features ###
                'fractures': None,
                'outlets': None,
                'inlets': None,
            },
            'model': {}
        }
        
        ### Create the grid
        self.__grid = Grid(**grid_parameters)
        
        ### Create the project log file
        self._create_log_file()
        
        ### Create the project YAML file
        self._export_project_file()
        
        ### Determine if it's a 2D or 3D project
        self._define_project_dimension()
        
        return None
    
    def __len__(self):
        return self.n_simulations
    
    def __getitem__(self, i):
        data = self._get_simulation_data(i)
        return data
    
    def _create_log_file(self) -> None:
        """
        Create the project log file.
        """
        
        ##### Initialization #####
    
        # Reset logging module
        root = logging.getLogger()
        list(map(root.removeHandler, root.handlers))
        list(map(root.removeFilter, root.filters))
            
        # Set logging output formats
        logging_entry = ' %(name)-30s | %(levelname)-8s | %(message)s'

        # Set logging level
        self._logging_levels = {
            0: logging.DEBUG,
            1: logging.INFO,
            2: logging.WARNING,
            3: logging.ERROR,
            4: logging.CRITICAL,
        }
        
        ##### Set logging file #####
        
        # Set the path
        outputs_directory = self.core['paths']['project_dir']
        log_filename = self.core['filenames']['log']
        log_file = outputs_directory + '/' + log_filename

        # Create new logging file
        logging.basicConfig(filename=log_file,
                            encoding='utf-8',
                            level=logging.INFO,
                            filemode="w",
                            format=logging_entry)
        
        # Print pyKasso 'logo'
        self.logger = logging.getLogger("♠")
        file = self._pckg_paths['misc_dir_location'] + 'log_logo.txt'
        with open(file, "r") as f:
            lines = f.readlines()
        for line in lines:
            self.logger.info(line.strip('\n'))

        # Print basic information in the log
        project_loc = os.getcwd() + '/' + outputs_directory
        project_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.logger = logging.getLogger("☼")
        self.logger.info("Directory location : " + project_loc)
        self.logger.info("Date log creation  : " + project_date)
        self.logger.info("pyKasso version    : v." + __version__)
        self.logger.info("")
        
        # Print dependencies versions
        packages = [
            ('numpy', 'mandatory'),
            ('pandas', 'mandatory'),
            ('agd', 'mandatory'),
            ('matplotlib', 'optional'),
            ('pyvista', 'optional'),
            ('karstnet', 'optional'),
            ('networkx', 'optional')
        ]
        self.logger.info("Dependencies versions:")
        for (package, mode) in packages:
            result = importlib.util.find_spec(package)
            if result is not None:
                msg = ' - ' + package + ' : ' + version(package)
            else:
                msg = ' - ' + package + ' : not installed'
            self.logger.info(msg)
        self.logger.info("")
        
        # Print grid parameters
        params = self.grid.get_grid_parameters()
        self.logger.info("Grid parameters:")
        for (key, value) in params.items():
            self.logger.info(" - {} : {}".format(key, value))
        self.logger.info("")
            
        return None
    
    def _export_project_file(self) -> None:
        """Export the project attributes in a YAML file."""
        project = self._return_project_status()
        # Set the path
        outputs_directory = self.core['paths']['project_dir']
        project_filename = self.core['filenames']['project']
        project_file = outputs_directory + '/' + project_filename
        with open(project_file, 'w') as handle:
            yaml.dump(project, handle, sort_keys=False)
        return None
    
    def _return_project_status(self) -> dict:
        """Update and return the dictionary describing the project."""
        project = {
            'name': self.name,
            'description': self.description,
            'creation_date': self.creation_date,
            'pykasso_version': self.pykasso_version,
            'grid': self.grid.get_grid_parameters(),
            'core': self.core,
            'n_simulations': self.n_simulations,
            'simulations': self.simulations,
        }
        return project

    ###########################
    ### GETTERS and SETTERS ###
    ###########################

    @property
    def name(self) -> str:
        """Return the name of the project."""
        return self.__name
    
    @property
    def description(self) -> str:
        """Return the description of the project."""
        return self.__description
    
    @description.setter
    def description(self, text: str = '') -> None:
        """Set the description of the project."""
        self.__description = text
        return None
        
    @property
    def creation_date(self) -> str:
        """Return the date of the creation of the project."""
        return self.__creation_date
    
    @property
    def pykasso_version(self) -> str:
        """Return the version of pykasso used during the creation of the
        project.
        """
        return self.__pykasso_version
    
    @property
    def grid(self) -> Grid:
        """Return the grid of the project."""
        return self.__grid
    
    @property
    def n_simulations(self) -> int:
        """Return the number of simulations already computed within the
        project.
        """
        return self.__n_simulations
    
    @property
    def simulations(self) -> list:
        """Return a list containing the locations of the simulations already
        computed within the project.
        """
        return self.__simulations
    
    @simulations.setter
    def simulations(self, locations: list) -> None:
        """Set the list containing the locations of the simulations already
        computed within the project.
        """
        self.__simulations = locations
        return None
    
    @property
    def dimension(self) -> str:
        """
        TODO
        """
        return self.__dimension
    
    ##############
    ### Others ###
    ##############
    
    def _define_project_dimension(self) -> None:
        """Determine if the project is in 2D or 3D."""
        counter = Counter(self.grid.shape)
        if counter[1] == 1:
            self.__dimension = '2D'
        else:
            self.__dimension = '3D'
        return None
    
    def _get_simulation_data(self, n: int) -> dict:
        """Return the data of computed simulation ``n``."""
        simulation_directory = self.simulations[n]
        simulation_data_path = simulation_directory + 'results.pickle'
        simulation_data = self._read_pickle(simulation_data_path)
        return simulation_data
    
    def _read_pickle(self, path: str) -> dict:
        """Read a pickle from a given path."""
        with open(path, 'rb') as handle:
            out = pickle.load(handle)
            return out
        
    def _increment_n_simulations(self) -> None:
        self.__n_simulations = self.__n_simulations + 1
        return None
        
    def get_simulations(self, ni: int = 0, nf: int = None) -> list:
        """TODO"""
        list_sims = list(range(self.n_simulations))
        selected_sims = list_sims[ni:nf]
        return selected_sims
    
    def get_first_simulations(self, n: int) -> list:
        """TODO"""
        first_sims = self.get_simulations(0, n)
        return first_sims
    
    def get_last_simulations(self, n: int) -> list:
        """TODO"""
        last_sims = self.get_simulations(-n, None)
        return last_sims
