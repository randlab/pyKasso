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

### External dependencies
import yaml

### Local dependencies
from pykasso.core.grid import Grid
from .._version import __version__


class Project:
    """Defining a pyKasso project in the application class.
    
    This class stores and defines all the variables needed for defining a
    project in pyKasso.
    
    Attributes
    ----------
    TODO
        -
    TODO
        -
    
    """
    
    def __init__(self, grid_parameters: dict, project_location: str = None,
                 example: str = None, force: bool = False) -> None:
        """
        Initialize a pyKasso project.

        Parameters
        ----------
        grid_parameters : dict
            _description_
        project_location : str, optional
            _description_, by default None
        example : str, optional
            _description_, by default None
        force : bool, optional
            _description_, by default False

        Returns
        -------
        _type_
            _description_

        Raises
        ------
        ValueError
            _description_
        """
        
        ### Prepare the input variables
    
        # Define some useful variables, paths and locations
        package_location = os.path.dirname(os.path.abspath(__file__))
        misc_dir_path = '\\..\\_misc\\'
        misc_dir_loc = package_location + misc_dir_path
        self._pckg_paths = {
            'package_location': package_location,
            'misc_dir_path': misc_dir_path,
            'misc_dir_location': misc_dir_loc,
        }
        current_location = os.getcwd()
        
        # Test type of parameters
        is_project_loc_defined = (project_location is not None)
        is_example_provided = (example is not None)
        
        ### Clean and test provided variables
        
        # Clean 'project_location' parameter
        if is_project_loc_defined:
            project_location = project_location.lower()
            project_location = project_location.replace('/', '\\')
            project_location = project_location.strip('\\')
        # Test validity of 'example' parameter
        if is_example_provided:
            example = example.lower()
            valid_examples = os.listdir(misc_dir_loc + 'cases\\')
            if example not in valid_examples:
                msg = ("'example' argument value is not valid. Valid examples"
                       "names : {}".format(valid_examples))
                raise ValueError(msg)
             
        ### Handle input cases
       
        # Case 1 : 'project_location' and 'example' are not provided
        if not is_project_loc_defined and not is_example_provided:
            date = datetime.datetime.now().strftime("%Y%m%d")
            project_dir = 'pyKasso_project_' + date
        
        # Case 2 : only 'project_location' is provided
        elif is_project_loc_defined and not is_example_provided:
            project_dir = project_location
        
        # Case 3 : only 'example' is provided
        elif not is_project_loc_defined and is_example_provided:
            project_dir = example
        
        # Case 4 : 'project_location' and 'example' are both provided
        elif is_project_loc_defined and is_example_provided:
            project_dir = project_location
        
        ### Install pyKasso's project directories
        
        # Define the project subdirectories and filenames
        self.core = {
            'locations': {
                'root': current_location,
                'project': current_location + '\\' + project_dir + '\\',
            },
            'paths': {
                'project_dir': project_dir + '\\',
                'inputs_dir': project_dir + '\\inputs\\',
                'outputs_dir': project_dir + '\\outputs\\',
            },
            'filenames': {
                'parameters': 'parameters.yaml',
                'project': 'project.yaml',
                'log': 'project.log',
            }
        }
        
        # Copy the files from queried example
        if is_example_provided:
            source = misc_dir_loc + 'cases\\' + example
            destination = current_location + '\\' + project_dir
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
        self.name = project_dir.split('\\')[-1]
        self.description = ''
        self.creation_date = datetime.datetime.now().strftime("%Y/%m/%d %H:%M")
        self.pykasso_version = __version__
        self.n_simulations = 0
        self.simulations = []
            
        ### Construct the memoization dictionary
        self._create_memoization_dict()
        
        ### Create the project log file
        self._create_log_file()
        
        ### Create the grid
        self._create_grid(grid_parameters)
        
        ### Create the project YAML file
        self._export_project_file()
        
        return None
    
    def __repr__(self) -> str:
        return str(self._return_project_status())
    
    # TODO
    def __str__(self) -> str:
        txt = ("Project"
               "\n - Name : "
               "\n - Description : ") 
    #            "\n[dx, dy, dz] : ({}, {}, {})"
    #            .format(self.x0, self.y0, self.z0,
    #                    self.nx, self.ny, self.nz,
    #                    self.dx, self.dy, self.dz))
        return txt
    
    def _create_memoization_dict(self) -> None:
        """
        Dictionary used for memory optimization/memoizing operations
        """
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
        return None
        
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
        log_file = outputs_directory + '\\' + log_filename

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
        project_loc = os.getcwd() + '\\' + outputs_directory
        project_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.logger = logging.getLogger("☼")
        self.logger.info("Directory location : " + project_loc)
        self.logger.info("Date log creation  : " + project_date)
        self.logger.info("pyKasso version    : v." + __version__)
        self.logger.info("")
        
        # Output dependencies versions
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
            
        return None
    
    def _create_grid(self, grid_parameters: dict) -> None:
        """Create and declare a Grid instance."""
        self.grid = Grid(**grid_parameters)
        return None
    
    def _export_project_file(self) -> None:
        """Export the project attributes in a YAML file."""
        project = self._return_project_status()
        # Set the path
        outputs_directory = self.core['paths']['project_dir']
        project_filename = self.core['filenames']['project']
        project_file = outputs_directory + '\\' + project_filename
        with open(project_file, 'w') as handle:
            yaml.dump(project, handle, sort_keys=False)
        return None
    
    def _return_project_status(self) -> dict:
        """TODO"""
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

    def _read_pickle(self, path: str) -> dict:
        """Read a pickle from a given path."""
        with open(path, 'rb') as handle:
            out = pickle.load(handle)
            return out
        
    def _get_simulation_data(self, n: int) -> dict:
        """Get the n computed karst network simulation."""
        simulation_directory = self.simulations[n]
        simulation_data_path = simulation_directory + 'results.pickle'
        simulation_data = self._read_pickle(simulation_data_path)
        return simulation_data
