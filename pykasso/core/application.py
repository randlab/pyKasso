"""
Module containing the application class.
"""

### Internal dependencies

### External dependencies

### Local dependencies
from pykasso.core.project import Project
from pykasso.model.sks import SKS
from pykasso.analysis.analysis import Analyzer
from pykasso.visualization.visualizer import Visualizer

from pykasso._utils.validation import (is_parameter_type_valid,
                                       is_key_in_dict,
                                       is_parameter_value_valid)


class Application():
    """Embedding of the pyKasso subpackages.
     
    This class provides the capability to access to the different pyKasso
    subpackages.
    
    Attributes
    ----------
    project
        -
    model
        -
    analyzer
        -
    visualizer
        -
        
    Examples
    --------
    >>> import pykasso as pk
    >>> app = pk.pykasso()
    """
    
    def __init__(self) -> None:
        """Initialize the application."""
        self.__project = None
        self.__model = None
        self.__analyzer = None
        self.__visualizer = None
        
    ######################
    ### MANAGE PROJECT ###
    ######################
        
    def new_project(self, project_name: str, grid_parameters: dict,
                    force: bool = True) -> None:
        """Create a new project.

        Parameters
        ----------
        project_name : str
            The name of the project.
        grid_parameters : dict
            The grid parameters.
        force : bool, optional
            TODO, by default True

        Examples
        --------
        >>> import pykasso as pk
        >>> app = pk.pykasso()
        >>> grid = {TODO}
        >>> app.new_project(TODO)
        """
        
        ### Input validation
        
        # Parameters type
        is_parameter_type_valid(project_name, 'project_name', (str))
        is_parameter_type_valid(grid_parameters, 'grid_parameters', (dict))
        
        # Grid parameters presence
        for p in ['x0', 'y0', 'z0', 'nx', 'ny', 'nz', 'dx', 'dy', 'dz']:
            is_key_in_dict(grid_parameters, 'grid_parameters', p)
            
        # Grid parameters values
        for p in ['x0', 'y0', 'z0', 'dx', 'dy', 'dz']:
            is_parameter_type_valid(grid_parameters[p], p, (int, float))
        for p in ['nx', 'ny', 'nz']:
            is_parameter_type_valid(grid_parameters[p], p, (int))
        for p in ['nx', 'ny', 'nz']:
            is_parameter_value_valid(grid_parameters[p], p, '>', 0)
        
        ### Initialization of the application
        
        # Set a project instance
        self.__project = Project(grid_parameters=grid_parameters,
                                 project_location=project_name,
                                 force=force)
        
        # Initialize the 'model' module
        self.__model = SKS(self.project)
        
        # Initialize the 'analysis' module
        self.__analyzer = Analyzer(self.project)
        
        # Initialize the 'visualisation' module
        self.__visualizer = Visualizer(self.project)
        
        return None
    
    def open_project(self) -> None:
        return None
    
    def save_project(self) -> None:
        return None
    
    def export_project(self) -> None:
        return None

    ###############
    ### GETTERS ###
    ###############

    @property
    def project(self):
        """Return the project class."""
        if self.__project is None:
            txt = "No project available yet. Please create or load a project."
            print(txt)
            return None
        else:
            return self.__project
        
    @property
    def model(self):
        """Return the model class."""
        if self.__model is None:
            txt = ("This feature is not available yet. Please create or load a"
                   " project first.")
            print(txt)
            return None
        else:
            return self.__model
        
    @property
    def analyzer(self):
        """Return the analyzer class."""
        if self.__analyzer is None:
            txt = ("This feature is not available yet. Please create or load a"
                   " project first.")
            print(txt)
            return None
        else:
            return self.__analyzer
        
    @property
    def visualizer(self):
        """Return the visualizer class."""
        if self.__visualizer is None:
            txt = ("This feature is not available yet. Please create or load a"
                   " project first.")
            print(txt)
            return None
        else:
            return self.__visualizer
