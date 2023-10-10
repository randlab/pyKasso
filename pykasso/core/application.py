"""
Module containing the application class.
"""

### Internal dependencies
import warnings

### External dependencies

### Local dependencies
from pykasso.core.project import Project
from pykasso.model.sks import SKS
from pykasso.analysis.analysis import Analyzer
from pykasso.visualization.visualizer import Visualizer

### Validation
from pykasso._utils.validation import (is_parameter_type_valid,
                                       is_key_in_dict,
                                       is_parameter_value_valid)

### Variables
from pykasso.core._namespaces import (GRID_PARAMETERS)


class Application():
    """
    Class modeling an application and embedding the pyKasso subpackages.
     
    This class manages a pyKasso project and provides access to the different
    subpackages by storing them as class attributes.
    
    Attributes
    ----------
    project
        Project class. TODO
    model
        Model class. TODO
    analyzer
        Analyzer class. TODO
    visualizer
        Visualizer class TODO
        
    Notes
    -----
    The attributes are set to ``None`` until a project is created or loaded.
        
    Examples
    --------
    This class can be instancied by using the public function ``pykasso()``.
    >>> import pykasso as pk
    >>> app = pk.pykasso()
    """
    
    def __init__(self) -> None:
        self.__project = None
        self.__model = None
        self.__analyzer = None
        self.__visualizer = None
        
    ######################
    ### MANAGE PROJECT ###
    ######################
        
    def new_project(self,
                    project_name: str,
                    grid_parameters: dict,
                    force: bool = True
                    ) -> None:
        """
        Create a new project.
        
        Instance a ``Project`` class within the ``project`` attribute and
        initialize the subpackages.

        Parameters
        ----------
        project_name : str
            The name of the project. A new directory is created if the
            argument points to a non-existant folder.
        grid_parameters : dict
            The dictionary containing the grid parameters.
        force : bool, optional
            If True, overwrite files in case of conflict when ``project_name``
            points to an already existing directory.
            Default is True.

        Examples
        --------
        >>> import pykasso as pk
        >>> app = pk.pykasso()
        >>> project_name = TODO
        >>> grid_parameters = {TODO}
        >>> app.new_project(TODO)
        """
        
        ### Input validation
        
        # Test 'project_name' type
        is_parameter_type_valid(parameter_name='project_name',
                                parameter_value=project_name,
                                valid_types=(str))
        
        # Test 'grid_parameters' type
        is_parameter_type_valid(parameter_name='grid_parameters',
                                parameter_value=grid_parameters,
                                valid_types=(dict))

        # Test 'Grid' parameters presence
        for parameter in GRID_PARAMETERS:
            is_key_in_dict(dictionary=grid_parameters,
                           dictionary_name='grid_parameters',
                           key_to_test=parameter)
            
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
    
    def open_project(self) -> NotImplementedError:
        """
        Not implemented yet.
        """
        msg = "Not implemented yet."
        raise NotImplementedError(msg)
    
    def save_project(self) -> NotImplementedError:
        """
        Not implemented yet.
        """
        msg = "Not implemented yet."
        raise NotImplementedError(msg)
    
    def export_project(self) -> NotImplementedError:
        """
        Not implemented yet.
        """
        msg = "Not implemented yet."
        raise NotImplementedError(msg)

    ###############
    ### GETTERS ###
    ###############

    @property
    def project(self) -> Project:
        """
        Return the project class.
        """
        if self.__project is None:
            msg = "No project available yet. Please create or load a project."
            print(msg)
            return None
        else:
            return self.__project
        
    @property
    def model(self) -> SKS:
        """
        Return the SKS model class.
        """
        if self.__model is None:
            msg = ("This feature is not available yet. Please create or load a"
                   " project first.")
            print(msg)
            return None
        else:
            return self.__model
        
    @property
    def analyzer(self) -> Analyzer:
        """
        Return the analyzer class.
        """
        if self.__analyzer is None:
            msg = ("This feature is not available yet. Please create or load a"
                   " project first.")
            print(msg)
            return None
        else:
            return self.__analyzer
        
    @property
    def visualizer(self) -> Visualizer:
        """
        Return the visualizer class.
        """
        if self.__visualizer is None:
            msg = ("This feature is not available yet. Please create or load a"
                   " project first.")
            print(msg)
            return None
        else:
            return self.__visualizer
