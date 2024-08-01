"""
Module containing the application class.
"""

### Internal dependencies
import warnings
import logging

### External dependencies

### Local dependencies
from pykasso.core.project import Project
from pykasso.model.sks import SKS
from pykasso.analysis.analysis import Analyzer
from pykasso.visualization.visualizer import Visualizer

### Validation
from pykasso._utils.validation import (
    is_variable_type_valid,
    is_key_in_dict,
    is_parameter_comparison_valid,
)

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
        Project class.
    model
        Model class.
    analyzer
        Analyzer class.
    visualizer
        Visualizer class
        
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
                    name: str,
                    grid_parameters: dict,
                    force: bool = True
                    ) -> None:
        """
        Create a new project.
        
        Instance a ``Project`` class within the ``project`` attribute and
        initialize the subpackages.

        Parameters
        ----------
        name : str
            The name of the project. A new directory is created if the
            argument points to a non-existant folder.
        grid_parameters : dict
            The dictionary containing the grid parameters.
        force : bool, optional
            If True, overwrite files in case of conflict when ``name``
            points to an already existing directory.
            Default is True.

        Examples
        --------
        >>> import pykasso as pk
        >>> app = pk.pykasso()
        >>> name = TODO
        >>> grid_parameters = {TODO}
        >>> app.new_project(TODO)
        """
        
        ### Input validation
        
        # Test 'name' type
        try:
            is_variable_type_valid(variable_name='name',
                                   variable_value=name,
                                   valid_types=(str))
        except TypeError:
            raise
        
        # Test 'grid_parameters' type
        try:
            is_variable_type_valid(variable_name='grid_parameters',
                                   variable_value=grid_parameters,
                                   valid_types=(dict))
        except TypeError:
            raise

        # Test 'Grid' parameters presence
        for parameter in GRID_PARAMETERS:
            try:
                is_key_in_dict(dictionary=grid_parameters,
                               dictionary_name='grid_parameters',
                               key_to_test=parameter)
            except KeyError:
                raise
            
        # Test if the values of attributes are of type int or float
        for parameter_name in ['x0', 'y0', 'z0', 'dx', 'dy', 'dz']:
            try:
                parameter_value = grid_parameters[parameter_name]
                is_variable_type_valid(variable_name=parameter_name,
                                       variable_value=parameter_value,
                                       valid_types=(int, float))
            except TypeError:
                raise

        # Test if the values of attributes are of type int
        for parameter_name in ['nx', 'ny', 'nz']:
            try:
                parameter_value = grid_parameters[parameter_name]
                is_variable_type_valid(variable_name=parameter_name,
                                       variable_value=parameter_value,
                                       valid_types=(int))
            except TypeError:
                raise

        # Test if the values of attributes are well upper 0
        for parameter_name in ['nx', 'ny', 'nz']:
            try:
                parameter_value = grid_parameters[parameter_name]
                is_parameter_comparison_valid(parameter_name=parameter_name,
                                              parameter_value=parameter_value,
                                              logical_test='>',
                                              compared_to=0)
            except ValueError:
                raise
            
        ### Initialization of the application
        
        # Set a project instance
        self.__project = Project(grid_parameters=grid_parameters,
                                 project_location=name,
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
