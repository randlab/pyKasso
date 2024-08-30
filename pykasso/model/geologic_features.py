"""
This module contains classes modeling the geological features.
"""

### External dependencies
import numpy as np
import pandas as pd

### Local dependencies
from pykasso._utils.datareader import DataReader
from pykasso.core._namespaces import (
    DEFAULT_FMM_COSTS,
    DEFAULT_FEATURE_PARAMETERS,
)

### Typing
from typing import Union
from pykasso.core.grid import Grid

##################
### Main class ###
##################


class GeologicFeature(DataReader):
    """
    Class modeling a geological feature.
    """
    
    def __init__(
        self,
        grid: Grid,
        feature: str,
        dim: int,
        *args,
        **kwargs,
    ) -> None:
        """
        Construct a geological feature.

        Parameters
        ----------
        grid : Grid
            TODO
        feature : str
            TODO
        dim : int
            TODO
        """
        
        # Initialization
        super().__init__(grid, *args, **kwargs)
        self.feature = feature
        self.dim = dim
        self.data_surface = None
        self.data_volume = None
        self.stats = None
        self.names = None
        self.costs = None
        self.model = None
        
        # Retrieve arguments from kwargs
        data = kwargs.get('data', None)
        axis = kwargs.get('axis', 'z')
        names = kwargs.get('names', {})
        costs = kwargs.get('costs', {})
        model = kwargs.get('model', {})
        
        # Set the data
        self.set_data(data, axis)
        self.compute_statistics()
        if self.feature not in ['topography', 'water_table', 'bedrock']:
            self.set_names(names)
            self.set_costs(costs)
            self.set_model(model)
            
    def overview(self) -> pd.DataFrame:
        """
        Return a pandas DataFrame describing each contained unit with its name,
        its cost, and if it will be considered during the simulation. Basic
        statistics are also described.
        """
        index = self.stats.index
        data = {
            'names': self.names.values(),
            'costs': self.costs.values(),
            'model': self.model.values(),
        }
        df = pd.DataFrame(data, index=index)
        df = pd.merge(df, self.stats, left_index=True, right_index=True)
        return df
            
    ###############
    ### SETTERS ###
    ###############
        
    def set_data(
        self,
        data: Union[None, str, np.ndarray],
        axis: str = 'z',
    ) -> None:
        """TODO"""
        # If no data is provdided
        if data is None:
            
            value = DEFAULT_FEATURE_PARAMETERS[self.feature]['nodata']

            if self.dim == 2:
                self.data_surface = self._get_data_full_2D(value)
            elif self.dim == 3:
                self.data_volume = self._get_data_full_3D(value)
        
        # Else
        elif isinstance(data, np.ndarray):
            if self.dim == 2:
                self.data_surface = data
            elif self.dim == 3:
                self.data_volume = data
        else:
            if self.dim == 2:
                self.data_surface = self.get_data_from_file(data,
                                                            False)
            elif self.dim == 3:
                self.data_volume = self.get_data_from_file(data,
                                                           True,
                                                           axis)
        return None
    
    def set_names(
        self,
        names: dict[int, str],
        default_name: str = 'item {}',
    ) -> None:
        """
        Assign names to items based on the provided ``names`` dictionary, with
        an optional default naming pattern.

        Parameters
        ----------
        names : dict[int, str]
            A dictionary where the keys are item indices (integers) and the
            values are the corresponding names (strings) to be assigned. This
            dictionary specifies which items should receive custom names.
        default_name : str, optional
            A format string used to generate default names for items not
            explicitly named in the ``names`` dictionary. The format string
            should include a placeholder (e.g., '{}') that will be replaced by
            the item's index. The default pattern is 'item {}'.
        
        Notes
        -----
        This function does not return a value. It rewrites the ``self.names``
        attribute with the new specified dictionary.
        """
        ids = self.stats.index
        names_df = {}
        for id in ids:
            names_df[id] = names.get(id, default_name.format(id))
        self.names = names_df
        return None
        
    def set_costs(
        self,
        costs: dict[int, str],
        default_cost: float = 0.5,
    ) -> None:
        """
        Assign costs to items based on the provided dictionary, with an
        optional default cost.

        Parameters
        ----------
        costs : dict[int, str]
            A dictionary where the keys are item indices (integers) and the
            values are the corresponding costs (floats) to be assigned. This
            dictionary specifies which items should receive custom costs.
        default_cost : float, optional
            The default cost to be applied to items not explicitly listed in
            the `costs` dictionary. The default value is 0,5.
        
        Notes
        -----
        This function does not return a value. It rewrites the ``self.costs``
        attribute with the new specified dictionary.
        """
        ids = self.stats.index
        costs_df = {}
        for id in ids:
            costs_df[id] = costs.get(id, default_cost)
        self.costs = costs_df
        return None
    
    def set_model(
        self,
        model: dict[int, str],
        default_model: bool = True,
    ) -> None:
        """
        Indicate if an item should be considered in the modelisation based on
        the provided dictionary, with an optional default setting.

        Parameters
        ----------
        model : dict[int, bool]
            A dictionary where the keys are item indices (integers) and the
            values are booleans indicating if the item is considered or not.
        default_model : bool, optional
            The default value to be applied to items not explicitly listed in
            the `model` dictionary. The default value is `True`.
        
        Notes
        -----
        This function does not return a value. It rewrites the ``self.model``
        attribute with the new specified dictionary.
        """
        ids = self.stats.index
        model_df = {}
        for id in ids:
            model_df[id] = model.get(id, default_model)
        self.model = model_df
        return None
    
    ###############
    ### GETTERS ###
    ###############
    
    def get_data_units(self, units: list[int]) -> np.ndarray:
        """
        Return a copy of the ``self.data_volume`` attribute only containing
        the specified units.
        
        Parameters
        ----------
        units: list[int]
            List of units to retrieve.
        
        Returns
        -------
        np.ndarray
        """
        data = np.empty(self.grid.shape) * np.nan
        test = np.isin(self.data_volume, units)
        data = np.where(test, self.data_volume, data)
        return data
    
    def get_data_model(self) -> np.ndarray:
        """
        Return a copy of the ``self.data_volume`` attribute corresponding of
        the state of the ``self.model`` attribute.
        """
        valid_ids = [id_ for (id_, boolean) in self.model.items() if boolean]
        geology = self.get_data_units(valid_ids)
        return geology
    
    #############
    ### OTHER ###
    #############
    
    def compute_statistics(self) -> None:
        """
        Populate the ``self.stats`` attribute with a pandas DataFrame
        containing statistics (counts and frequency) on the data.
        
        Returns
        -------
        None
        """
        values, counts = np.unique(self.data_volume, return_counts=True)
        values = values.astype('int')
        stats = {
            'counts': counts,
            'freq': counts / self.grid.nodes,
            'volume': counts * self.grid.node_volume,
        }
        self.stats = pd.DataFrame(data=stats, index=values)
        return None
    

#######################
### Main subclasses ###
#######################


class Surface(GeologicFeature):
    """
    Subclass modeling a two dimensional geological feature.
    """
    
    def __init__(
        self,
        grid: Grid,
        feature: str,
        *args,
        **kwargs,
    ) -> None:
        dim = 2
        super().__init__(grid, feature, dim, *args, **kwargs)
    
    def _surface_to_volume(self, condition: str, grid: Grid) -> np.ndarray:
        """
        Convert a two dimensional array in a three dimensional array.
        """
        k = grid.get_k(self.data_surface)
        data_volume = np.zeros((grid.nx, grid.ny, grid.nz))
        for z in range(grid.nz):
            data_volume[:, :, z] = z
            if condition == '>=':
                test = data_volume[:, :, z] >= k
            elif condition == '=':
                test = data_volume[:, :, z] == k
            elif condition == '<=':
                test = data_volume[:, :, z] <= k
            data_volume[:, :, z] = np.where(test, 1, 0)
        return data_volume
    
    def compute_statistics(self) -> None:
        """
        Populate the ``self.stats`` attribute with a pandas DataFrame
        containing statistics (counts and frequency) on the data.
        
        Returns
        -------
        None
        """
        values, counts = np.unique(self.data_surface, return_counts=True)
        values = values.astype('int')
        stats = {
            'counts': counts,
            'freq': counts / self.grid.nodes,
            'surface': counts * self.grid.node_area,
        }
        self.stats = pd.DataFrame(data=stats, index=values)
        return None
    

#####################
### 3D Subclasses ###
#####################


class Geology(GeologicFeature):
    """
    Class modeling the geologic model.
    """
    
    def __init__(
        self,
        grid: Grid,
        *args,
        **kwargs,
    ) -> None:
        feature = 'geology'
        dim = 3
        super().__init__(grid, feature, dim, *args, **kwargs)
        
    def set_names(
        self,
        names: dict[int, str],
        default_name: str = 'unit {}',
    ) -> None:
        """
        Assign names to geologic units based on the provided ``names``
        dictionary, with an optional default naming pattern.

        Parameters
        ----------
        names : dict[int, str]
            A dictionary where the keys are geologic unit indices (integers)
            and the values are the corresponding names (strings) to be
            assigned. This dictionary specifies which geologic unit should
            receive custom names.
        default_name : str, optional
            A format string used to generate default geologic unit names for
            items not explicitly named in the ``names`` dictionary. The format
            string should include a placeholder (e.g., '{}') that will be
            replaced by the item's index. The default pattern is 'family {}'.
        
        Notes
        -----
        This function does not return a value. It rewrites the ``self.names``
        attribute with the new specified dictionary.
        """
        return super().set_names(names, default_name)
    
    def set_costs(
        self,
        costs: dict[int, str],
        default_cost: float = DEFAULT_FMM_COSTS['geology'],
    ) -> None:
        """
        Assign costs to geologic units based on the provided dictionary, with an
        optional default cost.

        Parameters
        ----------
        costs : dict[int, str]
            A dictionary where the keys are geologic unit indices (integers)
            and the values are the corresponding costs (floats) to be assigned.
            This dictionary specifies which geologic units should receive
            custom costs.
        default_cost : float, optional
            The default cost to be applied to geologic units not explicitly
            listed in the `costs` dictionary. The default values are taken
            from the `DEFAULT_FMM_COSTS['geology']` dictionary.
        
        Notes
        -----
        This function does not return a value. It rewrites the ``self.costs``
        attribute with the new specified dictionary.
        """
        return super().set_costs(costs, default_cost)
    
    def set_model(
        self,
        model: dict[int, str],
        default_model: bool = True,
    ) -> None:
        """
        Indicate if a geologic unit should be considered in the modelisation
        based on the provided dictionary, with an optional default setting.

        Parameters
        ----------
        model : dict[int, bool]
            A dictionary where the keys are geologic unit indices (integers)
            and the values are booleans indicating if the item is considered or
            not.
        default_model : bool, optional
            The default value to be applied to geologic units not explicitly
            listed in the `model` dictionary. The default value is `True`.
        
        Notes
        -----
        This function does not return a value. It rewrites the ``self.model``
        attribute with the new specified dictionary.
        """
        model.setdefault(0, True)
        return super().set_model(model, default_model)
        

class Faults(GeologicFeature):
    """
    Class modeling the faults model.
    """
    
    def __init__(
        self,
        grid: Grid,
        *args,
        **kwargs,
    ) -> None:
        feature = 'faults'
        dim = 3
        super().__init__(grid, feature, dim, *args, **kwargs)
        
    def set_names(
        self,
        names: dict[int, str],
        default_name: str = 'fault {}',
    ) -> None:
        """
        Assign names to fault items based on the provided ``names`` dictionary
        , with an optional default naming pattern.

        Parameters
        ----------
        names : dict[int, str]
            A dictionary where the keys are fault item indices (integers) and
            the values are the corresponding names (strings) to be assigned.
            This dictionary specifies which geologic unit should receive
            custom names.
        default_name : str, optional
            A format string used to generate default fault item names for
            items not explicitly named in the ``names`` dictionary. The format
            string should include a placeholder (e.g., '{}') that will be
            replaced by the item's index. The default pattern is 'family {}'.
        
        Notes
        -----
        This function does not return a value. It rewrites the ``self.names``
        attribute with the new specified dictionary.
        """
        return super().set_names(names, default_name)

    def set_costs(
        self,
        costs: dict[int, str],
        default_cost: float = DEFAULT_FMM_COSTS['faults'],
    ) -> None:
        """
        Assign costs to fault items based on the provided dictionary, with an
        optional default cost.

        Parameters
        ----------
        costs : dict[int, str]
            A dictionary where the keys are fault item indices (integers) and
            the values are the corresponding costs (floats) to be assigned.
            This dictionary specifies which geologic units should receive
            custom costs.
        default_cost : float, optional
            The default cost to be applied to fault items not explicitly
            listed in the `costs` dictionary. The default values are taken
            from the `DEFAULT_FMM_COSTS['faults']` dictionary.
        
        Notes
        -----
        This function does not return a value. It rewrites the ``self.costs``
        attribute with the new specified dictionary.
        """
        return super().set_costs(costs, default_cost)
    
    def set_model(
        self,
        model: dict[int, str],
        default_model: bool = True,
    ) -> None:
        """
        Indicate if a fault item should be considered in the modelisation
        based on the provided dictionary, with an optional default setting.

        Parameters
        ----------
        model : dict[int, bool]
            A dictionary where the keys are fault item indices (integers) and
            the values are booleans indicating if the item is considered or
            not.
        default_model : bool, optional
            The default value to be applied to fault items not explicitly
            listed in the `model` dictionary. The default value is `True`.
        
        Notes
        -----
        This function does not return a value. It rewrites the ``self.model``
        attribute with the new specified dictionary.
        """
        model.setdefault(0, False)
        return super().set_model(model, default_model)

# DOCSTRING
# Surface.compute_statistics.__doc__ = GeologicFeature.compute_statistics.__doc__
