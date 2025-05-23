"""
This module contains classes modeling the geological features.
"""

### External dependencies
import numpy as np
import pandas as pd

### Local dependencies
from ..._utils.datareader import DataReader
from ...core._namespaces import DEFAULT_FEATURE_PARAMETERS

### Typing
from typing import Union
from ...core.grid import Grid


class GeologicFeature(DataReader):
    """
    Class modeling a geological feature.
    """
    
    def __init__(
        self,
        grid: Grid,
        feature: str,
        dim: int,
        default_fmm_cost: float,
        *args,
        **kwargs,
    ) -> None:
        """
        Construct a geological feature.

        Parameters
        ----------
        grid : Grid
            pyKasso's ``Grid`` of the model.
        feature : str
            Define the type of geological feature.
            
            Available 2D geological features:
                - ``'topography'``
                - ``'water_table'``
                - ``'bedrock'``

            Available 3D geological features:
                 - ``'geology'``
                - ``'faults'``
                - ``'fractures'``
        dim : int
            Define whether the geological feature corresponds to a 2D or 3D dataset.
            
        default_fmm_cost : float
            Define the default fast-marching method cost value.
        """
        
        # Initialization
        super().__init__(grid, *args, **kwargs)
        self.feature = feature
        self.dim = dim
        self.default_fmm_cost = default_fmm_cost
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
            self.set_costs(costs, self.default_fmm_cost)
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
        """        
        """
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
        default_name : str, default: 'item {}'
            A format string used to generate default names for items not
            explicitly named in the ``names`` dictionary. The format string
            should include a placeholder (e.g., '{}') that will be replaced by
            the item's index.
        
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
        default_cost: float = None,
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
        default_cost : float, default: 0.5
            The default cost to be applied to items not explicitly listed in
            the ``costs`` dictionary.
        
        Notes
        -----
        This function does not return a value. It rewrites the ``self.costs``
        attribute with the new specified dictionary.
        """
        # Retrieve default fmm cost
        if default_cost is None:
            default_cost = self.default_fmm_cost
        
        # Assign costs
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
        default_model : bool, default: True
            The default value to be applied to items not explicitly listed in
            the ``model`` dictionary.
        
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
        # data = np.empty(self.grid.shape) * np.nan # ISSUES with plotting
        data = np.zeros(self.grid.shape)
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
