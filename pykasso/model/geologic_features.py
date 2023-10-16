"""
This module contains classes modeling the geological features.
"""

### External dependencies
import numpy as np
import pandas as pd

### Local dependencies
from pykasso._utils.datareader import DataReader
from pykasso.core._namespaces import (DEFAULT_FMM_COSTS,
                                      DEFAULT_FEATURE_PARAMETERS)

### Typing
from typing import Union

##################
### Main class ###
##################

####################### TODO ##########################
# # Processes the image
# if self.label in ['faults', 'fractures']:
#     npy_image = (npy_image[:, :] == 0) * 1
# else:
#     npy_image = npy_image + 1000
#     n_colors = np.unique(npy_image)
#     for i, color in enumerate(np.flip(n_colors)):
#         npy_image = np.where(npy_image == color, i + 1, npy_image)
####################### TODO ##########################


class GeologicFeature(DataReader):
    """
    Class modeling a geological feature.
    """
    
    def __init__(self,
                 feature: str,
                 dim: int,
                 *args,
                 **kwargs,
                 ) -> None:
        """
        Construct a geological feature.

        Parameters
        ----------
        feature : str
            _description_
        dim : int
            _description_
        """
        super().__init__(*args, **kwargs)
        
        self.feature = feature
        self.dim = dim
        self.data_surface = None
        self.data_volume = None
        self.stats = None
        self.names = None
        self.costs = None
        self.model = None
        
        self.set_data(kwargs.get('data'), kwargs.get('axis'))
        self.compute_statistics()
        if self.feature not in ['topography', 'water_level', 'bedrock']:
            self.set_names(kwargs.get('names'))
            self.set_costs(kwargs.get('costs'))
            self.set_model(kwargs.get('model'))
        
    def set_data(self,
                 data: Union[None, str, np.ndarray],
                 axis: str = None,
                 ) -> None:
        """TODO"""
        # If no data is provdided
        if data is None:
            
            value = DEFAULT_FEATURE_PARAMETERS[self.feature]['nodata']

            if self.dim == 2:
                self.data_surface = self._get_data_full_2D(value)
            elif self.dim == 3:
                self.data_volume = self._get_data_full_3D(value)

            return None
        
        # Else
        if isinstance(data, np.ndarray):
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
    
    def compute_statistics(self) -> None:
        """
        Compute the statistics (counts and frequency) on the data.
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
    
    def set_names(self, names: dict[int, str]) -> None:
        """TODO"""
        ids = self.stats.index
        names_df = {}
        default_name = DEFAULT_FEATURE_PARAMETERS[self.feature]['name']
        for id in ids:
            names_df[id] = names.get(id, default_name.format(id))
        self.names = names_df
        return None
        
    def set_costs(self, costs: dict[int, str]) -> None:
        """TODO"""
        ids = self.stats.index
        costs_df = {}
        default_cost = DEFAULT_FMM_COSTS[self.feature]
        for id in ids:
            costs_df[id] = costs.get(id, default_cost)
        self.costs = costs_df
        return None
    
    def set_model(self, model: dict[int, str]) -> None:
        """TODO"""
        ids = self.stats.index
        if self.feature in ['faults', 'fractures']:
            model.setdefault(0, False)
        default_model = DEFAULT_FEATURE_PARAMETERS[self.feature]['model']
        model_df = {}
        for id in ids:
            model_df[id] = model.get(id, default_model)
        self.model = model_df
        return None
    
    def get_geologic_units(self, units: list[int]) -> np.ndarray:
        """TODO"""
        data = np.zeros(self.grid.shape)
        data = np.where(np.isin(self.data_volume, units), self.data_volume, 0)
        return data
    
    def get_geologic_model(self) -> np.ndarray:
        """TODO"""
        valid_ids = [id_ for (id_, boolean) in self.model.items() if boolean]
        geology = self.get_geologic_units(valid_ids)
        return geology
    
    @property
    def overview(self) -> pd.DataFrame:
        """TODO"""
        index = self.stats.index
        names = self.names.values()
        costs = self.costs.values()
        model = self.model.values()
        data = {
            'names': names,
            'costs': costs,
            'model': model,
        }
        df = pd.DataFrame(data, index=index)
        df = pd.merge(df, self.stats, left_index=True, right_index=True)
        return df

#######################
### Main subclasses ###
#######################


class Surface(GeologicFeature):
    """
    Subclass modeling a two dimensional geological feature.
    """
    
    def __init__(self, feature, *args, **kwargs):
        dim = 2
        super().__init__(feature, dim, *args, **kwargs)
    
    def _surface_to_volume(self, condition, grid) -> None:
        """Converts a two dimensional array in a three dimensional array."""
        k = grid.get_k(self.data_surface)
        data_volume = np.zeros((grid.nx, grid.ny, grid.nz))
        for z in range(grid.nz):
            data_volume[:, :, z] = z
            if condition == '>=':
                test = data_volume[:, :, z] >= k
                data_volume[:, :, z] = np.where(test, 1, 0)
            elif condition == '=':
                test = data_volume[:, :, z] == k
                data_volume[:, :, z] = np.where(test, 1, 0)
            elif condition == '<=':
                test = data_volume[:, :, z] <= k
                data_volume[:, :, z] = np.where(test, 1, 0)
        self.data_volume = data_volume
        return None
    
    def compute_statistics(self) -> None:
        """
        Compute the statistics (counts and frequency) on the data.
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
    
    def __init__(self,
                 *args,
                 **kwargs,
                 ) -> None:
        feature = 'geology'
        dim = 3
        super().__init__(feature, dim, *args, **kwargs)


class Faults(GeologicFeature):
    """
    Class modeling the faults model.
    """
    
    def __init__(self,
                 *args,
                 **kwargs,
                 ) -> None:
        feature = 'faults'
        dim = 3
        super().__init__(feature, dim, *args, **kwargs)
