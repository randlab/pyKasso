"""
This module contains classes modeling the geological features.
"""

### Internal dependencies
import sys

### External dependencies
import numpy as np
import pandas as pd

### Local dependencies
from pykasso._utils.readers import DataReader

sks = sys.modules['pykasso.model.sks']

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
    """Class modeling a geological feature."""
    
    def __init__(self, label: str, dim: int, *args, **kwargs) -> None:
        """
        Constructs a geological feature.

        Parameters
        ----------
        label : str
            _description_
        dim : int
            _description_
        """
        super().__init__(*args, **kwargs)
        
        self.label = label
        self.dim = dim
        self.data_surface = None
        self.data_volume = None
        
    def set_data(self, data, grid, axis: str = None) -> None:
        
        if isinstance(data, np.ndarray):
            if self.dim == 2:
                self.data_surface = data
            elif self.dim == 3:
                self.data_volume = data
        else:
            if self.dim == 2:
                self.data_surface = self._set_data_from_file(grid,
                                                             data,
                                                             False)
            elif self.dim == 3:
                self.data_volume = self._set_data_from_file(grid,
                                                            data,
                                                            True,
                                                            axis)
        return None
    
    def to_file(self, name: str = None, data: str = ''):
        
        grid = sks.ACTIVE_PROJECT['settings']['grid']

        
        path = sks.ACTIVE_PROJECT['simulation_locations'][-1]
        
        if name is not None:
            extension = name.split('.')[-1]
            # if extension == 'a':
            #     pass
            # elif extension == 'a':
            #     pass
            # elif extension == 'a':
            #     pass
            # elif extension == 'a':
            #     pass
            # elif extension == 'a':
            #     pass
        else:
            extension = ''
            filename = self.label + extension
        return None
    
    def _save_to_gslib(self):
        grid = sks.ACTIVE_PROJECT['settings']['grid']
        x0, y0, z0 = grid['x0'], grid['y0'], grid['z0']
        nx, ny, nz = grid['nx'], grid['ny'], grid['nz']
        dx, dy, dz = grid['dx'], grid['dy'], grid['dz']
        comment = ("x0: {}, y0: {}, z0:, {}, nx: {}, ny: {}, nz: {}, dx: {},"
                   "dy: {}, dz: {}").format(x0, y0, z0, nx, ny, nz, dx, dy, dz)
        n_var = 1
        
        
        pass
    
    def _save_to_np_pickle(self) -> None:
        # np.save()
        return None
    
    def _save_to_csv(self) -> None:
        return None
    
    def _save_to_asc(self) -> None:
        return None
    
    def _save_to_grd(self) -> None:
        return None
                
#######################
### Main subclasses ###
#######################


class Surface(GeologicFeature):
    """Subclass modeling a two dimensional geological feature."""
    
    def __init__(self, label, *args, **kwargs):
        dim = 2
        super().__init__(label, dim, *args, **kwargs)
    
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

   
class Volume(GeologicFeature):
    """Subclass modeling a three dimensional geological feature."""
    
    def __init__(self, label, *args, **kwargs):
        dim = 3
        super().__init__(label, dim, *args, **kwargs)
    
    def _set_costs(self, costs: dict) -> None:
        self.costs = costs
        return None
        
    def _compute_statistics(self, grid):
        """Computes the statistics (counts and frequency) on the data."""
        values, counts = np.unique(self.data_volume, return_counts=True)
        stats = {
            'counts': counts,
            'freq': counts / grid.nodes,
            'volume': counts * grid.node_volume,
        }
        self.stats = pd.DataFrame(data=stats, index=values)
        return None
    
        
#####################
### 3D Subclasses ###
#####################


class Geology(Volume):
    """Class modeling the geologic model."""
    
    def __init__(self, *args, **kwargs):
        label = 'geology'
        super().__init__(label, *args, **kwargs)


class Faults(Volume):
    """Class modeling the faults model."""
    
    def __init__(self, *args, **kwargs):
        label = 'faults'
        super().__init__(label, *args, **kwargs)
