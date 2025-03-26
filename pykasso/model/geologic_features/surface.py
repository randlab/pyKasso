import numpy as np
import pandas as pd

from .geologicfeature import GeologicFeature
from ...core.grid import Grid


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
