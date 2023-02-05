"""
Module modeling the domain.
"""

### Local dependencies
from .geologic_features import Surface

### External dependencies
import numpy as np
from scipy.ndimage import binary_dilation
from matplotlib.path import Path


class Domain():
    """
    Class modeling the extension of the domain within the grid.
    """
    def __init__(self, grid, data=None, delimitation=None, topography=None, bedrock=None) -> None:
        """ 
        TODO
        """
        # Initialization
        self.grid        = grid
        self.data_volume = np.zeros_like(grid.data_volume)
        
        if delimitation is not None: self.delimitation = delimitation
        if topography   is not None: self.topography   = topography
        if bedrock      is not None: self.bedrock      = bedrock
        
        ### Computes domain extension
        
        # Volumetric domain extension
        if data is None:
            if delimitation is not None: self.data_volume += delimitation.data_volume
            if topography   is not None: self.data_volume += topography.data_volume
            if bedrock      is not None: self.data_volume += np.logical_not(bedrock.data_volume).astype(int)
            self.data_volume = np.where(self.data_volume == self.data_volume.max(), 1, 0)
        else:
            # TODO - mieux définir le domaine, ai-je le droit de le couper avec bedrock par exemple ?
            self.data_volume = data
        
        # Surfacic domain extension
        self.data_surface = self.data_volume.max(axis=2)
        
        # Computes domain border
        k = np.zeros((3,3,3), dtype=int); k[:,1,1], k[1,:,1], k[1,1,:] = 1,1,1
        self.data_volume_border = binary_dilation(self.data_volume==0, k, border_value=1) & self.data_volume
        
        # Computes domain surfaces
        # TODO - zu verbessern
        i,j,k = np.indices(self.data_volume.shape)
        test = self.data_volume.astype('bool')
        self.surface_up   = k.max(axis=2, initial=-1, where=test)
        self.surface_down = k.min(axis=2, initial=0 , where=test)
    
    def is_coordinate_in_domain(self, x:float, y:float, z:float) -> bool:
        if self.grid.is_inbox(x,y,z):
            i, j, k = self.grid.get_indices(x,y,z)
            return bool(self.data_volume[i,j,k])
        else:
            return False
        
    def is_coordinate_on_border(self, x:float, y:float, z:float) -> bool:
        if self.grid.is_inbox(x,y,z):
            i, j, k = self.grid.get_indices(x,y,z)
            return bool(self.data_volume_border[i,j,k])
        else: 
            return False
        
    def is_coordinate_2D_valid(self, x:float, y:float) -> bool:
        if self.grid.path.contains_point((x,y)):
            return bool(self.data_surface[x,y])
        else:
            return False
        
    def get_elevation_from_surface(self, x:float, y:float, surface:str) -> float:
        if surface == 'up':
            return self._get_elevation_from_upper_surface(x, y)
        elif surface == 'down':
            return self._get_elevation_from_lower_surface(x, y)
        else:
            return None
    
    def _get_elevation_from_upper_surface(self, x:float, y:float) -> float:
        i, j = self.grid.get_indices(x,y)
        k = self.surface_up[i,j]
        return self.grid.Z[i,j,k] + self.grid.dz/2 - self.grid.dz/100

    def _get_elevation_from_lower_surface(self, x:float, y:float) -> float:
        i, j = self.grid.get_indices(x,y)
        k = self.surface_down[i,j]
        return self.grid.Z[i,j,k] - self.grid.dz/2 + self.grid.dz/100
    
    
class Delimitation():
    """
    Class modeling the catchment delimitation of the studied domain.
    """
    def __init__(self, vertices:list, grid) -> None:
        """
        Creates a mask.
        This class is designed to describe a particular catchment delimitation.

        Parameters
        ----------
        data : list | str
            List of vertices from Python list or datafile location.
        name : str, optional
            Name of the mask.
        """
        self.vertices = vertices
        self._set_path()
        self._set_mask(grid)
            
    def __str__(self) -> str:
        l0 = "Mask \n"
        l1 = "{} vertices : \n{}".format(len(self.vertices), self.vertices)
        return l0+l1
         
    def _set_path(self) -> None:
        """ 
        Sets the polygon with a matplotlib Path
        """
        path_vertices = self.vertices.copy()
        path_vertices.append(path_vertices[0])
        self.path = Path(path_vertices)
        return None
    
    def _set_mask(self, grid) -> None:
        """ 
        Sets the mask array with numpy-array
        """
        row, col = np.indices((grid.nx, grid.ny))
        pts = np.column_stack((grid.X[row, col, 0].flatten(), grid.Y[row, col, 0].flatten()))
        msk = self.path.contains_points(pts).reshape((grid.nx, grid.ny)).astype(int)
        self.data_volume = np.repeat(msk[:, :, np.newaxis], grid.nz, axis=2)
        return None


class Topography(Surface):
    """
    Class modeling the topography of the studied domain.
    """
    def __init__(self, data, grid, **kwargs):
        """
        TODO
        """
        label = 'topography'
        super().__init__(label, data, grid, **kwargs)
        self._surface_to_volume('<=', grid)
        
        
class Bedrock(Surface):
    """
    Class modeling the bedrock elevation of the studied domain.
    """
    def __init__(self, data, grid, **kwargs):
        label = 'bedrock_elevation'
        super().__init__(label, data, grid, **kwargs)
        self._surface_to_volume('<=', grid)