"""
Module to model the domain.
"""
### Local dependencies
from .geologic_features import Surface

### External dependencies
import numpy as np
import numpy.ma as ma
from scipy.ndimage import binary_dilation
from matplotlib.path import Path

### Typing
from typing import Union
from pykasso.typing import Grid, Delimitation, Topography, Bedrock, WaterLevel


class Domain():
    """
    Class modeling the extension of the domain within the grid.
    """
    def __init__(self, grid:Grid, data:np.ndarray=None, delimitation:Delimitation=None, topography:Topography=None, bedrock:Bedrock=None, water_level:WaterLevel=None) -> None:
        """ """
        # Initialization
        self.grid        = grid
        self.data_volume = np.zeros_like(grid.data_volume)
        
        if delimitation is not None: self.delimitation = delimitation
        if topography   is not None: self.topography   = topography
        if bedrock      is not None: self.bedrock      = bedrock
        if water_level  is not None: self.water_level  = water_level
        
        ### Computes domain extension
        # Volumetric domain extension
        if data is None:
            if delimitation is not None: self.data_volume += delimitation.data_volume
            if topography   is not None: self.data_volume += topography.data_volume
            if bedrock      is not None: self.data_volume += np.logical_not(bedrock.data_volume).astype(int)
            self.data_volume = np.where(self.data_volume == self.data_volume.max(), 1, 0)
        else:
            # TODO - mieux définir le domaine, ai-je le droit de le couper avec bedrock par exemple ?
            # Montrer des exemples
            self.data_volume = data
        # Apparent surface from domain according to dimension
        self.data_surfaces = {
            'x' : self.data_volume.max(axis=0),
            'y' : self.data_volume.max(axis=1),
            'z' : self.data_volume.max(axis=2),
        } 
        
        ### Computes borders from domain
        self.data_borders = {}
        # Full borders
        k = np.zeros((3,3,3), dtype=int); k[:,1,1], k[1,:,1], k[1,1,:] = 1,1,1
        self.data_borders['xyz'] = binary_dilation(self.data_volume==0, k, border_value=1) & self.data_volume
        # x and y directions borders
        k = np.zeros((3,3), dtype=int); k[:,1], k[1,:] = 1,1
        self.data_borders['xy'] = np.zeros_like(self.data_volume)
        for z in range(self.grid.nz):            
            self.data_borders['xy'][:,:,z] = binary_dilation(self.data_volume[:,:,z]==0, k, border_value=1) & self.data_volume[:,:,z]
        
        ### Computes specific surfaces from domain
        self.faces = {}
        I, J, K = np.indices(self.data_volume.shape)
        domain = self.data_volume.astype('bool')
        face_up   = K.max(axis=2, initial=-1, where=domain)
        face_down = K.min(axis=2, initial=self.grid.nz, where=domain)
        i_, j_ = np.indices(self.data_surfaces['z'].shape)
        i_, j_ = i_.flatten(), j_.flatten()
        k_up, k_down = face_up.flatten(), face_down.flatten()
        nodes_up, nodes_down = list(zip(i_,j_,k_up)), list(zip(i_,j_,k_down))
        valid_nodes_up   = [(i,j,k) for (i,j,k) in nodes_up   if k in list(range(self.grid.nz))]
        valid_nodes_down = [(i,j,k) for (i,j,k) in nodes_down if k in list(range(self.grid.nz))]
                                  
        i,j,k = zip(*valid_nodes_up)
        volume_up = np.zeros_like(self.data_volume) 
        volume_up[i,j,k] = 1
        self.faces['up'] = volume_up
        
        i,j,k = zip(*valid_nodes_down)
        volume_down = np.zeros_like(self.data_volume) 
        volume_down[i,j,k] = 1
        self.faces['down'] = volume_down
        
        if self.water_level is not None:
            self._set_water_table()
        
    # TODO - names ?
    def _set_water_table(self):
        
        water_volume  = self.water_level.data_volume
        water_surface = self.water_level.data_surface
        
        self.phreatic = {
            'vadose_zone'         : np.logical_and(self.data_volume, np.logical_not(water_volume)),
            'phreatic_zone'       : np.logical_and(self.data_volume, water_volume),
            'water_level'         : np.logical_and(self.data_volume, water_surface),
            'water_level_borders' : np.logical_and(self.data_borders['xy'], water_volume),
        }
        
        # I, J, K = np.indices(self.data_volume.shape)
        # domain = self.phreatic['vadose_zone'].astype('bool')
        # face = K.max(axis=2, initial=-1, where=domain)
        # i_, j_ = np.indices(self.data_surfaces['z'].shape)
        # i_, j_ = i_.flatten(), j_.flatten()
        # k_ = face.flatten()
        # nodes = list(zip(i_,j_,k_))
        # valid_nodes = [(i,j,k) for (i,j,k) in nodes if k in list(range(self.grid.nz))]
        
        # i,j,k = zip(*valid_nodes)
        # volume = np.zeros_like(self.data_volume) 
        # volume[i,j,k] = 1
        # self.phreatic['water_level'] = volume
        
        return None
        
        
        
    # def is_coordinate_in_domain(self, x:float, y:float, z:float) -> bool:
    #     if self.grid.is_inbox(x,y,z):
    #         i, j, k = self.grid.get_indices(x,y,z)
    #         return bool(self.data_volume[i,j,k])
    #     else:
    #         return False
        
    # def is_coordinate_on_border(self, x:float, y:float, z:float) -> bool:
    #     if self.grid.is_inbox(x,y,z):
    #         i, j, k = self.grid.get_indices(x,y,z)
    #         return bool(self.data_volume_border[i,j,k])
    #     else: 
    #         return False
        
    # def is_coordinate_2D_valid(self, x:float, y:float) -> bool:
    #     if self.grid.path.contains_point((x,y)):
    #         return bool(self.data_surface[x,y])
    #     else:
    #         return False
        
    # def get_elevation_from_xy_coordinate(self, x:float, y:float, surface:str='up') -> float:
    #     if surface == 'up':
    #         return self._get_elevation_from_surface_up(x, y)
    #     elif surface == 'down':
    #         return self._get_elevation_from_surface_down(x, y)
    
    # def _get_elevation_from_surface_up(self, x:float, y:float) -> float:
    #     i, j = self.grid.get_indices(x,y)
    #     k = self.surfaces['up'][i,j]
    #     return self.grid.Z[i,j,k] + self.grid.dz/2 - self.grid.dz/100

    # def _get_elevation_from_surface_down(self, x:float, y:float) -> float:
    #     i, j = self.grid.get_indices(x,y)
    #     k = self.surfaces['down'][i,j]
    #     return self.grid.Z[i,j,k] - self.grid.dz/2 + self.grid.dz/100
    
    
class Delimitation():
    """
    Class modeling the catchment delimitation of the studied domain.
    """
    def __init__(self, vertices:list, grid:Grid) -> None:
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
    
    def _set_mask(self, grid:Grid) -> None:
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
    def __init__(self, data:Union[str, np.ndarray], grid:Grid, **kwargs) -> None:
        label = 'topography'
        super().__init__(label, data, grid, **kwargs)
        self.data_volume = self._surface_to_volume('<=', grid)
        
        
class Bedrock(Surface):
    """
    Class modeling the bedrock elevation of the studied domain.
    """
    def __init__(self, data:Union[str, np.ndarray], grid:Grid, **kwargs) -> None:
        label = 'bedrock_elevation'
        super().__init__(label, data, grid, **kwargs)
        self.data_volume = self._surface_to_volume('<=', grid)
        
        
class WaterLevel(Surface):
    """
    TODO
    """
    def __init__(self, data:Union[str, np.ndarray], grid:Grid, **kwargs) -> None:
        label = 'water_level'
        super().__init__(label, data, grid, **kwargs)
        self.data_volume  = self._surface_to_volume('<=', grid)
        self.data_surface = self._surface_to_volume('=', grid)