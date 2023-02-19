"""
# TODO
# - Documentation
# - Gestion des cartes de densité de probabilités pour le tirage des points
"""

import logging

### External dependencies
import numpy as np

### Typing
from pykasso._typing import Domain, Geology

class PointManager():
    """"""
    def __init__(self, rng, mode:str, domain:Domain, geology:Geology, geologic_ids:list) -> None:
        """"""
        self.rng     = rng
        self.mode    = mode
        self.domain  = domain
        self.geology = geology
        probabilistic_domain = self._get_domain()
        self.probability_map = self._calculate_probability_map(probabilistic_domain)
        self.geologic_ids    = self._controls_geologic_ids_validity(geologic_ids)
    
    def _get_domain(self):
        """"""
        ## TODO - numbers instead of names ?  
        if self.mode == 'uniform':
            data = self.domain.data_volume
        elif self.mode == 'surface_up':
            data = self.domain.faces['up']
        elif self.mode == 'surface_down':
            data = self.domain.faces['down']
        elif self.mode == 'borders':
            data = self.domain.borders['domain']
        elif self.mode == 'side_borders':
            data = self.domain.borders['domain_sides']
        elif self.mode == 'water_level_surface':
            data = self.domain.phreatic['water_level_surface']
        elif self.mode == 'water_level_surface_borders':
            data = self.domain.borders['water_level_surface_sides']
        elif self.mode == 'water_level_borders':
            data = self.domain.phreatic['water_level_borders']     
        return data
    
    def _calculate_probability_map(self, array:np.ndarray) -> list:
        """"""
        i,j,k = np.indices(array.shape)
        i,j,k = i.flatten(), j.flatten(), k.flatten(), 
        state = array.flatten()
        nodes = list(zip(state, i, j, k))
        probability_map = [(i,j,k) for (state,i,j,k) in nodes if state == 1]
        return probability_map
    
    def _generate_coordinates(self, size:int=1) -> np.ndarray:
        """"""
        indices = self.rng.choice(self.probability_map, size=size)
        i, j, k = zip(*indices)
        i, j, k = np.array(i), np.array(j), np.array(k)
        x = self.domain.grid.xmin + (i + self.rng.random()) * self.domain.grid.dx
        y = self.domain.grid.ymin + (j + self.rng.random()) * self.domain.grid.dy
        z = self.domain.grid.zmin + (k + self.rng.random()) * self.domain.grid.dz
        return np.dstack((x,y,z))[0]
        
    def _generate_3D_coordinate_from_2D_coordinate(self, coordinate:tuple) -> tuple:
        """"""
        x, y = coordinate
        i, j = self.domain.grid.get_indices(x, y)
        new_probability_map = [(i_, j_, k_) for (i_, j_, k_) in self.probability_map if ((i_ == i) and (j_ == j))]
        # if len(new_probability_map) == 0 # TODO
        i_, j_, k_ = self.rng.choice(new_probability_map)
        z = self.domain.grid.zmin + (k_ + self.rng.random()) * self.domain.grid.dz
        return (x, y, z)
    
    def _controls_geologic_ids_validity(self, geologic_ids:list) -> list:
        """ 
        Controls validity of 'geologic_ids'.
        """
        if geologic_ids is None:
            return None
        else:
            values = self.geology.stats.index.to_list()
            validated_geology_ids = [geologic_id for geologic_id in geologic_ids if geologic_id in values]
            # TODO - log
            if len(validated_geology_ids) == 0:
                pass
                # TODO - log - erreur
                return None
            return validated_geology_ids
        
    def _is_coordinate_2D_valid(self, coordinate:tuple) -> bool:
        """Checks if the 2D point is valid."""
        x, y = coordinate
        if self.domain.is_coordinate_2D_valid(x,y):
            return True  
        else:
            return False
        
    def _is_coordinate_3D_valid(self, coordinate:tuple) -> bool:
        """Checks if the 3D points is valid."""
        x, y, z = coordinate
            
        # Case 1 - No geology
        if self.geologic_ids is None:
            return self.domain.is_coordinate_in_domain(x,y,z)
        
        # Case 2 - With geology
        else:
            i, j, k = self.domain.grid.get_indices(x,y,z)
            return self.domain.is_coordinate_in_domain(x,y,z) and (self.geology.data_volume[i,j,k] in self.geologic_ids)  