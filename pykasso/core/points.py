"""
# TODO
# - Documentation
# - Gestion des cartes de densité de probabilités pour le tirage des points
"""

import sys
import logging

### Typing
from pykasso._typing import Domain, Geology

### External dependencies
import numpy as np

# sks  = sys.modules['pykasso.core.sks']
# import pykasso.visualization as pkv
# pkv.show_array(probabilistic_domain_geology, ghost=1)

class PointManager():
    """"""
    def __init__(self, rng, mode:str, domain:Domain, geology:Geology, geologic_ids:list) -> None:
        """"""
        logging.getLogger('points.')
        self.rng     = rng
        self.mode    = mode
        self.domain  = domain
        self.geology = geology
        self.geologic_ids = self._controls_geologic_ids_validity(geologic_ids)
        self.valid_domain_volume  = self._get_domain()
        self.valid_domain_surface = np.sum(self.valid_domain_volume, axis=2) > 0 
        self.probability_map = self._calculate_probability_map(self.valid_domain_volume)
        
    def _controls_geologic_ids_validity(self, geologic_ids:list) -> list:
        """ 
        Controls validity of 'geologic_ids'.
        """
        if geologic_ids is None:
            return None
        else:
            values = self.geology.stats.index.to_list()
            validated_geology_ids = []
            
            for geologic_id in geologic_ids:
                if geologic_id in values:
                    validated_geology_ids.append(geologic_id)
                else:
                    logging.warning("Declared geologic id #{} is not present in geology model.".format(geologic_id))
                
            if len(validated_geology_ids) == 0:
                logging.warning("No declared geologic id is present in the geologic model, geologic constraints are ignored.".format(geologic_id))
                return None
            else:
                return validated_geology_ids

    def _get_domain(self):
        """"""
        modes = {
            # domain
            "domain"         : "self.domain.data_volume",
            "domain_surface" : "self.domain.faces['up']",
            "domain_bottom"  : "self.domain.faces['down']",
            # borders
            "domain_borders"         : "self.domain.borders['domain_full']",
            "domain_borders_sides"   : "self.domain.borders['domain_sides']",
            "domain_borders_surface" : "self.domain.borders['domain_up']",
            "domain_borders_bottom"  : "self.domain.borders['domain_down']", 
            # phreatic
            "vadose"                   : "self.domain.phreatic['vadose_zone']",
            "vadose_borders"           : "self.domain.borders['vadose_zone']",
            "phreatic"                 : "self.domain.phreatic['phreatic_zone']",
            "phreatic_surface"         : "self.domain.phreatic['water_level_surface']",
            "phreatic_borders_surface" : "self.domain.borders['water_level_surface']",
        }  
        try:
            probabilistic_domain = eval(modes[self.mode])
        except:
            # TODO
            print('points.py - _get_domain - ERROR')
            raise
        
        # TODO
        if self.geologic_ids is not None:
            domain_geology = np.zeros_like(probabilistic_domain)

            for geologic_id in self.geologic_ids:
                domain_geology = np.where(self.geology.data_volume == geologic_id, 1, domain_geology)
            
            probabilistic_domain_geology = np.logical_and(probabilistic_domain, domain_geology)
            
            # if not probabilistic_domain_geology.all():
            #     # TODO
            #     print('ERROR')
            # else:
            probabilistic_domain = probabilistic_domain_geology
        
        return probabilistic_domain
    
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
    
    def _is_coordinate_2D_valid(self, coordinate:tuple) -> bool:
        """Checks if the 2D point is valid."""
        x, y = coordinate
        if self.domain.is_coordinate_2D_valid(x,y):
            i, j = self.domain.grid.get_indices(x, y)
            return self.valid_domain_surface[i,j]
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