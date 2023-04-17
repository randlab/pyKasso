"""
This module contains a point manager class designed for complex conditionnal
point generation.
"""

### Internal dependencies
import logging

### External dependencies
import numpy as np
from shapely.geometry import Point

### Typing
from pykasso._typing import Domain, Geology, RandomNumberGenerator


class PointGenerator():
    """
    Class modeling a point manager tool for point generation and management.
    The tool is not stored in the memory and is erased after usage.
    """
    
    def __init__(self, rng: RandomNumberGenerator, mode: str, domain: Domain,
                 geology: Geology, geologic_ids: list) -> None:
        """_summary_

        Parameters
        ----------
        rng : RNG
            A random number generator from numpy :
            https://numpy.org/devdocs/reference/random/generator.html
        mode : str
            _description_
        domain : Domain
            _description_
        geology : Geology
            _description_
        geologic_ids : list
            _description_
        """
        
        # TODO
        logging.getLogger('points.')
        
        # Initialization
        self.rng = rng
        self.mode = mode
        self.domain = domain
        self.geology = geology
        
        # Controls dictionary validity of geological constraints
        self.geologic_ids = self._controls_geologic_ids_validity(geologic_ids)
        
        #
        self.valid_domain_volume = self._get_domain()
        
        #
        test = np.sum(self.valid_domain_volume, axis=2) > 0
        self.valid_domain_surface = test
        
        #
        self.probability_map = self._calculate_probability_map(
            self.valid_domain_volume
        )
        
    def _controls_geologic_ids_validity(self, geologic_ids: list) -> list:
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
                    msg = ("Declared geologic id #{} is not present in "
                           "geology model.".format(geologic_id))
                    logging.warning(msg)
                
            if len(validated_geology_ids) == 0:
                msg = ("No declared geologic id is present in the geologic "
                       "model, geologic constraints are ignored.")
                logging.warning(msg)
                return None
            else:
                return validated_geology_ids

    def _get_domain(self):
        """"""
        # try:
        probabilistic_domain = self.domain.get_subdomain(self.mode)
        # except:
        # TODO
        # print('points.py - _get_domain - ERROR')
        # raise
        
        # TODO
        if self.geologic_ids is not None:
            domain_geology = np.zeros_like(probabilistic_domain)

            for geologic_id in self.geologic_ids:
                domain_geology = np.where(
                    self.geology.data_volume == geologic_id,
                    1,
                    domain_geology
                )
            test = np.logical_and(probabilistic_domain, domain_geology)
            probabilistic_domain_geology = test
            
            # if not probabilistic_domain_geology.all():
            #     # TODO
            #     print('ERROR')
            # else:
            probabilistic_domain = probabilistic_domain_geology
        
        return probabilistic_domain
    
    def _calculate_probability_map(self, array: np.ndarray) -> list:
        """"""
        i, j, k = np.indices(array.shape)
        i, j, k = i.flatten(), j.flatten(), k.flatten()
        state = array.flatten()
        nodes = list(zip(state, i, j, k))
        probability_map = [(i, j, k) for (state, i, j, k)
                           in nodes if state == 1]
        return probability_map
    
    def _generate_coordinates(self, size: int = 1) -> np.ndarray:
        """"""
        indices = self.rng.choice(self.probability_map, size=size)
        i, j, k = zip(*indices)
        i, j, k = np.array(i), np.array(j), np.array(k)
        x = (self.domain.grid.xmin + (i + self.rng.random())
             * self.domain.grid.dx)
        y = (self.domain.grid.ymin + (j + self.rng.random())
             * self.domain.grid.dy)
        z = (self.domain.grid.zmin + (k + self.rng.random())
             * self.domain.grid.dz)
        return np.dstack((x, y, z))[0]
        
    def _generate_3D_coord_from_2D_coord(self, coordinate: tuple) -> tuple:
        """"""
        x, y = coordinate
        i, j = self.domain.grid.get_indices(x, y)
        new_probability_map = ([(i_, j_, k_) for (i_, j_, k_)
                                in self.probability_map
                                if ((i_ == i) and (j_ == j))])
        # if len(new_probability_map) == 0 # TODO = erreur
        i_, j_, k_ = self.rng.choice(new_probability_map)
        z = (self.domain.grid.zmin + (k_ + self.rng.random())
             * self.domain.grid.dz)
        return (x, y, z)
    
    def _is_coordinate_2D_valid(self, coordinate: tuple) -> bool:
        """Checks if the 2D point is valid."""
        x, y = coordinate
        if self.domain.is_coordinate_2D_valid(x, y):
            i, j = self.domain.grid.get_indices(x, y)
            out = self.valid_domain_surface[i, j]
            return out
        else:
            return False
        
    def _is_coordinate_3D_valid(self, coordinate: tuple) -> bool:
        """Checks if the 3D points is valid."""
        x, y, z = coordinate
            
        # Case 1 - No geology
        if self.geologic_ids is None:
            return self.domain.is_coordinate_in_domain(x, y, z)
        
        # Case 2 - With geology
        else:
            i, j, k = self.domain.grid.get_indices(x, y, z)
            test = self.domain.is_coordinate_in_domain(x, y, z)
            list_ = self.geology.data_volume[i, j, k] in self.geologic_ids
            out = (test and list_)
            return out
