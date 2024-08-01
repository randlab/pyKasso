"""
This module contains a point manager class designed for complex conditionnal
point generation.
"""

### Internal dependencies
import logging

### External dependencies
import numpy as np
import pandas as pd

### Typing
from pykasso._typing import Domain, Geology, RandomNumberGenerator


class PointGenerator():
    """
    Class modeling a point manager tool for point generation and management.
    The tool is not stored in the memory and is erased after usage.
    """
    
    def __init__(self, rng: RandomNumberGenerator, domain: Domain,
                 geology: Geology, subdomain: str, geologic_ids: list) -> None:
        """_summary_

        Parameters
        ----------
        rng : RNG
            A random number generator from numpy :
            https://numpy.org/devdocs/reference/random/generator.html
        subdomain : str
            _description_
        domain : Domain
            _description_
        geology : Geology
            _description_
        geologic_ids : list
            _description_
        """
        # Todo
        # logging.getLogger('points.')
        
        ### Initialization
        self.rng = rng
        self.grid = domain.grid
        
        ### Compute domain validity
        # Retrieve subdomain
        self.data_volume = domain.get_subdomain(subdomain)
        
        # Cross with geological constraints
        if geologic_ids is not None:
            geologic_ids = self._controls_geologic_ids(geology, geologic_ids)
            if len(geologic_ids) > 0:
                dom_geol = np.zeros_like(self.data_volume)
                dom_geol = np.where(
                    np.isin(geology.data_volume, geologic_ids),
                    1,
                    dom_geol
                )
                self.data_volume = np.logical_and(self.data_volume, dom_geol)
        
        # Calculate domain surface validity
        self.data_surface = (np.sum(self.data_volume, axis=2) > 0)
        
        # Retrieve the indices where point generation is allowed
        self.valid_cells = self._get_valid_cells(self.data_volume)
        
    def _controls_geologic_ids(self, geology, geologic_ids: list) -> list:
        """
        Controls validity of 'geologic_ids'.
        """
        values = geology.stats.index.to_list()
        validated_geology_ids = []
        
        for geologic_id in geologic_ids:
            if geologic_id in values:
                validated_geology_ids.append(geologic_id)
            else:
                msg = ("Declared geologic id #{} is not present in "
                       "geology model.".format(geologic_id))
                logging.warning(msg)
            
        if len(validated_geology_ids) == 0:
            msg = ("None of the geologic ids declared are present in the "
                   "geologic model, geologic constraints are ignored.")
            logging.warning(msg)
            return []
        else:
            return validated_geology_ids
    
    def _get_valid_cells(self, array: np.ndarray) -> list:
        """
        TODO
        """
        i, j, k = np.indices(array.shape)
        i, j, k = i.flatten(), j.flatten(), k.flatten()
        test = array.flatten()
        i = i[test == 1]
        j = j[test == 1]
        k = k[test == 1]
        
        valid_cells = pd.DataFrame({
            'i': i,
            'j': j,
            'k': k,
        })
        # valid_cells = list(zip(i, j, k))
        # nodes = list(zip(test, i, j, k))
        # valid_cells = [(i, j, k) for (test, i, j, k) in nodes if test == 1]
        return valid_cells
    
    def _generate_points(self, size: int = 1) -> np.ndarray:
        """
        TODO
        """
        indices = self.rng.choice(self.valid_cells, size=size)
        i, j, k = zip(*indices)
        i, j, k = np.array(i), np.array(j), np.array(k)
        x = (self.grid.xmin + (i + self.rng.random()) * self.grid.dx)
        y = (self.grid.ymin + (j + self.rng.random()) * self.grid.dy)
        z = (self.grid.zmin + (k + self.rng.random()) * self.grid.dz)
        return np.dstack((x, y, z))[0]
        
    def _3D_point_from_2D_point(self, point: tuple) -> tuple:
        """
        TODO
        """
        x, y = point
        i, j = self.grid.get_indices(point)
        i, j = int(i), int(j) 
        new_valid_cells = self.valid_cells[(self.valid_cells['i'] == i)
                                           & (self.valid_cells['j'] == j)]
        i_, j_, k_ = self.rng.choice(new_valid_cells)
        z = (self.grid.zmin + (k_ + self.rng.random()) * self.grid.dz)
        return (x, y, z)
    
    def _is_point_valid(self, point: tuple) -> bool:
        """
        Check if 2D or 3D point is valid.
        """
        if self.grid.is_inbox(point):
            if len(point) == 2:
                i, j = self.grid.get_indices(point)
                out = self.data_surface[i, j]
            elif len(point) == 3:
                i, j, k = self.grid.get_indices(point)
                out = self.data_volume[i, j, k]
            return out
        else:
            return False
