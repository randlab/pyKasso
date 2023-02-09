"""
# TODO
# - Documentation
# - Gestion des cartes de densité de probabilités pour le tirage des points
"""

import logging

### External dependencies
import numpy as np

### Typing
from pykasso.typing import Domain, Geology

class PointManager():
    """"""
    def __init__(self, rng, settings:dict, domain:Domain, geology:Geology, piezometry,geologic_ids:list) -> None:
        """"""
        self.rng        = rng
        self.settings   = settings
        self.domain     = domain
        self.geology    = geology
        self.piezometry = piezometry
        self.probability_domain = self._evaluate_probability_domain()
        self.probability_map = self._calculate_probability_map()
        self.geologic_ids    = self._controls_geologic_ids_validity(geologic_ids)
    
    def _evaluate_probability_domain(self):
        """"""
        if self.settings['mode']   == 'surface_up':
            data = self.domain.surfaces['up']
        elif self.settings['mode'] == 'surface_down':
            data = self.domain.surfaces['down']
        elif self.settings['mode'] == 'uniform':
            data = self.domain.data_volume
        elif self.settings['mode'] == 'xyz_borders':
            data = self.domain.data_borders['xyz']
        elif self.settings['mode'] == 'xy_borders':
            data = self.domain.data_borders['xy']
            
        if self.settings['use_piezometry'] == True :
            if self.piezometry is not None:
                data = data & self.piezometry.data_volume
                
        return None
    
    def _calculate_probability_map(self, array:np.ndarray) -> list:
        """"""
        i,j,k = np.indices(array.shape)
        i,j,k = i.flatten(), j.flatten(), k.flatten(), 
        state = array.flatten()
        nodes = list(zip(*state, (i,j,k)))
        probability_map = [(i,j,k) for (state, (i,j,k)) in nodes if state == 1]
        return probability_map
    
    def _generate_coordinate(self, probability_map:list) -> tuple:
        """"""
        i,j,k = self.rng.choice(probability_map)
        x = self.domain.grid.xmin + (i + self.rng.random()) * self.domain.grid.dx
        y = self.domain.grid.ymin + (j + self.rng.random()) * self.domain.grid.dy
        z = self.domain.grid.zmin + (k + self.rng.random()) * self.domain.grid.dz
        return (x,y,z)
        
    # def _evaluates_modes(self) -> dict:
    #     """"""
    #     # Defines environement for 'eval' function
    #     env = {
    #         'rng'    : self.rng,
    #         'grid'   : self.domain.grid,
    #         'domain' : self.domain,
    #     }
        
    #     # x case
    #     if self.modes['x'] == 'uniform':
    #         x_lambda_function = 'lambda : grid.xmin + grid.nx * rng.random() * grid.dx'
    #         x_function = eval(x_lambda_function, env)
    #     else:
    #         x_function = eval(x_lambda_function, env)
        
    #     # y case
    #     if self.modes['y'] == 'uniform':
    #         y_lambda_function = 'lambda : grid.ymin + grid.ny * rng.random() * grid.dy'
    #         y_function = eval(y_lambda_function, env)
    #     else:
    #         y_function = eval(y_lambda_function, env)
        
    #     # z case
    #     if self.modes['z'] == 'uniform':
    #         z_lambda_function = 'lambda : grid.zmin + grid.nz * rng.random() * grid.dz'
    #         z_function = eval(z_lambda_function, env)
    #     elif self.modes['z'] == 'surface_up':
    #         z_function = self.domain._get_elevation_from_surface_up
    #     elif self.modes['z'] == 'surface_down':
    #         z_function = self.domain._get_elevation_from_surface_down
    #     else:
    #         z_function = eval(z_lambda_function, env)
        
    #     return {'x' : x_function, 'y' : y_function, 'z' : z_function}


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

    # def _generate_coordinate(self) -> tuple:
    #     """ """
    #     # Intializes some variables
    #     is_coordinate_valid = False
    #     iteration = 0
    #     iteration_limit = 10000
        
    #     # Tries to generate a valid coordinate
    #     while not is_coordinate_valid:
    #         coordinate = self._generate_3D_coordinate()
    #         is_coordinate_valid = self._is_coordinate_3D_valid(coordinate)
            
    #         # TODO
    #         # security
    #         if iteration > iteration_limit:
    #             # LOG
    #             print(coordinate)
    #             raise ValueError('_generate_coordinate - iteration_limit exceeded')

    #         iteration += 1
    #     return coordinate

    # def _generate_3D_coordinate(self) -> tuple:
    #     """"""
    #     # Generates randomly a coordinate
    #     x = self.functions['x']()
    #     y = self.functions['y']()
    #     if self.modes['z'] in ['surface_up', 'surface_down']:
    #         z = self.functions['z'](x,y)
    #     else:
    #         z = self.functions['z']()
    #     return (x,y,z)
        
    # def _generate_3D_coordinate_from_2D_coordinate(self, coordinate:tuple) -> tuple:
    #     """"""
    #     x, y = coordinate
    #     if self.modes['z'] in ['surface_up', 'surface_down']:
    #         z = self.functions['z'](x,y)
    #     else:
    #         z = self.functions['z']()
    #     return (x,y,z)
        
    # def _is_coordinate_2D_valid(self, coordinate:tuple) -> bool:
    #     """Checks if the 2D point is valid"""
    #     x, y = coordinate
    #     if self.domain.is_coordinate_2D_valid(x,y):
    #         return True  
    #     else:
    #         return False
        
    # def _is_coordinate_3D_valid(self, coordinate:tuple) -> bool:
    #     """ 
    #     Checks if the points are well located ... 
    #     Two cases of coordinate validation:
    #     1 - No geology
    #     2 - With geology
    #     """
    #     x, y, z = coordinate
            
    #     # Case 1 - No geology
    #     if self.geologic_ids is None:
    #         return self.domain.is_coordinate_in_domain(x,y,z)
        
    #     # Case 2 - With geology
    #     else:
    #         i, j, k = self.domain.grid.get_indices(x,y,z)
    #         return self.domain.is_coordinate_in_domain(x,y,z) and (self.geology.data_volume[i,j,k] in self.geologic_ids)  