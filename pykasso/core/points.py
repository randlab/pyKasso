"""
TODO
"""

import logging
import numpy as np

# TODO
# - Documentation
# - Validation des variables
# - Gestion des cartes de densité de probabilités pour le tirage des points


def generate_coordinate(lambda_functions, rng, grid, mask, geology, xi=None, yi=None, zi=None, geologic_ids=None):
    """
    TODO

    Generates random points according to the parameters.

    Four cases of points generation:
    1 - No geology, no mask
    2 - No geology, mask
    3 - Geology, no mask
    4 - Geology, mask
    """

    # Evaluates the function
    env = {
        'rng'  : rng,
        'grid' : grid,
        # TODO - define mask limits
        # 'mask' : mask,
        # 'geology' : geology
    }
    x_func = eval(lambda_functions['x'], env)
    y_func = eval(lambda_functions['y'], env)
    z_func = eval(lambda_functions['z'], env)

    # Controls validity of 'geologic_ids'
    if geologic_ids is not None:
        values = geology.stats.index.to_list()
        validated_geology_ids = []
        for geologic_id in geologic_ids:
            if geologic_id in values:
                validated_geology_ids.append(geologic_id)
            else:
                pass
                # TODO - log

        if len(validated_geology_ids) == 0:
            pass
            # TODO - log - erreur


    # Tries to generate a valid coordinate
    coordinate_is_valid = False
    iteration = 0
    iteration_limit = 10000
    while not coordinate_is_valid:
        
        # TODO
        # security
        if iteration > iteration_limit:
            # LOG
            raise ValueError

        # Generates a coordinate
        if xi == None:
            x = x_func()
        else:
            x = xi

        if yi == None:
            y = y_func()
        else:
            y = yi

        if zi == None:
            z = z_func()
        else:
            z = zi


        # 1 - No geology, no mask
        if geologic_ids is None: 
            if mask is None:
                coordinate_is_valid = is_point_valid((x,y,z), grid, mask)
        
        # 2 - No geology, mask
            else:
                if int(mask.polygon.contains_point((x,y))) and is_point_valid((x,y,z), grid, mask):
                    coordinate_is_valid = True

        
        else:
            if is_point_valid((x,y,z), grid, mask):
                i, j, k = grid.get_i(x), grid.get_j(y), grid.get_k(z)

        # 3 - Geology, no mask
                if (mask is None):
                    if geology.data[i][j][k] in validated_geology_ids:
                        coordinate_is_valid = True

        # 4 - Geology, mask
                else:
                    if (geology.data[i][j][k] in validated_geology_ids) and (int(mask.polygon.contains_point((x,y)))):
                        coordinate_is_valid = True

        iteration += 1

    return [x,y,z]


def is_point_valid(point, grid, mask):
    """
    TODO
    Checks if the points are well located inside the grid, or well inside the mask if one is provided.
    """
    if len(point) == 2:

        x, y = point

        if mask is None:
            test = grid.path.contains_point((x, y))
        
        else:
            test = grid.path.contains_point((x, y)) & mask.polygon.contains_point((x, y))

    if len(point) == 3:

        x, y, z = point

        if mask is None:
            test = grid.is_inbox(x, y, z)

        else:
            test = grid.is_inbox(x, y, z) & mask.polygon.contains_point((x, y))

    return test


# def get_lambda_function(dim, keyword):
#     """
#     TODO
#     """
#     lambda_functions = {
#         'x' : {
#             'surface' : ,
#             'bottom'  : ,
#         },
#         'y' : {

#         },
#         'z' : {

#         }
#     }

#     return lambda_function