"""
TODO
"""

import numpy as np

# TODO
# - Documentation
# - Gestion des cartes de densité de probabilités pour le tirage des points
# - Fonction z

def generate_random_points(number, rng, grid, mask, geology, constraints=None):
    """
    TODO

    Generates random points according to the parameters.

    Four cases of points generation:
    1 - No geology, no mask
    2 - No geology, mask
    3 - Geology, no mask
    4 - Geology, mask
    """

    ### Lambda functions
    random_x = lambda : grid.x0 - grid.dx/2 + grid.nx * rng.random() * grid.dx
    random_y = lambda : grid.y0 - grid.dy/2 + grid.ny * rng.random() * grid.dy
    rand_x, rand_y = [], []
    validated_points = 0
    validated_geology = []
    surface = geology.surface

    # TODO - changer de position + logging ?
    ### Inspects validity of geological constraints
    # if constraints is not None:
        # for geological_id in constraints['geology']:
            # pass
            # if geological_id in 



    # TODO (1))
    # dictionnaire des geology présentes à la surface
    # for gid in self.geology_ID:
    #     assert gid in surf, "/!\\ - create_point_feature() - ERROR : geologic id n°{} is not on the surface geologic model.".format(gid)

    # 1 - No geology, no mask
    if (constraints is None) and (mask is None):
        rand_x = [random_x() for x in range(number)]
        rand_y = [random_y() for y in range(number)]

    # 2 - No geology, mask
    elif (constraints is None):
        while validated_points < number:
            x, y = random_x(), random_y()
            if int(mask.polygon.contains_point((x,y))):
                rand_x.append(x)
                rand_y.append(y)
                validated_points += 1

    # 3 - Geology, no mask
    elif (mask is None):
        while validated_points < number:
            x, y = random_x(), random_y()
            if surface[grid.get_i(x)][grid.get_j(y)] in constraints:
                rand_x.append(x)
                rand_y.append(y)
                # rand_z.append(z)
                validated_points += 1

    # 4 - Geology, mask
    else:
        while validated_points < number:
            x, y = random_x(), random_y()
            if (surface[grid.get_i(x)][grid.get_j(y)] in constraints) and (int(mask.polygon.contains_point((x,y)))):
                rand_x.append(x)
                rand_y.append(y)
                validated_points += 1

    return list(zip(rand_x, rand_y))


def inspect_points_validity(points, grid, mask):
    """
    TODO
    Checks if the points are well located inside the grid, or well inside the mask if one is provided.
    Prints out only when an issue is encountered.
    """

    if mask is None:
        tests = grid.path.contains_points(points)
    
    else:
        tests = grid.path.contains_points(points) & mask.polygon.contains_points(points)

    validated_points = [point for (point, test) in zip(points, tests) if test == True]

    return validated_points