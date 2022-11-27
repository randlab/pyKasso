"""
TODO
"""

import numpy as np

def generate_random_points(number, grid, mask, geology, constraints=None):
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
    random_x = lambda : grid.x0 - grid.dx/2 + grid.nx * np.random.random() * grid.dx
    random_y = lambda : grid.y0 - grid.dy/2 + grid.ny * np.random.random() * grid.dy
    # random_z = lambda : grid.z0 - grid.dz/2 + grid.nz * np.random.random() * grid.dz
    # fixed_z  = lambda : grid.zmax - 1
    # rand_x, rand_y, rand_z = [], [], []
    rand_x, rand_y = [], []
    validated_points = 0
    surface = geology.surface

    # TODO
    # dictionnaire des geology présentes à la surface
    # for gid in self.geology_ID:
    #     assert gid in surf, "/!\\ - create_point_feature() - ERROR : geologic id n°{} is not on the surface geologic model.".format(gid)

    # ##### z points option generation ####TODO####
    # z_option = 'random' # or 'surface'
    # if z_option == 'random':
    #     z_func = random_z
    # else:
    #     z_func = fixed_z
    # ###################################

    # 1 - No geology, no mask
    if (constraints is None) and (mask is None):
        rand_x = [random_x() for x in range(number)]
        rand_y = [random_y() for y in range(number)]
        # rand_z = [z_func()   for z in range(number)]

    # 2 - No geology, mask
    elif (constraints is None):
        while validated_points < number:
            # x, y, z = random_x(), random_y(), z_func()
            x, y = random_x(), random_y()
            if int(mask.polygon.contains_point((x,y))):
                rand_x.append(x)
                rand_y.append(y)
                # rand_z.append(z)
                validated_points += 1

    # 3 - Geology, no mask
    elif (mask is None):
        while validated_points < number:
            # x, y, z = random_x(), random_y(), z_func()
            x, y = random_x(), random_y()
            if surface[grid.get_i(x)][grid.get_j(y)] in constraints:
                rand_x.append(x)
                rand_y.append(y)
                # rand_z.append(z)
                validated_points += 1

    # 4 - Geology, mask
    else:
        while validated_points < number:
            # x, y, z = random_x(), random_y(), z_func()
            x, y = random_x(), random_y()
            if (surface[grid.get_i(x)][grid.get_j(y)] in constraints) and (int(mask.polygon.contains_point((x,y)))):
                rand_x.append(x)
                rand_y.append(y)
                # rand_z.append(z)
                validated_points += 1

    return list(zip(rand_x, rand_y))

def inspect_points_validity(points):
    """
    TODO
    Checks if the points are well located inside the grid, or well inside the mask if one is provided.
    Prints out only when an issue is encountered.
    """

    return None











    # from IPython.display import display, HTML

    # # Analyse points locations outside domain delimitations
    # df.loc[:, 'out_domain'] = (df['x'] > self.grid.xmax) | (df['x'] < self.grid.xmin) | (df['y'] > self.grid.ymax) | (df['y'] < self.grid.ymin) | (df['z'] > self.grid.zmax) | (df['z'] < self.grid.zmin)
    # df_out_domain = df[df['out_domain'] == True]
    # df_in_domain  = df[df['out_domain'] == False]
    # nb_out_domain = len(df_out_domain)
    # nb_in_domain  = len(df_in_domain)

    # # Set return value
    # df = df_in_domain[['label', 'x', 'y', 'z']]

    # # Analyse points locations outside mask delimitations but inside domain delimitations
    # if (self.mask is not None) and (nb_in_domain != 0):
    #     points = list(zip(df_in_domain['x'], df_in_domain['y']))
    #     test = pd.Series(self.mask.polygon.contains_points(points))
    #     df_in_domain['out_mask'] = ~test

    #     df_out_mask = df_in_domain[df_in_domain['out_mask'] == True]
    #     df_in_mask  = df_in_domain[df_in_domain['out_mask'] == False]
    #     nb_out_mask = len(df_out_mask)
    #     nb_in_mask  = len(df_in_mask)

    #     df = df_in_mask[['label', 'x', 'y', 'z']]

    #     if (nb_out_mask > 0):
    #         print('/!\\ _inspect_points() - WARNING : Point(s) inside domain but outside mask:')
    #         display(df_out_mask)

    #     if (nb_in_mask == 0):
    #         print('/!\\ _inspect_points() - ERROR : All the points are outside mask:')
    #         sys.exit()

    # if (nb_out_domain > 0):
    #     print('/!\\ _inspect_points() - WARNING : Point(s) out of domain:')
    #     display(df_out_domain)

    # if (nb_in_domain == 0):
    #     print('/!\\ _inspect_points() - ERROR : All the points are out of domain:')
    #     sys.exit()

    # return df