#############################################################################################################################
### pyvista ###
###############

# def _voxelize_fractures_pyvista(grid, fractures):
#     """
#     TODO
#     """

#     ##############################
#     ### 1 - Construct the rays ###
#     ##############################
#     slice_x = slice(0, grid.nx)
#     slice_y = slice(0, grid.ny) 
#     slice_z = slice(0, grid.nz)
#     val_s, val_f = 0, -1
#     SLICES = [ # xz, zy, yx
#         [[slice_x, val_s  , slice_z], [slice_x, val_f  , slice_z]],
#         [[val_s  , slice_y, slice_z], [val_f  , slice_y, slice_z]],
#         [[slice_x, slice_y,   val_s], [slice_x, slice_y,   val_f]]
#     ]
#     DIMS = ['X', 'Y', 'Z']
#     POINTS_START = []
#     POINTS_FINAL = []
#     # RAYS = []

#     for (slice_start, slice_final) in SLICES:

#         points_start = []
#         points_final = []

#         for dim in DIMS:
#             i_s, j_s, k_s = slice_start
#             coordinates = grid._get_property(dim)[i_s, j_s, k_s]
#             coordinates = coordinates.flatten(order='F')
#             points_start.append(coordinates)

#             i_f, j_f, k_f = slice_final
#             coordinates = grid._get_property(dim)[i_f, j_f, k_f]
#             coordinates = coordinates.flatten(order='F')
#             points_final.append(coordinates)

#         points_start = [list(tpl) for tpl in zip(*points_start)]
#         points_final = [list(tpl) for tpl in zip(*points_final)]

#         POINTS_START.append(points_start)
#         POINTS_FINAL.append(points_final)

#         # for (start, final) in zip(points_start, points_final):
#         #     RAYS.append(pv.Line(start, final))

#     # flat the lists
#     POINTS_START = [item for sublist in POINTS_START for item in sublist]
#     POINTS_FINAL = [item for sublist in POINTS_FINAL for item in sublist]
#     # convert in numpy arrays
#     POINTS_START = np.array(POINTS_START)
#     POINTS_FINAL = np.array(POINTS_FINAL)

#     #####################################################
#     ### 2 - Transform fractures into pyvista polygons ###
#     #####################################################
#     POLYGONS = []
#     polygon_definition = 8 # TODO - Memory settings
#     for fracture in fractures:
#         x, y, z = fracture.get_position()
#         a, b, c = fracture.get_normal()
#         rad     = fracture.radius
#         POLYGONS.append(pv.Polygon(center=(x, y, z), radius=rad, normal=(a, b, c), n_sides=polygon_definition))

#     ###############################
#     ### 3 - Intersection points ###
#     ###############################
#     POINTS = []
#     for polygon in POLYGONS:
#         for (start, final) in zip(POINTS_START, POINTS_FINAL):
#             points, ind = polygon.ray_trace(start, final)
#             if points.size != 0:
#                 POINTS.append(points)

#     ############################
#     ### 4 - Retrieve indices ###
#     ############################
#     INDICES = []
#     for p in POINTS:
#         x, y, z = p[0]
#         i, j, k = grid.get_i(x), grid.get_j(y), grid.get_k(z)
#         INDICES.append([i, j, k])

#     # Remove duplicates
#     INDICES = np.array(INDICES)
#     INDICES = np.unique(INDICES, axis=0)
#     i, j, k = zip(*INDICES)

#     ############################
#     ### 5 - Fracturation map ###
#     ############################
#     fracs_array = np.zeros((grid.nx, grid.ny, grid.nz))
#     fracs_array[i, j, k] = 1

#     return fracs_array