#     ######################
#     ### Karst Analysis ###
#     ######################

#     def compare_stats(self, mean=False):
#         """
#         Compare statistics between reference indicators and calculated networks.

#         TODO
#         param 'iteration=0'

#         """
#         indicators = ['cpd', 'cv degree', 'cv length', 'orientation entropy', 'length entropy', 'aspl', 'mean degree', 'mean length', 'correlation vertex degree']
#         stats = pd.DataFrame(columns=indicators)

#         for (i, karst_network) in enumerate(self.SIMULATIONS):
#             stats.loc[i] = karst_network.stats

#         # Apply style
#         def _bg_color(x, min_val, max_val):
#             if (x < min_val) or (x > max_val):
#                 return 'background-color: red'
#             else:
#                 return 'background-color: green'

#         display(stats.style.applymap(_bg_color, min_val = self.reference_statistics['cpd']['min'], max_val = self.reference_statistics['cpd']['max'], subset = ['cpd'])\
#         .applymap(_bg_color, min_val = self.reference_statistics['cv degree']['min'], max_val = self.reference_statistics['cv degree']['max'], subset = ['cv degree'])\
#         .applymap(_bg_color, min_val = self.reference_statistics['cv length']['min'], max_val = self.reference_statistics['cv length']['max'], subset = ['cv length'])\
#         .applymap(_bg_color, min_val = self.reference_statistics['orientation entropy']['min'], max_val = self.reference_statistics['orientation entropy']['max'], subset = ['orientation entropy'])\
#         .applymap(_bg_color, min_val = self.reference_statistics['length entropy']['min'], max_val = self.reference_statistics['length entropy']['max'], subset = ['length entropy'])\
#         .applymap(_bg_color, min_val = self.reference_statistics['aspl']['min'], max_val = self.reference_statistics['aspl']['max'], subset = ['aspl'])\
#         .applymap(_bg_color, min_val = self.reference_statistics['mean degree']['min'], max_val = self.reference_statistics['mean degree']['max'], subset = ['mean degree'])\
#         .applymap(_bg_color, min_val = self.reference_statistics['mean length']['min'], max_val = self.reference_statistics['mean length']['max'], subset = ['mean length'])\
#         .applymap(_bg_color, min_val = self.reference_statistics['correlation vertex degree']['min'], max_val = self.reference_statistics['correlation vertex degree']['max'], subset = ['correlation vertex degree']))
#         return None


#     def compute_average_paths(self, mask=0):
#         """
#         TODO
#         Compute the mean of all the simulations.
#         """

#         # Calculate the average from all the simulations
#         karst_maps = []
#         for karst_simulation in self.SIMULATIONS:
#             data = karst_simulation.maps['karst'][-1]
#             karst_maps.append(data)
#         karst_prob = sum(karst_maps)/len(karst_maps)

#         self.karst_prob = karst_prob
#         return karst_prob

#     def show_average_paths(self):
#         """
#         todo
#         """
#         ### Call the plotter
#         p = pv.Plotter(notebook=False)

#         ### Construct the grid
#         vtk = pv.UniformGrid()
#         vtk.dimensions = np.array((self.GRID.nx, self.GRID.ny, self.GRID.nz)) + 1
#         vtk.origin     = (self.GRID.x0 - self.GRID.dx/2, self.GRID.y0 - self.GRID.dy/2, self.GRID.z0 - self.GRID.dz/2)
#         vtk.spacing    = (self.GRID.dx, self.GRID.dy, self.GRID.dz)

#         vtk['values'] = self.karst_prob.flatten(order="F")

#         mesh = vtk.cast_to_unstructured_grid()
#         ghosts = np.argwhere(vtk['values'] < 1.0)
#         mesh.remove_cells(ghosts)
#         p.add_mesh(mesh, show_edges=False)

#         ### Plotting
#         # p.add_title(feature)
#         p.add_axes()
#         bounds = p.show_bounds(mesh=vtk)
#         p.add_actor(bounds)
#         p.show(cpos='xy')

#         return None