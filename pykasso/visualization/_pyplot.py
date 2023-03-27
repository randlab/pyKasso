"""
TODO
"""

### External dependencies
import numpy as np
import matplotlib.pyplot as plt

### Local dependencies
from .._utils import ProjectReader


class PyplotVisualizer(ProjectReader):
    """
    TODO
    """
    
    def __init__(self, project_directory: str, *args, **kwargs):
        """
        TODO
        """
        super().__init__(project_directory, *args, **kwargs)
        
    

# def _show_data(environment, feature, settings, show=True):
#     """ """
#     nx, ny, nz = environment.grid.shape
    
#     # Get the grid
#     grid = np.ones((nx, ny, nz))
    
#     ax = plt.figure().add_subplot(projection='3d')
#     ax.voxels(grid, edgecolor='k')

#     plt.show()
#     return None

# def _show_grid(environment):
#     return None

# def _show_mask(environment):
#     return None

# def _show_geology(environment):
#     return None



# #     def show_catchment(self, label='geology', title=None, cmap='binary'):
# #         """
# #         Show the entire study domain.
# #
# #         Parameters
# #         ----------
# #         label : str, optional
# #             Data to show : 'geology', 'topography', 'orientationx', 'orientationy' 'faults' or 'fractures'.
# #             By default : 'geology'.
# #         title : str, optional
# #             Title of the plot. If 'None', 'data' becomes the label.
# #         cmap : str, optional
# #             Color map, 'binary' by default.
# #         """
# #         import matplotlib.patches as mp
# #         fig, ax1 = plt.subplots()
# #         #if title is None:
# #         #    title = label
# #         #fig.suptitle(title, fontsize=16)
# #
# #         # Load data to show
# #         try:
# #             data = [data.data[:,:,0] for data in self.geology.data if data.label==label][-1]
# #         except:
# #             print('no data for indicated label parameter')
# #
# #         im1 = ax1.imshow(data.T, origin="lower", extent=self.grid.extent, cmap=cmap)
# #
# #         fig.colorbar(im1, ax=ax1)
# #         if self.settings['data_has_mask']:
# #             import matplotlib.patches as mp
# #             p1 = mp.PathPatch(self.mask.polygon, lw=2, fill=0, edgecolor='red', label='mask')
# #             ax1.add_patch(p1)
# #
# #         # Points
# #         for pts in self.points.points:
# #             x, y = zip(*pts.points)
# #             ax1.plot(x, y, 'o', label=pts.points_key)
# #
# #         ax1.set_aspect('equal', 'box')
# #         plt.legend(loc='upper right')
# #         plt.show()
# #         return fig