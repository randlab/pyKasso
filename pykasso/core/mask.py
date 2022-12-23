"""
Module to model the mask.
"""

import sys
import copy
import numpy as np
from matplotlib.path import Path

class Mask():
    """
    Class modeling the catchment delimitation of the studied domain.
    """

    def __init__(self, data):
        """
        Creates a mask.
        This class is designed to describe a particular catchment delimitation.

        Parameters
        ----------
        data : list | str
            List of vertices from Python list or datafile location.
        name : str, optional
            Name of the mask.
        """
        self.validated_vertices = None
        self.polygon = None
        self.mask = None

        # If needed, loads the data
        if isinstance(data, str):
            self.vertices = np.genfromtxt(data)
        elif isinstance(data, list):
            self.vertices = np.array(data)
        else:
            # TODO
            raise
         

    def validate_vertices(self, grid):
        """
        TODO
        """
        # controls if the vertices are well defined inside the grid limits
        vertices = self.vertices
        validated_vertices = [k for k,(x,y) in enumerate(vertices) if (grid.path.contains_point((x,y)) == True)]

        if len(validated_vertices) == 0:
            ## TODO
            print("MASK - ERROR : No vertices validated inside the grid limits.")
            sys.exit()

        if len(validated_vertices) < len(vertices):
            # TODO if debug:
            print('MASK - WARNING : {}/{} vertices validated inside the grid limits.'.format(len(validated_vertices), len(vertices)))

        self.validated_vertices = vertices[validated_vertices]

        # sets the polygon with a matplotlib Path
        dc_vertices = copy.deepcopy(vertices).tolist()
        dc_vertices.append(dc_vertices[0])
        self.polygon = Path(dc_vertices)

        # sets the mask array with numpy-array
        row, col = np.indices((grid.nx, grid.ny))
        pts = np.column_stack((grid.X[row, col, 0].flatten(), grid.Y[row, col, 0].flatten()))
        msk = self.polygon.contains_points(pts).reshape((grid.nx, grid.ny)).astype(int)
        self.mask = np.repeat(msk[:, :, np.newaxis], grid.nz, axis=2)

        return None

    def __repr__(self):
        txt = "Name : " + str(self.name) + " Nbr of vertices : " + str(len(self.validated_vertices))
        return txt