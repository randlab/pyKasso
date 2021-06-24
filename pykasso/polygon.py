from .functions       import opendatafile, loadpoints
from matplotlib.path  import Path

import sys
import copy
import numpy             as np
import matplotlib.pyplot as plt

class Polygon():
    """
    Class modeling the catchment delimitation of the studied domain.
    """

    def __init__(self, grid):
        """
        Creates a polygon manager on a Grid instance.

        Parameters
        ----------
        grid : Grid instance
            The polygon must be set on the studied grid.

        Examples
        --------
        >>> grid = pk.Grid(0, 0, 0, 10, 10, 10, 1, 1, 1)
        >>> poly = pk.Polygon(grid)
        """
        self.grid     = grid
        self.vertices = None
        self.polygon  = None
        self.mask     = None

    def set_polygon(self, vertices):
        """
        Creates a polygon from vertices coordinates.
        The polygon is stored in the 'polygon' attribute.
        It creates a mask and stores it in the 'mask' attribute.
        It stores the polygon surface in a 'surface' attribute.

        Parameters
        ----------
        vertices : str || array
            Location of the datafile, or array of the vertices coordinates like [[x1,y1],[x2,y2],[xn,yn]].

        Examples
        --------
        >>> poly = pk.Polygon(grid)
        >>> poly.set_polygon([[0,0], [0,10], [10,0]])
        >>> poly.set_polygon("polygon.csv")
        >>> print(poly.polygon)
        >>> print(poly.mask)
        """
        if isinstance(vertices, str):
            text     = opendatafile(vertices)
            vertices = loadpoints(text)

        if len(vertices) < 3:
            print("- set_polygon() - Error : Not enough vertices to create a polygon (3 minimum).")
        else:
            self.vertices = vertices
            dc_vertices = copy.deepcopy(vertices)
            dc_vertices.append(vertices[0])
            self.polygon = Path(dc_vertices)
            self.mask    = self._set_mask()
            self._inspect_polygon()
        return None

    def _set_mask(self):
        """
        If a polygon is set, cells outside the polygon can be hidden with a mask.
        """
        mask = np.zeros((self.grid.nx, self.grid.ny, self.grid.nz), dtype=np.int_)
        for y in range(self.grid.ny):
            for x in range(self.grid.nx):
                mask[x][y][0] = not(int(self.polygon.contains_point((self.grid.X[x][y][0],self.grid.Y[x][y][0]))))
        for z in range(self.grid.nz):
            mask[:,:,z] = mask[:,:,0]
        return mask

    def _inspect_polygon(self):
        """
        Checks if the vertices of the polygon are located inside the grid or not.
        If all the vertices are outside the grid, the polygon is reset to none.
        """
        unvalidated_vertices = []

        for k,(x,y) in enumerate(self.vertices):
            if not int(self.grid.path.contains_point((x,y))):
                unvalidated_vertices.append(k)

        if len(unvalidated_vertices) == len(self.vertices):
            print('- inspect_polygon() - ERROR : 0 vertice inside the grid limits.')
            self.remove_polygon()
            sys.exit()

        if len(unvalidated_vertices) > 0:
            print('- inspect_polygon() - WARNING : {}/{} vertices not inside the grid limits.'.format(len(unvalidated_vertices), len(self.vertices)))
            for vertex in unvalidated_vertices:
                print('- vertex {} : {}'.format(vertex, self.vertices[vertex]))

        return None

    def remove_polygon(self):
        """
        Removes the polygon. Sets the 'polygon', the 'vertices' and the 'mask' attributes to none type.

        Examples
        --------
        >>> poly.remove_polygon()
        """
        self.polygon  = None
        self.vertices = None
        self.mask     = None

    def show(self):
        """
        Shows the delimitation of the grid and, if a polygon is present, displays its limits.

        Returns
        -------
        result : pyplot

        Examples
        --------
        >>> poly.show()
        """
        import matplotlib.patches as mp

        if self.vertices is not None:
            fig, ax = plt.subplots()
            fig.suptitle('polygon', fontsize=16)

            p1 = mp.PathPatch(self.grid.path, lw=2, fill=0, edgecolor='red')
            p2 = mp.PathPatch(self.polygon,   lw=2, fill=0, edgecolor='orange')
            ax.add_patch(p1)
            ax.add_patch(p2)

            plt.xlim((self.grid.xlimits[0] - self.grid.dx/2, self.grid.xlimits[1] + self.grid.dx/2))
            plt.ylim((self.grid.ylimits[0] - self.grid.dy/2, self.grid.ylimits[1] + self.grid.dy/2))

            ax.set_aspect('equal', 'box')
            plt.legend(('polygon','grid limits'), loc='center left', bbox_to_anchor=(1, 0.5))
            return fig
        else:
            print('- show() - ERROR : no polygon to show.')
            return None
