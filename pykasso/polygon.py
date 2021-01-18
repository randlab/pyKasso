from .functions      import opendatafile, loadpoints
from matplotlib.path import Path

import numpy             as np
import matplotlib.pyplot as plt

class Polygon():
    """
	Create a polygon manager : a class for managing a polygon as an area study delimitor.

	Parameters
	----------
	grid : Grid()
		Polygon() class needs a Grid() object as argument.
	"""

    def __init__(self, grid):
        self.grid    = grid
        self.polygon = None
        self.mask    = None

    def set_polygon(self, vertices):
        """
        Create a polygon from vertices coordinates.

        Parameters
        ----------
        vertices : str || list
            Location of the datafile or list of the vertices coordinates like [[x1,y1],[x2,y2],[xn,yn]].
        """
        if isinstance(vertices, str):
            text     = opendatafile(vertices)
            vertices = loadpoints(text)
        
        if len(vertices) < 3:
            print("- set_polygon() - Error : Not enough vertices to create a polygon (3 minimum).")
        else:
            self.polygon = vertices
            self.mask = self._set_mask()
        return None

    def _set_mask(self):
        """
        Set the mask.
        """
        mask = np.zeros((self.grid.nx, self.grid.ny))
        for y in range(self.grid.ny):
            for x in range(self.grid.nx):
                mask[x][y] = not int(Path(self.polygon).contains_point((self.grid.X[x][y][0],self.grid.Y[x][y][0])))
        return mask

    def inspect_polygon(self):
        """
        Check if the vertices of the polygon are located inside the grid or not.
        If all the vertices are outside the grid, the polygon is rested to none.

        Returns
        -------
            List of vertices out of the grid else 'None'.
        """
        if self.polygon is not None:
            unvalidated_vertices = []
            xlimits = [self.grid.xlimits[0], self.grid.xlimits[0], self.grid.xlimits[1], self.grid.xlimits[1]]
            ylimits = [self.grid.ylimits[0], self.grid.ylimits[1], self.grid.ylimits[1], self.grid.ylimits[0]]
            for k,(x,y) in enumerate(self.polygon):
                if not int(Path(list(zip(xlimits, ylimits))).contains_point((x,y))):
                    unvalidated_vertices.append(k+1)
            validated_vertices = len(self.polygon) - len(unvalidated_vertices)

            if len(unvalidated_vertices) == len(self.polygon):
                print('- inspect_polygon() - Warning : 0 vertices inside the grid limits. Polygon is reseted to none.'.format(len(unvalidated_vertices),len(self.polygon)))
                self.polygon = None
                self.mask    = None
                return None

            if len(unvalidated_vertices) > 0:
                print('- inspect_polygon() - Warning : {} vertices not inside the grid limits on {} vertices.'.format(len(unvalidated_vertices),len(self.polygon)))
                for vertex in unvalidated_vertices:
                    print('- vertice {}'.format(vertex))
                return unvalidated_vertices
        else:
            print('- inspect_polygon() - Error : no polygon to inspect. Please set a polygon.')
        return None

    def clean_polygon(self):
        """
        Remove the polygon.
        """
        self.polygon = None
        self.mask    = None

    def show(self):
        """
        Show the delimitation of the grid and, if a polygon is present, display its limits.
        """
        if self.polygon is not None:
            closed_polygon = self.polygon[:]
            closed_polygon.append(closed_polygon[0])
            x,y = zip(*closed_polygon)

            fig, ax = plt.subplots()
            fig.suptitle('Show polygon', fontsize=16)
            ax.plot(x,y, color='black')

            xlimits = [self.grid.xlimits[0], self.grid.xlimits[0], self.grid.xlimits[1], self.grid.xlimits[1], self.grid.xlimits[0]]
            ylimits = [self.grid.ylimits[0], self.grid.ylimits[1], self.grid.ylimits[1], self.grid.ylimits[0], self.grid.ylimits[0]]
            
            ax.plot(xlimits, ylimits, color='red')
            ax.set_aspect('equal', 'box')
            plt.legend(('polygon','grid limits'), loc='center left', bbox_to_anchor=(1, 0.5))
            plt.show()
            return None
        else:
            print('- show_Polygon() - Error : no polygon to show.')
            return None