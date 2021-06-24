from .functions      import opendatafile, loadpoints
from matplotlib.path import Path

import sys
import numpy             as np
import numpy.ma          as ma
import matplotlib.pyplot as plt

class PointManager():
    """
    Class modeling the inlets and the outlets of the studied karst network's domain.
    """

    def __init__(self, grid, polygon, geology=None):
        """
        Creates a point manager on a Grid instance.
        This class is designed to handle points and to transform them into inlets or outlets.

        Parameters
        ----------
        grid : Grid instance
            The point manager must be set on the studied grid.
        polygon : Polygon instance
            The point manager needs a Polygon instance even if it's empty.
        geology : GeologyManager instance, optional
	       The point manager can use geologic data as parameters.

        Examples
        --------
        >>> grid = pk.Grid(0, 0, 0, 10, 10, 10, 1, 1, 1)
        >>> poly = pk.Polygon(grid)
        >>> geol = pk.GeologyManager(grid)
        >>> pts  = pk.PointManager(grid, poly, geol)
        """
        self.points  = {}
        self.grid    = grid
        self.polygon = polygon
        self.geology = geology

    def set_points(self, points_key, points):
        """
        Sets new points for the indicated points key.

        Parameters
        ----------
        points_key : str
            Type of points : 'inlets' or 'outlets'.
        points : str || array
            Location of the datafile or array of the points coordinates like [[x1,y1],[x2,y2],[xn,yn]].

        Examples
        --------
        >>> pts.set_points("inlets", [[0,0], [10,10]])
        >>> pts.set_points("outlets", "outlets.csv")
        """
        if isinstance(points, str):
            text = opendatafile(points)
            self.points[points_key] = loadpoints(text)
        else:
            self.points[points_key] = points
        return None

    def generate_points(self, points_key, points_number, geological_IDs=None):
        """
        Generates random points on the grid according to the parameters.

        Parameters
        ----------
        points_key : str
            Type of points : 'inlets' or 'outlets'.
        points_number : int
            Number of points to generate.
        geological_IDs : array, optional
            List of geologic facies where points can be created.

        Examples
        --------
        >>> pts.generate_points("inlets", 20)
        >>> pts.generate_points("outlets", 3, geological_IDs=[0,2])
        """
        random_x = lambda : self.grid.x0 - self.grid.dx/2 + self.grid.nx * np.random.random() * self.grid.dx
        random_y = lambda : self.grid.y0 - self.grid.dy/2 + self.grid.ny * np.random.random() * self.grid.dy

        # With geology, without polygon
        if geological_IDs is None:
            if self.polygon.polygon is None:
                rand_x = [random_x() for x in range(points_number)]
                rand_y = [random_y() for y in range(points_number)]
                self.points[points_key] = list(zip(rand_x, rand_y))
        # Without geology, with polygon
            else:
                validated_inlets = 0
                rand_x, rand_y = [], []
                while validated_inlets < points_number:
                    x, y = random_x(), random_y()
                    if int(self.polygon.polygon.contains_point((x, y))):
                        rand_x.append(x)
                        rand_y.append(y)
                        validated_inlets += 1
                self.points[points_key] = list(zip(rand_x, rand_y))

        # With geology
        else:
            if self.geology is None:
                print('- generate_points() - ERROR : no geology to consider in order to generate points.')
                sys.exit()
            if self.geology.data["geology"]["stats"] is None:
                self.geology.compute_stats_on_data("geology")
            # Keep only present ids
            ids = self.geology.data["geology"]["stats"]["ID"]
            new_ids = [id for id in geological_IDs if (id in ids)]
            if not new_ids:
                 print('- generate_points() - ERROR : none of the geology ids provided are valid.')
                 sys.exit()
            validated_inlets = 0
            rand_x, rand_y = [], []
            # Without polygon
            if self.polygon.polygon is None:
                while validated_inlets < points_number:
                    x, y = random_x(), random_y()
                    if self.geology.data["geology"]["data"][self.grid.get_i(x)][self.grid.get_j(y)][0] in new_ids:
                        rand_x.append(x)
                        rand_y.append(y)
                        validated_inlets += 1
                self.points[points_key] = list(zip(rand_x,rand_y))
            # With polygon
            else:
                while validated_inlets < points_number:
                    x, y = random_x(), random_y()
                    if self.geology.data["geology"]["data"][self.grid.get_i(x)][self.grid.get_j(y)][0] in new_ids:
                        if int(self.polygon.polygon.contains_point((x,y))):
                            rand_x.append(x)
                            rand_y.append(y)
                            validated_inlets += 1
                self.points[points_key] = list(zip(rand_x,rand_y))
        return None

    def composite_points(self, points_key, points, points_number, geological_IDs=None):
        """
        Generates random points and combines it with user indicated points.

        Parameters
        ----------
        points_key : str
            Type of points : 'inlets' or 'outlets'.
        points : str || array
            Location of the datafile or list of the points coordinates.
        points_number : int
            Number of points to generate.
        geological_IDs : array, optional
            List of geologic facies where points can be created.

        Examples
        --------
        >>> pts.composite_points("inlets", "inlets.csv", 20)
        >>> pts.composite_points("outlets", [[25,25]], 3, geological_IDs=[1])
        """
        self.generate_points(points_key, points_number, geological_IDs)
        if isinstance(points, str):
            text = opendatafile(points)
            other_points = loadpoints(text)
            self.points[points_key] += other_points
        else:
            self.points[points_key] += points
        return None

    def inspect_points(self):
        """
        Checks if the points are well located inside the grid, or well inside the polygon if one is provided.
        Prints out only when an issue is encountered.

        Examples
        --------
        >>> pts.inspect_points()
        """
        if self.points is None:
            print('- inspect_points() - ERROR : no points to inspect.')
            sys.exit()
        else:
            for key in self.points:
                if self.polygon.polygon is not None:
                    logicals = np.logical_and(self.polygon.polygon.contains_points(self.points[key]), self.grid.path.contains_points(self.points[key]))
                else:
                    logicals = self.grid.path.contains_points(self.points[key])

                for i, value in enumerate(logicals):
                    if value is True:
                        print('- inspect_points() - WARNING : point n°{} {} removed because not inside inside the polygon or the domain.'.format(i, self.points[key][i]))

                self.points[key] = [list(point) for (point, logical) in zip(self.points[key], logicals) if logical]

                if not self.points[key]:
                    print('- inspect_points() - ERROR : all the {} have been removed '.format(key))
                    sys.exit()
            return None

    def show(self):
        """
        Shows the delimitations of the grid, of the polygon (if present) and of the locations of the points (if present).

        Examples
        --------
        >>> pts.show()
        """
        import matplotlib.patches as mp

        fig, ax = plt.subplots()
        fig.suptitle('Points', fontsize=16)

        # Geology
        if self.geology is not None:
            try:
                d = self.geology.data['geology']['img'][:,:,0]
                plt.imshow(d, extent=self.grid.extent, cmap='gray_r')
            except:
                None

        # Grid limits
        p1 = mp.PathPatch(self.grid.path, lw=2, fill=0, edgecolor='red', label='grid limits')
        ax.add_patch(p1)

        # Polygon
        if self.polygon.polygon is not None:
            p2 = mp.PathPatch(self.polygon.polygon, lw=2, fill=0, edgecolor='orange', label='polygon')
            ax.add_patch(p2)

        # Points
        for key in self.points:
            x, y = zip(*self.points[key])
            ax.plot(x, y, 'o', label=key)

        plt.xlim((self.grid.xlimits[0] - self.grid.dx/2, self.grid.xlimits[1] + self.grid.dx/2))
        plt.ylim((self.grid.ylimits[0] - self.grid.dy/2, self.grid.ylimits[1] + self.grid.dy/2))

        ax.set_aspect('equal', 'box')
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.show()
        return fig
