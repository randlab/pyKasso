from .functions      import opendatafile, loadpoints
from matplotlib.path import Path

import sys
import numpy             as np
import numpy.ma          as ma
import matplotlib.pyplot as plt

class PointManager():
    """
    Create a point manager : a class for managing points and transform them to inlets or outlets.

    Parameters
    ----------
    grid : Grid()
        PointManager() class needs a Grid() object as argument.
    polygon : Polygon()
        PointManager() class needs a Polygon() object as argument.
	geology : GeologyManager()
		PointManager() class needs a GeologyManager() object as argument.
    """

    def __init__(self, grid, polygon, geology=None):
        self.points  = {}
        self.grid    = grid
        self.polygon = polygon
        self.geology = geology

    def set_points(self, points_key, points):
        """
        Set new points.

        Parameters
        ----------
        points_key : str
            Type of points : 'inlets' or 'outlets'.
        points : str || list
            Location of the datafile or list of the points coordinates like [[x1,y1],[x2,y2],[xn,yn]].
        """
        if isinstance(points, str):
            text = opendatafile(points)
            self.points[points_key] = loadpoints(text)
        else:
            self.points[points_key] = points
        return None

    def generate_points(self, points_key, points_number, geological_IDs=None):
        """
        Generate random points on the grid, according to the parameters.

        Parameters
        ----------
        points_key : str
            Type of points : 'inlets' or 'outlets'.
        points_number : int
            Number of points to generate.
        geological_IDs : list
            List of geologic facies where points can be created.
        """
        if geological_IDs is None:
            if self.polygon.polygon is None:
                rand_x = [self.grid.x0 - self.grid.dx/2 + self.grid.nx * np.random.random() * self.grid.dx for x in range(points_number)]
                rand_y = [self.grid.y0 - self.grid.dy/2 + self.grid.ny * np.random.random() * self.grid.dy for y in range(points_number)]
                self.points[points_key] = list(zip(rand_x, rand_y))
            else:
                validated_inlets = 0
                rand_x = []
                rand_y = []
                while validated_inlets < points_number:
                    x = self.grid.x0 - self.grid.dx/2 + self.grid.nx*np.random.random() * self.grid.dx
                    y = self.grid.y0 - self.grid.dy/2 + self.grid.ny*np.random.random() * self.grid.dy
                    if int(Path(self.polygon.polygon).contains_point((x, y))):
                        rand_x.append(x)
                        rand_y.append(y)
                        validated_inlets += 1
                self.points[points_key] = list(zip(rand_x, rand_y))

        # Case when geological_IDs is indicated
        else:
            if self.geology is None:
                print('- generate_points() - Error : no geology to consider in order to generate points.')
                sys.exit()
            if self.geology.data["geology"]["stats"] is None:
                self.geology.compute_stats_on_data("geology")
            # Keep only present ids
            ids = self.geology.data["geology"]["stats"]["ID"]
            new_ids = [id for id in geological_IDs if (id in ids)]
            if not new_ids:
                 print('- generate_points() - Error : none of the ids provided are valid.')
                 sys.exit()
            validated_inlets = 0
            rand_x = []
            rand_y = []
            # Without polygon
            if self.polygon.polygon is None:
                while validated_inlets < points_number:
                    x = self.grid.x0 - self.grid.dx/2 + self.grid.nx*np.random.random() * self.grid.dx
                    y = self.grid.y0 - self.grid.dy/2 + self.grid.ny*np.random.random() * self.grid.dy
                    if self.geology.data["geology"]["data"][self.grid.get_i(x)][self.grid.get_j(y)][0] in new_ids:
                        rand_x.append(x)
                        rand_y.append(y)
                        validated_inlets += 1
                self.points[points_key] = list(zip(rand_x,rand_y))
            # With polygon
            else:
                while validated_inlets < points_number:
                    x = self.grid.x0 - self.grid.dx/2 + self.grid.nx*np.random.random() * self.grid.dx
                    y = self.grid.y0 - self.grid.dy/2 + self.grid.ny*np.random.random() * self.grid.dy
                    if self.geology.data["geology"]["data"][self.grid.get_i(x)][self.grid.get_j(y)][0] in new_ids:
                        if int(Path(self.polygon.polygon).contains_point((x,y))):
                            rand_x.append(x)
                            rand_y.append(y)
                            validated_inlets += 1
                self.points[points_key] = list(zip(rand_x,rand_y))
        return None

    def composite_points(self, points_key, points, points_number):
        """
        Generate random points on the grid, and add indicated points.

        Parameters
        ----------
        points_key : string
            Type of points : 'inlets' or 'outlets'.
        points : string or list
            Location of the datafile or list of the points coordinates.
        points_number : integer
            Number of points to generate.
        """
        self.generate_points(points_key, points_number)
        if isinstance(points, str):
            text = opendatafile(points)
            other_points = loadpoints(text)
            self.points[points_key] += other_points
        else:
            self.points[points_key] += points
        return None

    def inspect_points(self):
        """
        Check if the points are well located inside the grid.
        If there is no print out, so everything is OK.
        """
        if self.points is None:
            print('- inspect_points() - Error : no points to inspect.')
            sys.exit()
        else:
            xlimits = [self.grid.xlimits[0], self.grid.xlimits[0], self.grid.xlimits[1], self.grid.xlimits[1]]
            ylimits = [self.grid.ylimits[0], self.grid.ylimits[1], self.grid.ylimits[1], self.grid.ylimits[0]]
            for key in self.points:
                mask = []
                unvalidated_points = []
                for k,(x,y) in enumerate(self.points[key]):
                    if self.polygon.polygon is not None:
                        a = not int(Path(self.polygon.polygon).contains_point((x, y)))
                    else:
                        a = not int(Path(list(zip(xlimits, ylimits))).contains_point((x, y)))
                    mask.append((a,a))
                    if a:
                        unvalidated_points.append(k+1)
                self.points[key] = ma.masked_array(self.points[key], mask=mask)
                if len(unvalidated_points) > 0:
                    print('- inspect_points() - Warning : {} {} on {} masked because not inside polygon.'.format(len(unvalidated_points),key,len(self.points[key])))
                    for point in unvalidated_points:
                        print('- point on line {}'.format(point))
                if len(unvalidated_points) == len(self.points[key]):
                    print('- inspect_points() - Error : all the "{}" are outside the domain.'.format(key))
                    sys.exit()
            return None

    def show(self):
        """
        Show the delimitation of the grid, of the polygon (if present) and of the locations of the points (if present).
        """
        fig, ax = plt.subplots()
        fig.suptitle('Show points', fontsize=16)
        # Geology
        if self.geology is not None:
            origin = None
            if 'img' in self.geology.data['geology']:
                d = self.geology.data['geology']['img']
                if self.polygon.mask is not None:
                    import numpy.ma as ma
                    mask = np.transpose(self.polygon.mask)
                    d = ma.MaskedArray(d, mask=mask)
            else:
                d = self.geology.data['geology']['data']
            if self.geology.data['geology']['mode'] in ['gslib', 'csv']:
                origin="lower"
            if self.polygon.mask is not None:
                import numpy.ma as ma
                d = ma.MaskedArray(d, mask=self.polygon.mask)
                #geology = np.flipud(np.transpose(geology, (1,0,2))) # we need to reverse transformations from geologymanager
            plt.imshow(d , extent=self.grid.extent, cmap='gray_r')#, origin=origin)

        # Grid limits
        xlimits = [self.grid.xlimits[0], self.grid.xlimits[0], self.grid.xlimits[1], self.grid.xlimits[1], self.grid.xlimits[0]]
        ylimits = [self.grid.ylimits[0], self.grid.ylimits[1], self.grid.ylimits[1], self.grid.ylimits[0], self.grid.ylimits[0]]
        ax.plot(xlimits, ylimits, color='red', label='grid limits')

        # Polygon
        if self.polygon.polygon is not None:
            closed_polygon = self.polygon.polygon[:]
            closed_polygon.append(closed_polygon[0])
            x, y = zip(*closed_polygon)
            ax.plot(x, y, color='black', label='polygon')

        # Points
        for key in self.points:
            x, y = zip(*self.points[key])
            ax.plot(x, y, 'o', label=key)

        plt.xlim((self.grid.xlimits[0] - self.grid.dx/2, self.grid.xlimits[1] + self.grid.dx/2))
        plt.ylim((self.grid.ylimits[0] - self.grid.dy/2, self.grid.ylimits[1] + self.grid.dy/2))
        ax.set_aspect('equal', 'box')
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.show()
        return None
