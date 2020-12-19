"""
pyKasso
=======

pyKasso is an object-oriented library intented for karst network simulation.

License
-------
Released under the MIT license:
   Copyright (C) 2020 Univertsity of Neuchâtel - CHYN
   François Miville <francois.miville@unine.ch>
"""

import os
import sys
import math

# Dependencies
import yaml
import mpmath
import numpy    as np
import numpy.ma as ma
import pandas   as pd
import skfmm
import karstnet as kn
import matplotlib.pyplot as plt
from matplotlib.path    import Path
from matplotlib.patches import Rectangle

import sympy as sp
import shapely.geometry as sg

import pyvista   as pv
import pyvistaqt as pvqt0

#####
# 1 #
###############################
# Create a work environnement #
###############################

def get_settings(example=False):
    """
    Provide the datafiles settings.

    Parameters
    ----------
    example : bool, optionnal
        If True, pyKasso will provide you the Betteraz's files example

    Examples
    --------
        >>> pk.get_settings(example=True)
    """

    # copying defaults file from source package
    import shutil
    import glob

    path = os.path.dirname(os.path.abspath(__file__)) + '/' + 'default_files' + '/'

    if example == True:

        # create inputs directory
        dst = 'inputs'

        if not os.path.exists(dst):
            try:
                os.makedirs(dst)
                print("Directory named '", dst, "' has been created.")
            except FileExistsError:
                raise

        path = path + 'betteraz' + '/'

        srcs = glob.glob(path + '*')
        dst  = dst + '/'

        for src in srcs:
            shutil.copy(src, dst)

    else:
        src = path + 'settings.yaml'
        dst  = 'settings.yaml'

        shutil.copy(src, dst)

    return None

def _opendatafile(file_location):
    """
    Private function to simply open a file.
    """
    try:
        file = open(file_location,"r")
        text = file.readlines()
        file.close()
    except:
        print("Error : impossible to open datafile.")
        raise
    return text

def _loadpoints(text):
    """
    Private function to load points in a text. Generally combined with '_opendatafile()'.
    """
    points = []
    try:
        for data_line in text:
            x, y = data_line.strip().split()
            points.append((float(x), float(y)))
    except ValueError:
        print("Dimensions of data frame does not match with number of variables. Missing values or values in excess.")
        raise
    except:
        print("Unexpected error:", sys.exc_info()[0])
        raise
    return points

def _extents(f):
    """
    Little private function to correctly locate data on plots.
    """
    delta = f[1] - f[0]
    return [f[0] - delta/2, f[-1] + delta/2]

#####
# 2 #
###################################
# pyKasso class for data modeling #
###################################

###########
class Grid():
    """
    Create a grid for the simulation.
    The calculation points are located in the center of the cells, and the origin coordinates are located at the left-bottom of the grid.

    Parameters
    ----------
    x0 : float
        x-coordinate origin.
    y0 : float
        y-coordinate origin.
    z0 : float
        z-coordinate origin.
    xnum : integer
        Number of cells on x dimension.
    ynum : integer
        Number of cells on y dimension.
    znum : integer
        Number of cells on z dimension.
    dx : float
        Width of a cell.
    dy : float
        Length of a cell.
    dz : float
        Height of a cell.

    Notes
    -----
    - On this version, dx, dy and dz must be similar. pyKasso does not support unstructured grid yet.
    - x0, y0 and z0 are considered to be in the center of the first node.
	"""

    def __init__(self, x0, y0, z0, xnum, ynum, znum, dx, dy, dz):
        self.x0   = x0
        self.y0   = y0
        self.z0   = z0
        self.xnum = xnum
        self.ynum = ynum
        self.znum = znum
        self.dx   = dx
        self.dy   = dy
        self.dz   = dz

        # Calculating the grid
        self.x = np.arange(self.x0,self.x0+(self.xnum)*self.dx,self.dx, dtype=np.float_)
        self.y = np.arange(self.y0,self.y0+(self.ynum)*self.dy,self.dy, dtype=np.float_)
        self.z = np.arange(self.z0,self.z0+(self.znum)*self.dz,self.dz, dtype=np.float_)
        self.X,self.Y,self.Z = np.meshgrid(self.x, self.y, self.z)

        # Calculating the limits
        self.xlimits = [self.x0-self.dx/2, self.x[-1]+self.dx/2]
        self.ylimits = [self.y0-self.dy/2, self.y[-1]+self.dy/2]
        self.zlimits = [self.z0-self.dz/2, self.z[-1]+self.dz/2]

        xlimits = [self.x0-self.dx/2, self.x0-self.dx/2,    self.x[-1]+self.dx/2, self.x[-1]+self.dx/2, self.x0-self.dx/2]
        ylimits = [self.y0-self.dy/2, self.y[-1]+self.dy/2, self.y[-1]+self.dy/2, self.y0-self.dy/2,    self.y0-self.dy/2]
        self.limits = list(zip(xlimits, ylimits))


################
class Polygon():
    """
	Create a polygon manager : a class for managing a polygon as an area study delimitor.

	Parameters
	----------
	grid : Grid()
		Polygon() class needs a Grid() object as an argument.
	"""

    def __init__(self, grid):
        self.polygon = None
        self.mask    = None
        self.grid    = grid

    def set_polygon(self, vertices):
        """
        Create a polygon from vertices coordinates.

        Parameters
        ----------
        vertices : string or list
            Location of the datafile or list of the vertices coordinates.
        """
        if isinstance(vertices, str):
            text     = _opendatafile(vertices)
            vertices = _loadpoints(text)

        if len(vertices) < 3:
            print(" set_polygon() - Error : Not enough vertices to create a polygon (3 minimum).")
            return None
        else:
            self.polygon = vertices
            self.mask    = self._set_mask()
        return None

    def _set_mask(self):
        """
        Set the mask.
        """
        mask = np.zeros((self.grid.ynum, self.grid.xnum))
        for y in range(self.grid.ynum):
            for x in range(self.grid.xnum):
                mask[y][x] = not int(Path(self.polygon).contains_point((self.grid.X[y][x][0],self.grid.Y[y][x][0])))
        return mask

    def inspect_polygon(self):
        """
        Check if the vertices of the polygon are located inside the grid or not.
        A polygon must be set.

        Returns
        -------
            List of vertices out of the grid.
        """
        if self.polygon is not None:
            unvalidated_vertices = []
            for k,(x,y) in enumerate(self.polygon):
                if not int(Path(self.grid.limits).contains_point((x,y))):
                    unvalidated_vertices.append(k+1)
            validated_vertices = len(self.polygon) - len(unvalidated_vertices)

            if len(unvalidated_vertices) == len(self.polygon):
                print('- inspect_polygon() - Warning : 0 vertices inside the grid limits. No polygon is set.'.format(len(unvalidated_vertices),len(self.polygon)))
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

    def _show(self):
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

            ax.plot(self.grid.xlimits,self.grid.ylimits, color='red')
            ax.set_aspect('equal', 'box')
            plt.legend(('polygon','grid limits'), loc='center left', bbox_to_anchor=(1, 0.5))
            plt.show()
            return None
        else:
            print('- show_Polygon() - Error : no polygon to show.')
            return None


#####################
class PointManager():
    """
    Create a point manager : a class for managing points and transform them to inlets or outlets.

    Parameters
    ----------
    grid : Grid()
        PointManager() class needs a Grid() object as an argument.
    polygon : Polygon(), optionnal
        PointManager() class needs a Polygon() object as an argument.
    """

    def __init__(self, grid, polygon, geology):
        self.points  = {}
        self.grid    = grid
        self.polygon = polygon
        self.geology = geology

    def set_points(self, points_key, points):
        """
        Set new points.

        Parameters
        ----------
        points_key : string
            Type of points : 'inlets' or 'outlets'.
        points : string or list
            Location of the datafile or list of the points coordinates.
        """
        if isinstance(points, str):
            text = _opendatafile(points)
            self.points[points_key] = _loadpoints(text)
        else:
            self.points[points_key] = points
        return None

    def generate_points(self, points_key, points_number, geological_IDs=None):
        """
        Generate random points on the grid, according to the parameters.

        Parameters
        ----------
        points_key : string
            Type of points : 'inlets' or 'outlets'.
        points_number : integer
            Number of points to generate.
        """
        get_X = lambda x : int(math.ceil((x - self.grid.x0 - self.grid.dx/2) / self.grid.dx))
        get_Y = lambda y : int(math.ceil((y - self.grid.y0 - self.grid.dy/2) / self.grid.dy))

        if self.geology.data["geology"]["mode"] is "null":
            geological_IDs=None

        if geological_IDs is None:
            if self.polygon.polygon is None:
                rand_x = [self.grid.x0 - self.grid.dx/2 + self.grid.xnum * np.random.random() * self.grid.dx for x in range(points_number)]
                rand_y = [self.grid.y0 - self.grid.dy/2 + self.grid.ynum * np.random.random() * self.grid.dy for y in range(points_number)]
                self.points[points_key] = list(zip(rand_x,rand_y))
            else:
                validated_inlets = 0
                rand_x = []
                rand_y = []
                while validated_inlets < points_number:
                    x = self.grid.x0 - self.grid.dx/2 + self.grid.xnum*np.random.random() * self.grid.dx
                    y = self.grid.y0 - self.grid.dy/2 + self.grid.ynum*np.random.random() * self.grid.dy
                    if int(Path(self.polygon.polygon).contains_point((x,y))):
                        rand_x.append(x)
                        rand_y.append(y)
                        validated_inlets += 1
                self.points[points_key] = list(zip(rand_x,rand_y))
        # Case when geological_IDs is indicated
        else:
            validated_inlets = 0
            rand_x = []
            rand_y = []
            # Without polygon
            if self.polygon.polygon is None:
                while validated_inlets < points_number:
                    x = self.grid.x0 - self.grid.dx/2 + self.grid.xnum*np.random.random() * self.grid.dx
                    y = self.grid.y0 - self.grid.dy/2 + self.grid.ynum*np.random.random() * self.grid.dy
                    if self.geology.surface_IDs[get_Y(y)][get_X(x)] in geological_IDs:
                        rand_x.append(x)
                        rand_y.append(y)
                        validated_inlets += 1
                self.points[points_key] = list(zip(rand_x,rand_y))
            # With polygon
            else:
                while validated_inlets < points_number:
                    x = self.grid.x0 - self.grid.dx/2 + self.grid.xnum*np.random.random() * self.grid.dx
                    y = self.grid.y0 - self.grid.dy/2 + self.grid.ynum*np.random.random() * self.grid.dy
                    if (self.geology.surface_IDs[get_Y(y)][get_X(x)] in geological_IDs) and int(Path(self.polygon.polygon).contains_point((x,y))):
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
            text = _opendatafile(points)
            other_points = _loadpoints(text)
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
            for key in self.points:
                mask = []
                unvalidated_points = []
                for k,(x,y) in enumerate(self.points[key]):
                    if self.polygon.polygon is not None:
                        a = not int(Path(self.polygon.polygon).contains_point((x,y)))
                    else:
                        a = not int(Path(self.grid.limits).contains_point((x,y)))
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

    def _show(self):
        """
        Show the delimitation of the grid, of the polygon (if present) and of the locations of the points (if present).
        """
        fig, ax = plt.subplots()
        fig.suptitle('Show points', fontsize=16)

        # Grid limits
        x,y = list(zip(*Path(self.grid.limits).to_polygons()[0]))
        ax.plot(x,y,color='red',label='grid limits')

        # Polygon
        if self.polygon.polygon is not None:
            closed_polygon = self.polygon.polygon[:]
            closed_polygon.append(closed_polygon[0])
            x,y = zip(*closed_polygon)
            ax.plot(x,y,color='black',label='basin')

        # Points
        for key in self.points:
            x,y = zip(*self.points[key])
            ax.plot(x,y,'o',label=key)
        ax.set_aspect('equal', 'box')
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.show()
        return None


#######################
class GeologyManager():
    """
    Create a geology manager : a class for managing geology, faults and fractures.

    Parameters
    ----------
    grid : Grid()
        GeologyManager() class needs a Grid() object as an argument.
    """

    def __init__(self, grid):
        self.data  = {}   #data model = {key : {data : var, stats : var, mode : var}}
        self.grid  = grid

    def set_data_null(self, data_key):
        """
        Set data to 'null' for a indicated data key.

        Parameters
        ----------
        data_key : string
            Type of data : 'geology', 'faults' or 'fractures' or 'karst'.
        """
        self.data[data_key] = {}
        self.data[data_key]['data']  = np.zeros((self.grid.znum, self.grid.ynum, self.grid.xnum))
        occurence = self.grid.znum * self.grid.ynum * self.grid.xnum
        volume    = occurence * self.grid.dx * self.grid.dy * self.grid.dz
        self.data[data_key]['stats'] = {'ID':0, 'occurence':occurence, 'frequency':1, 'volume':volume}
        self.data[data_key]['mode']  = 'null'
        return None

    def set_data(self, data_key, datafile_location):
        """
        Set data from a gslib datafile.
        The datafile must have one column.
        For 'faults', 'fractures' and 'karst', in the gslib datafile :
        0 - absence
        1 - presence

        Parameters
        ----------
        data_key : string
            Type of data : 'geology', 'faults', 'fractures' or 'karst'.
        datafile_location : string
            Path of the gslib datafile.
        """
        self.data[data_key]          = {}
        self.data[data_key]['data']  = self._fill(datafile_location)
        self.data[data_key]['stats'] = self._compute_stats(data_key)
        self.data[data_key]['mode']  = 'import'
        return None

    def _fill(self, datafile_location):
        """
        GSLIB Reader.
        """
        # Try to open datafile
        text = _opendatafile(datafile_location)

        # Control if second line is well an integer
        try:
            nvar = int(text[1])
        except:
            print("- set_data() - Error : Second line of gslib file must be an integer.")
            raise

        # Control if second line is 1
        if nvar is not 1:
           sys.exit("- set_data() - Error : gslib file must have only one column.")

        # Control if size declared match with gslib's size file
        data_lines = text[2+nvar:]
        if len(data_lines) != (self.grid.xnum*self.grid.ynum*self.grid.znum):
            sys.exit("- set_data() - Error : Dimensions declared does not match with gslib file's size.")

        # Get values from text files
        data = np.zeros(len(data_lines))
        try:
            for data_line, k in zip(data_lines, range(len(data_lines))):
                data[k] = data_line.strip().split()[0]
        except:
            print("- set_data() - Unexpected error:", sys.exc_info()[0])
            raise

        return data.reshape((self.grid.znum, self.grid.ynum, self.grid.xnum))

    def _compute_stats(self, data_key):
        """
        Compute statistics on the imported data.
        """
        stats = {}

        # Count occurencies
        for z in range(self.grid.znum):
            for y in range(self.grid.ynum):
                for x in range(self.grid.xnum):
                    value = self.data[data_key]['data'][z][y][x]
                    try:
                        if stats.get(value) == None:
                            stats[value] = 1
                        else:
                            stats[value] += 1
                    except:
                        print(value)

        # Construct stats
        ID        = []
        occurence = []
        frequency = []
        volume    = []
        for k in stats:
            ID.append(k)
            occurence.append(stats[k])
            frequency.append(100*stats[k]/(self.grid.xnum*self.grid.ynum*self.grid.znum))
            volume.append(stats[k]*self.grid.dx*self.grid.dy*self.grid.dz)
        return {'ID':ID, 'occurence':occurence, 'frequency':frequency, 'volume':volume}

    def _check_geological_IDs(self, IDs):
        for ID in IDs:
            if ID not in self.data["geology"]["stats"]["ID"]:
                sys.exit("Detected geological IDs do not match which indicated ones.")
        self.geological_IDs = IDs
        return None

    def _compute_surface(self):
        """
        Compute the surface IDs and the digital elevation model.
        """
        if self.data["geology"]["mode"] is "null":
            self.DEM         = np.full((self.grid.ynum, self.grid.xnum), self.grid.znum * self.grid.dz)
            self.surface_IDs = np.full((self.grid.ynum, self.grid.xnum), 0)
        elif self.data["geology"]["mode"] is "import":
            self.DEM         = np.zeros((self.grid.ynum,self.grid.xnum))
            self.surface_IDs = np.zeros((self.grid.ynum,self.grid.xnum))
            for y in range(self.grid.ynum):
                for x in range(self.grid.xnum):
                    for z in range(self.grid.znum):
                        if self.data["geology"]["data"][z][y][x] == self.geological_IDs[0]:
                            if z == 0:
                                self.surface_IDs[y][x] = self.geological_IDs[0]
                                break
                        else:
                            self.DEM[y][x] = z + 1
                            self.surface_IDs[y][x] = self.data["geology"]["data"][z][y][x]
            self.DEM = self.DEM * self.grid.dz + self.grid.z0
        return None

    def generate_fractures(self, fractures_densities, fractures_min_orientation, fractures_max_orientation, alpha, fractures_min_length, fractures_max_length):
        """
        Generate fractures.

        Parameters
        ----------
        fractures_densities : list
            Fracture densities for each fracture family.
        fractures_min_orientation : list
            Fractures minimum orientation for each fracture family.
        fractures_max_orientation : list
            Fractures maximum orientation for each fracture family.
        alpha : float
            Degree of power law.
        fractures_min_length : list
            The minimum lenght of the fractures. For all the families.
        fractures_max_length : list
            The maximum lenght of the fractures. For all the families.
        """

        self.fractures = []
        fracture_id = 0

        ## Redefine fracturation domain
        lenmax = max(fractures_max_length)
        xd0, xd1, yd0, yd1, zd0, zd1 = self.grid.xlimits[0], self.grid.xlimits[1], self.grid.ylimits[0], self.grid.ylimits[1], self.grid.zlimits[0], self.grid.zlimits[1]
        Lx, Ly, Lz = xd1-xd0, yd1-yd0, zd1-zd0

        shiftx = min(Lx/2,lenmax/2)
        shifty = min(Ly/2,lenmax/2)
        shiftz = min(Lz/2,lenmax/2)

        Lex = 2 * shiftx + Lx
        Ley = 2 * shifty + Ly
        Lez = 2 * shiftz + Lz

        area = Lex * Ley
        xmin = self.grid.xlimits[0] - shiftx
        ymin = self.grid.ylimits[0] - shifty
        zmin = self.grid.zlimits[0] - shiftz


        ## Total numbers of fractures for each family
        fractures_numbers = np.array(fractures_densities) * area


        ## Calculate fractures in each family
        for frac_family in range(len(fractures_densities)):

            # Define all the constants required to randomly draw the length in distribution
            palpha = (1-alpha[frac_family])
            invpalpha = 1/palpha
            fmina = fractures_min_length[frac_family]**palpha
            frangea = fractures_max_length[frac_family]**palpha - fmina

            # Define all the constants required for the orientation distribution
            min_angle = fractures_min_orientation[frac_family]
            max_angle = fractures_max_orientation[frac_family]
            if min_angle > max_angle:
                min_angle = min_angle - 360

            mean_angle = (max_angle + min_angle)  / 2
            std_angle  = (max_angle - mean_angle) / 3
            mean_angle = math.radians(mean_angle)
            std_angle  = math.radians(std_angle)

            # Computes the kappa value for the Von Mises orientation distribution
            func  = lambda k : std_angle**2 - 1 + mpmath.besseli(1,k) / mpmath.besseli(0,k)
            kappa = mpmath.findroot(func,1)

            # Generate poisson number for each fracture family
            real_frac_number = np.random.poisson(fractures_numbers[frac_family])

            # Loop over the individual fractures
            for i in range(1, real_frac_number+1):

                # FRACTURE CENTER LOCATION from triple uniform distribution
                xm = xmin + np.random.random() * Lex
                ym = ymin + np.random.random() * Ley
                zm = zmin + np.random.random() * Lez

                # FRACTURE LENGHT from truncated power law distribution
                u = np.random.rand()
                frac_length = ( fmina + u * frangea )**invpalpha
                radius = frac_length/2

                # FRACTURE ORIENTATION from von Mises distribution
                frac_orientation =  np.random.vonmises(mean_angle, kappa)

                # Calculate normal Vector
                x = np.sin(frac_orientation)
                y = np.cos(frac_orientation)
                a = -y
                b = x
                c = 0
                n = [a,b,c]

                # Store calculated fracture
                self.fractures.append(Fracture(fracture_id, frac_family, [xm,ym,zm], radius, frac_orientation, n))
                fracture_id += 1

        return None

    def discretize_fractures(self):
        """
        Discretize generated fractures.
        """
        self.intersection_points = self.intersect_fractures()
        return None

    def intersect_fractures(self):
        """
        Discretize generated fractures.
        """
        self.data['fractures'] = {}
        self.data['fractures']['data'] = np.zeros((self.grid.ynum, self.grid.xnum, self.grid.znum))

        for fracture in self.fractures:
            intersection_points = []
            x, y, z = fracture.get_position()
            a, b, c = fracture.get_normal()
            fracture_plane = sp.Plane((x, y, z), normal_vector=(a, b, c))

            x0 = self.grid.x[0]
            y0 = self.grid.y[0]
            for z in self.grid.z:
                z_plane = sp.Plane((x0, y0, z), normal_vector=(0, 0, 1))
                intersection_line = fracture_plane.intersection(z_plane)
                if not intersection_line:
                    continue
                else :
                    intersection_line = intersection_line[0]
                    x1, y1, z1 = intersection_line.points[0]
                    x2, y2, z2 = intersection_line.points[1]
                    x1 = round(x1,2)
                    y1 = round(y1,2)
                    x2 = round(x2,2)
                    y2 = round(y2,2)
                    a = (y1-y2)/(x1-x2)
                    b = y1 - a*x1

                    x1 = self.grid.x[0]
                    y1 = a*x1 + b
                    x2 = self.grid.x[-1]
                    y2 = a*x2 + b

                    intersection_line = sg.LineString([(x1,y1), (x2,y2)])

                for x in self.grid.x:
                    #x_line = sp.Line3D(sp.Point3D(x, y0, z), sp.Point3D(x, y0+1, z))
                    x_line = sg.LineString([(x, y0), (x, self.grid.y[-1])])
                    intersection_point = intersection_line.intersection(x_line)
                    if not hasattr(intersection_point, 'y'):
                        continue
                    else:
                        y = intersection_point.y
                    intersection_points.append([x,round(y,2),z])

        return intersection_points



    #def generate_fractures(self, fractures_densities, fractures_min_orientation, fractures_max_orientation, alpha, fractures_min_length, fractures_max_length):
        #if self.grid.znum == 1:
        #    self.generate_2D_fractures(fractures_densities, fractures_min_orientation, fractures_max_orientation, alpha, fractures_min_length, fractures_max_length)
        #else:
        #    self.generate_3D_fractures()

    def generate_2D_fractures(self, fractures_densities, fractures_min_orientation, fractures_max_orientation, alpha, fractures_min_length, fractures_max_length):
        """
        Generate fractures.

        Parameters
        ----------
        fractures_densities : list
            Fracture densities for each fracture family.
        fractures_min_orientation : list
            Fractures minimum orientation for each fracture family.
        fractures_max_orientation : list
            Fractures maximum orientation for each fracture family.
        alpha : float
            Degree of power law.
        fractures_min_length : list
            The minimum lenght of the fractures. For all the families.
        fractures_max_length : list
            The maximum lenght of the fractures. For all the families.
        """
        self.data['fractures'] = {}
        self.data['fractures']['data'] = np.zeros((self.grid.ynum, self.grid.xnum, self.grid.znum))

        self.fractures = {}
        fracture_id = 0

        # Redefine fracturation domain
        lenmax = max( fractures_max_length )

        xd0, xd1, yd0, yd1 = self.grid.xlimits[0], self.grid.xlimits[1], self.grid.ylimits[0], self.grid.ylimits[1]
        Lx, Ly = xd1-xd0, yd1-yd0

        shiftx = min(Lx/2,lenmax/2)
        shifty = min(Ly/2,lenmax/2)
        Lex = 2 * shiftx + Lx
        Ley = 2 * shifty + Ly

        area = Lex * Ley
        xmin = self.grid.xlimits[0] - shiftx
        ymin = self.grid.ylimits[0] - shifty

        # Numbers of fractures
        fractures_numbers = np.array(fractures_densities) * area

        for frac_family in range(len(fractures_densities)):
            fracs = []

            # Define all the constants required to randomly draw the length in distribution
            palpha = (1-alpha[frac_family])
            invpalpha = 1/palpha
            fmina = fractures_min_length[frac_family]**palpha
            frangea = fractures_max_length[frac_family]**palpha - fmina

            # Define all the constants required for the orientation distribution
            min_angle = fractures_min_orientation[frac_family]
            max_angle = fractures_max_orientation[frac_family]
            if min_angle > max_angle:
                min_angle = min_angle - 360

            mean_angle = (max_angle + min_angle)  / 2
            std_angle  = (max_angle - mean_angle) / 3
            mean_angle = math.radians(mean_angle)
            std_angle  = math.radians(std_angle)

            # Computes the kappa value for the Von Mises orientation distribution
            func  = lambda k : std_angle**2 - 1 + mpmath.besseli(1,k) / mpmath.besseli(0,k)
            kappa = mpmath.findroot(func,1)

            # Generate poisson number for each fracture family
            real_frac_number = np.random.poisson(fractures_numbers[frac_family])

            # Loop over the individual fractures
            for i in range(1, real_frac_number+1):

                # FRACTURE CENTER LOCATION from double uniform distribution
                xm = xmin + np.random.random() * Lex
                ym = ymin + np.random.random() * Ley

                # FRACTURE LENGHT from truncated power law distribution
                u = np.random.rand()
                frac_length = ( fmina + u * frangea )**invpalpha
                mid_length = frac_length / 2

                # FRACTURE ORIENTATION from von Mises distribution
                frac_orientation =  np.random.vonmises(mean_angle, kappa)
                sinor = np.sin( frac_orientation )
                cosor = np.cos( frac_orientation )

                # Extremities
                dx1 = mid_length * sinor
                dy1 = mid_length * cosor
                dx2, dy2 = -dx1, -dy1

                x1, y1, x2, y2 = xm+dx1, ym+dy1, xm+dx2, ym+dy2

                # This part allows to crop all fractures generated outside the domain
                # We swap the point to ensure that point 1 as the lowest x
                if x2 < x1 :
                    x1, x2 = x2, x1
                    y1, y2 = y2, y1

                # Slope of fracture
                a = (y1-y2)/(x1-x2)

                # We initialize a flag to detect if fracture is whithin the domain
                flagin = True

                # Truncation along the x axis
                if x1>xd1 or x2<xd0 :
                    flagin = False
                else:
                    if x1<xd0 :
                        x1, y1 = xd0, a*(xd0-x1) + y1
                    if x2 > xd1 :
                        x2, y2 = xd1, a*(xd1-x1) + y1

                # Vertical swap: we want y2>y1
                if y2 < y1 :
                    x1, x2 = x2, x1
                    y1, y2 = y2, y1

                # Truncation along the y axis
                if y2<yd0 or y1 > yd1 :
                    flagin = False
                else :
                    if y1 < yd0 :
                        x1, y1 = (yd0-y1)/a + x1, yd0
                    if y2 > yd1 :
                        x2, y2 = (yd1-y1)/a + x1, yd1

                # We keep only the fractures within the domain
                if flagin:
                    fracs.append(Fracture(fracture_id, frac_family, [x1, y1, x2, y2], frac_length, frac_orientation))
                    fracture_id += 1

            # Complete fractures dictionary
            self.fractures[frac_family] = {'frac_nbr' : real_frac_number, 'fractures' : fracs}

        # Draw fractures maps
        self._generate_fractures_maps()
        return None

    def _generate_fractures_maps(self):
        """
        Draw the fractures map according to generated fracturation.
        """
        fractures_maps = np.zeros((len(self.fractures)+1,self.grid.ynum,self.grid.xnum))

        for frac_family in self.fractures:
            for fracture in self.fractures[frac_family]["fractures"]:
                x0, y0, x1, y1 = fracture.position
                d = np.sqrt((x1-x0)**2 + (y1-y0)**2)
                D = np.sqrt((self.grid.dx/2)**2 + (self.grid.dy/2)**2)
                n = int(d / D)
                xs = np.linspace(x0, x1, n)
                ys = np.linspace(y0, y1, n)

                for (x,y) in zip(xs, ys):
                    X = int(math.ceil((x - self.grid.x0 - self.grid.dx/2) / self.grid.dx))
                    Y = int(math.ceil((y - self.grid.y0 - self.grid.dy/2) / self.grid.dy))
                    fractures_maps[frac_family][Y][X] = 1

            self.fractures[frac_family]['frac_map'] = fractures_maps[frac_family]
        self.data['fractures']['data'] = (sum([fractures_maps[i] for i in range(len(self.fractures))]) > 0).astype(int)
        self.data['fractures']['mode'] = 'generate'
        return None

    def generate_3D_fractures(self):
        sys.exit("no such function yet")

    def show_fractures(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for frac_family in self.fractures:
            for frac in self.fractures[frac_family]["fractures"]:
                x0 = frac.position[0]
                y0 = frac.position[1]
                x1 = frac.position[2]
                y1 = frac.position[3]
                plt.plot([x0, x1], [y0, y1] ,'orange')
        ax.set_aspect('equal', adjustable='box')
        plt.show()
        return None

    def show(self, frac_family=None):
        """
        Show the geology, faults, fractures and initial karst network maps.
        """
        grid = self._get_pyvista_grid()
        grid.cell_arrays["values"] = self.data["geology"]["data"].flatten(order="F")

        p = pv.Plotter(shape=(2, 2))

        p.subplot(0, 0)
        p.add_text("Geology", font_size=30)
        p.add_mesh(grid)

        p.subplot(0, 1)
        p.add_text("Faults", font_size=30)
        p.add_mesh(pv.Cube(), show_edges=True, color="tan")

        p.subplot(1, 0)
        p.add_text("Fractures", font_size=30)
        sphere = pv.Sphere()
        p.add_mesh(sphere, scalars=sphere.points[:, 2])
        p.add_scalar_bar("Z")
        # plotter.add_axes()
        p.add_axes(interactive=True)

        p.subplot(1, 1)
        p.add_text("Karst network", font_size=30)
        p.add_mesh(pv.Cone(), color="g", show_edges=True)
        p.show_bounds(all_edges=True)

        # Display the window
        p.show()

    def _get_pyvista_grid(self):
        """
        Produce the pyvista grid.
        """
        values = self.data["geology"]["data"]
        grid = pv.UniformGrid()
        # Set the grid dimensions: shape + 1 because we want to inject our values on the CELL data
        grid.dimensions = np.array(values.shape) + 1
        # Edit the spatial reference
        grid.origin  = (self.grid.x0-self.grid.dx/2, self.grid.y0-self.grid.dy/2, self.grid.z0-self.grid.dz/2)
        grid.spacing = (self.grid.dx, self.grid.dy, self.grid.dz)
        return grid

def _extents(f):
    """
    Little private function to correctly locate data on plots.
    """
    delta = f[1] - f[0]
    return [f[0] - delta/2, f[-1] + delta/2]


#################
class Fracture():
    """
    A class for modeling fractures as objects.

    Parameters
    ----------
    ID : integer
        Fracture id
    family : integer
        Fracture family id
    position : list
        [x0, y0, x1, y1]
    radius : float
        Radius of the fracture
    orientation : float
        Orientation of the fracture
    normal : List
        Normal vector [a, b, c]
    """

    def __init__(self, ID, family, position, radius, orientation, normal):

        self.ID          = ID
        self.family      = family
        self.position    = position
        self.radius      = radius
        self.orientation = orientation
        self.normal      = normal

    def __repr__(self):
        return '[id:{}, fam.:{}, x:{}, y:{}, z:{}, rad:{}, or.:{}] \n'.format(self.ID, self.family, round(self.position[0],2), round(self.position[1],2), round(self.position[2],2), round(self.radius,2), round(self.orientation,2))

    def get_ID(self):
        return self.ID

    def get_family(self):
        return self.family

    def get_position(self):
        return self.position

    def get_radius(self):
        return self.radius

    def get_orientation(self):
        return self.orientation

    def get_normal(self):
        return self.normal

#####
# 3 #
##############################################
# pyKasso class for karst network simulation #
##############################################

class SKS():
    """
    A super-class to manage all the previous data class and to simulate karst networks.
    """

    def __init__(self, yaml_settings_file = None, rand_seed = None):
        """
        Construct a SKS class according to the specified settings datafile.

        Parameters
        ----------
        yaml_settings_file : string
            YAML settings file location.

        Examples
        --------
            >>> catchment = pk.SKS()
            >>> catchment = pk.SKS('exemple.yaml')
        """
        if yaml_settings_file is None:
            yaml_settings_file = os.path.dirname(os.path.abspath(__file__)) + '/' + 'default_files/settings.yaml'

        try:
            with open(yaml_settings_file, 'r') as stream:
                try:
                    settings = yaml.safe_load(stream)
                except yaml.YAMLError as exc:
                    print(exc)
                    raise
        except:
            print("/!\\ Error : unable to read the YAML settings file!")
            raise

        self.settings = settings

        self.settings['fractures_numbers'] = [int(self.settings['xnum']*self.settings['ynum']*self.settings['dx']*self.settings['dy']*float(i)) for i in self.settings['fractures_densities']]

        geology_velocity = []
        for elem in self.settings['geology_velocity']:
            geology_velocity.append(self.settings[elem])
        self.settings['geology_velocity'] = geology_velocity

        if self.settings['rand_seed'] > 0:
            np.random.seed(self.settings['rand_seed'])

        if rand_seed != None:
            np.random.seed(rand_seed)
            self.settings['rand_seed'] = rand_seed
        self._get_reference_statistics()
        self._load_data()
        self.karst_simulations = []

    ################
    # get_ methods #
    ################

    def get_x0(self):
        return self.settings['x0']

    def get_y0(self):
        return self.settings['y0']

    def get_z0(self):
        return self.settings['z0']

    def get_xnum(self):
        return self.settings['xnum']

    def get_ynum(self):
        return self.settings['ynum']

    def get_znum(self):
        return self.settings['znum']

    def get_dx(self):
        return self.settings['dx']

    def get_dy(self):
        return self.settings['dy']

    def get_dz(self):
        return self.settings['dz']

    def get_data_has_polygon(self):
        return self.settings['data_has_polygon']

    def get_polygon_data(self):
        return self.settings['polygon_data']

    def get_inlets_mode(self):
        return self.settings['inlets_mode']

    def get_inlets_data(self):
        return self.settings['inlets_data']

    def get_inlets_number(self):
        return self.settings['inlets_number']

    def get_outlets_mode(self):
        return self.settings['outlets_mode']

    def get_outlets_data(self):
        return self.settings['outlets_data']

    def get_outlets_number(self):
        return self.settings['outlets_number']

    def get_geological_mode(self):
        return self.settings['geological_mode']

    def get_geological_datafile(self):
        return self.settings['geological_datafile']

    def get_faults_mode(self):
        return self.settings['faults_mode']

    def get_faults_datafile(self):
        return self.settings['faults_datafile']

    def get_fractures_mode(self):
        return self.settings['fractures_mode']

    def get_fractures_datafile(self):
        return self.settings['fractures_datafile']

    def get_fractures_densities(self):
        return self.settings['fractures_densities']

    def get_fractures_min_orientation(self):
        return self.settings['fractures_min_orientation']

    def get_fractures_max_orientation(self):
        return self.settings['fractures_max_orientation']

    def get_alpha(self):
        return self.settings['alpha']

    def get_fractures_min_length(self):
        return self.settings['fractures_min_length']

    def get_fractures_max_length(self):
        return self.settings['fractures_max_length']

    def get_code_out(self):
        return self.settings['code_out']

    def get_code_aquifere(self):
        return self.settings['code_aquifere']

    def get_code_aquiclude(self):
        return self.settings['code_aquiclude']

    def get_code_faults(self):
        return self.settings['code_faults']

    def get_code_fractures(self):
        return self.settings['code_fractures']

    def get_code_conduits(self):
        return self.settings['code_conduits']

    def get_geology_id(self):
        return self.settings['geology_id']

    def get_geology_velocity(self):
        return self.settings['geology_velocity']

    def get_importance_factor(self):
        return self.settings['importance_factor']

    def get_rand_seed(self):
        return self.settings['rand_seed']

    ##################################
    # get_ methods for uploaded data #
    ##################################

    def get_polygon_vertices(self):
        """
        Get the polygon vertices as a list.
        """
        if self.polygon.polygon is None:
            print("No polygon set.")
            return None
        else:
            return self.polygon.polygon

    def get_inlets(self):
        """
        Get the inlets as a list.
        """
        return self.inlets

    def get_outlets(self):
        """
        Get the outlets as a list.
        """
        return self.outlets

    def get_geology(self):
        """
        Get the geological data as a numpy-array.
        """
        return self.geology.data['geology']['data']

    def get_faults(self):
        """
        Get the faults data as a numpy-array.
        """
        return self.geology.data['faults']['data']

    def get_fractures(self, fractures_family=None):
        """
        Get the fractures as a numpy-array.

        Parameter
        ---------
        fractures_family : integer
            First family is 0
        """
        if fractures_family is None:
            return self.geology.data['fractures']['data']
        else:
            return self.geology.fractures[fractures_family]['frac_map']

    def get_fractures_numbers(self):
        """
        Get the number of fractures.

        Return
        ------
        dictionnary
            A dictionnary with the number of fractures given by user, and with the number of fractures calculated.
        """
        user = self.fractures_numbers
        model = []
        for frac_family in self.geology.fractures:
            model.append(self.geology.fractures[frac_family]['frac_nbr'])
        fracs = {'user' : user, 'model' : model}
        return fracs

    def get_mask(self):
        """
        Get the numpy-array mask calculated if a polygon is given.
        """
        if self.mask is None:
            return 'No mask to return.'
        else:
            return self.mask

    ################
    # set_ methods #
    ################

    def set_data_has_polygon(self, data_has_polygon):
        """
        Define if study area will be focused on a polygon delimitations.

        Parameter
        ---------
        data_has_polygon : bool
            If true, a polygon will be required.
        """
        self.settings['data_has_polygon'] = data_has_polygon
        return None

    def set_polygon_data(self, polygon_data):
        """
        Define the polygon datafile path.
        Usefull only when data_has_polygon is true.

        Parameter
        ---------
        polygon_data : string or list
            Polygon datafile path or list of vertices coordinates.
        """
        self.settings['polygon_data'] = polygon_data
        return None

    def set_inlets_mode(self, inlets_mode):
        """
        Define the inlets mode.

        Parameter
        ---------
        inlets_mode : string
            'random'    - Full random points
            'import'    - Import points
            'composite' - Add n random points to imported points
        """
        self.settings['inlets_mode'] = inlets_mode
        return None

    def set_inlets_data(self, inlets_data):
        """
        Define the inlets datafile path.
        Usefull only when inlets mode is on 'import' or 'composite'.

        Parameter
        ---------
        inlets_data : string or list
            Inlets datafile path or list of inlets coordinates.
        """
        self.settings['inlets_data'] = inlets_data
        return None

    def set_inlets_number(self, inlets_number):
        """
        Define the number of inlets to generate.
        Usefull only when inlets mode is on 'random' or 'composite'.

        Parameter
        ---------
        inlets_number : string
            Number of inlets to generate.
        """
        self.settings['inlets_number'] = inlets_number
        return None

    def set_outlets_mode(self, outlets_mode):
        """
        Define the outlets mode.

        Parameter
        ---------
        outlets_mode : string
            'random'    - Full random points
            'import'    - Import points
            'composite' - Add n random points to imported points
        """
        self.settings['outlets_mode'] = outlets_mode
        return None

    def set_outlets_data(self, outlets_data):
        """
        Define the outlets datafile path.
        Usefull only when outlets mode is on 'import' or 'composite'.

        Parameter
        ---------
        outlets_data : string or list
            Outlets datafile path or list of outlets coordinates.
        """
        self.settings['outlets_data'] = outlets_data
        return None

    def set_outlets_number(self, outlets_number):
        """
        Define the number of outlets to generate.
        Usefull only when outlets mode is on 'random' or 'composite'.

        Parameter
        ---------
        outlets_number : string
            Number of outlets to generate.
        """
        self.settings['outlets_number'] = outlets_number
        return None

    def set_geological_mode(self, geological_mode):
        """
        Define the geological mode.

        Parameter
        ---------
        geological_mode : string
            'null'   - No geology
            'import' - Import geology
            'image'  - Import geology via image
        """
        self.settings['geological_mode'] = geological_mode
        return None

    def set_geological_datafile(self, geological_datafile):
        """
        Define the geological datafile path.
        Usefull only when geological mode is on 'import' or 'image'.

        Parameter
        ---------
        geological_datafile : string
            Geological datafile path.
        """
        self.settings['geological_datafile'] = geological_datafile
        return None

    def set_faults_mode(self, faults_mode):
        """
        Define the mode for the faults.

        Parameter
        ---------
        faults_mode : string
            'null'   - No geology
            'import' - Import geology
            'image'  - Import geology via image
        """
        self.settings['faults_mode'] = faults_mode
        return None

    def set_faults_datafile(self, faults_datafile):
        """
        Define the faults datafile path.
        Usefull only when the mode for faults is on 'import' or 'image'.

        Parameter
        ---------
        faults_datafile : string
            Faults datafile path.
        """
        self.settings['faults_datafile'] = faults_datafile
        return None

    def set_fractures_mode(self, fractures_mode):
        """
        Define the mode for the fractures.

        Parameter
        ---------
        fractures_mode : string
            'null'   - No geology
            'import' - Import geology
            'image'  - Import geology via image
            'random' - Generate random fractures
        """
        self.settings['fractures_mode'] = fractures_mode
        return None

    def set_fractures_datafile(self, fracture_datafile):
        """
        Define the fractures datafile path.
        Usefull only when the mode for fractures is on 'import' or 'image'.

        Parameter
        ---------
        fracture_datafile : string
            Fractures datafile path.
        """
        self.fracture_datafile = fracture_datafile
        return None

    def set_fractures_densities(self, fractures_densities):
        """
        Define the fractures densitiy for each fractures family.
        Usefull only when the mode for fractures is on 'random'.

        Parameter
        ---------
        fractures_densities : list
            One density for each family.
        """
        self.settings['fractures_densities'] = fractures_densities
        self.settings['fractures_numbers'] = [self.settings['xnum']*self.settings['ynum']*self.settings['dx']*self.settings['dy']*i/10**6 for i in self.settings['fractures_densities']]
        return None

    def set_fractures_min_orientation(self, fractures_min_orientation):
        """
        Define the minimum orientation of the fracturation for each fractures family.
        Usefull only when the mode for fractures is on 'random'.

        Parameter
        ---------
        fractures_min_orientation : list
            One orientation for each family.
        """
        self.settings['fractures_min_orientation'] = fractures_min_orientation
        return None

    def set_fractures_max_orientation(self, fractures_max_orientation):
        """
        Define the maximum orientation of the fracturation for each fractures family.
        Usefull only when the mode for fractures is on 'random'.

        Parameter
        ---------
        fractures_max_orientation : list
            One orientation for each family.
        """
        self.settings['fractures_max_orientation'] = fractures_max_orientation
        return None

    def set_alpha(self, alpha):
        """
        Define the ???.
        Usefull only when the mode for fractures is on 'random'.

        Parameter
        ---------
        alpha : float
        """
        self.settings['alpha'] = alpha
        return None

    def set_fractures_min_length(self, fractures_min_length):
        """
        Define the minimum lenght for all the fractures.
        Usefull only when the mode for fractures is on 'random'.

        Parameter
        ---------
        fractures_min_length : float
        """
        self.settings['fractures_min_length'] = fractures_min_length
        return None

    def set_fractures_max_length(self, fractures_max_length):
        """
        Define the maximum lenght for all the fractures.
        Usefull only when the mode for fractures is on 'random'.

        Parameter
        ---------
        fractures_max_length : float
        """
        self.settings['fractures_max_length'] = fractures_max_length
        return None

    def set_code_out(self, code_out):
        """
        Define the fast-marching value for the outside of the study area.
        The value must be low to avoid some unreallistic conduits.

        Parameter
        ---------
        code_out : float
        """
        self.settings['code_out'] = code_out
        return None

    def set_code_aquifere(self, code_aquifere):
        """
        Define the fast-marching value for the aquifere cells.

        Parameter
        ---------
        code_aquifere : float
        """
        self.settings['code_aquifere'] = code_aquifere
        return None

    def set_code_aquiclude(self, code_aquiclude):
        """
        Define the fast-marching value for the aquiclude cells.

        Parameter
        ---------
        code_aquiclude : float
        """
        self.settings['code_aquiclude'] = code_aquiclude
        return None

    def set_code_faults(self, code_faults):
        """
        Define the fast-marching value for the faults cells.

        Parameter
        ---------
        code_faults : float
        """
        self.settings['code_faults'] = code_faults
        return None

    def set_code_fractures(self, code_fractures):
        """
        Define the fast-marching value for the fractures cells.

        Parameter
        ---------
        code_fractures : float
        """
        self.settings['code_fractures'] = code_fractures
        return None

    def set_code_conduits(self, code_conduits):
        """
        Define the fast-marching value for the conduits cells.

        Parameter
        ---------
        code_conduits : float
        """
        self.settings['code_conduits'] = code_conduits
        return None

    def set_geology_id(self, geology_id):
        """
        Define the geology id (from geology datafile) to consider in the simulation.
        Only needed in 'import' mode.

        Parameter
        ---------
        geology_id : list
        """
        self.settings['geology_id'] = geology_id
        return None

    def set_geology_velocity(self, geology_velocity):
        """
        Define the velocities to apply to each geology id.

        Parameter
        ---------
        geology_velocity : list
        """
        self.settings['geology_velocity'] = geology_velocity
        return None

    def set_importance_factor(self, importance_factor):
        """
        Define the importance factor, and so the number of karstic conduits generation.

        Parameter
        ---------
        importance_factor : list
        """
        self.settings['importance_factor'] = importance_factor
        return None

    def set_rand_seed(self, rand_seed):
        """
        Define the random seed.
        May help for reproduicity.

        Parameter
        ---------
        rand_seed : integer
        """
        self.settings['rand_seed'] = rand_seed
        if self.settings['rand_seed'] == 0:
            np.random.seed()
        else:
            np.random.seed(self.settings['rand_seed'])
        return None

    def increment_rand_seed(self):
        """
        Increment by one the value of the random seed.
        """
        self.settings['rand_seed'] += 1
        np.random.seed(self.settings['rand_seed'])
        return None

    ###################
    # update_ methods #
    ###################

    def update_polygon(self):
        """
        Update the polygon settings.
        """
        if self.settings['data_has_polygon']:
            self.polygon.set_polygon(self.settings['polygon_data'])
            self.mask = self.polygon.mask
            self.polygon.inspect_polygon()
        else:
            self.polygon.polygon = None
            self.polygon.mask = None
            self.mask = self.polygon.mask
        return None

    def update_inlets(self):
        """
        Update the inlets settings.
        """
        if self.settings['inlets_mode'] == 'random':
            self.points.generate_points('inlets', self.settings['inlets_number'])
        elif self.settings['inlets_mode'] == 'import':
            self.points.set_points('inlets', self.settings['inlets_data'])
        elif self.settings['inlets_mode'] == 'composite':
            self.points.composite_points('inlets', self.settings['inlets_data'], self.settings['inlets_number'])
        else:
            print('/!\\ Error : unrecognized inlets setting.')
            sys.exit()
        self.points.inspect_points()
        self.inlets = self.points.points['inlets'][:]
        return None

    def update_outlets(self):
        """
        Update the outlets settings.
        """
        if self.settings['outlets_mode'] == 'random':
            self.points.generate_points('outlets', self.settings['outlets_number'])
        elif self.settings['outlets_mode'] == 'import':
            self.points.set_points('outlets', self.settings['outlets_data'])
        elif self.settings['outlets_mode'] == 'composite':
            self.points.composite_points('outlets', self.settings['outlets_data'], self.settings['outlets_number'])
        else:
            print('/!\\ Error : unrecognized outlets setting.')
            sys.exit()
        self.points.inspect_points()
        self.outlets = self.points.points['outlets'][:]
        return None

    def update_geology(self):
        """
        Update the geology settings.
        """
        if self.settings['geological_mode'] == 'null':
            self.geology.set_data_null('geology')
        elif self.settings['geological_mode'] == 'import':
            self.geology.set_data('geology', self.settings['geological_datafile'])
        elif self.settings['geological_mode'] == 'image':
            self.geology.set_data_from_image('geology', self.settings['geological_datafile'])
        else:
            print('/!\\ Error : unrecognized geological mode.')
            sys.exit()
        self.geology.compute_stats_on_data()
        if self.settings['data_has_polygon']:
            self.geology_masked['geology'] = ma.MaskedArray(self.geology.data['geology']['data'], mask=self.mask)
        return None

    def update_faults(self):
        """
        Update the faults settings.
        """
        if self.settings['faults_mode'] == 'null':
            self.geology.set_data_null('faults')
        elif self.settings['faults_mode'] == 'import':
            self.geology.set_data('faults', self.settings['faults_datafile'])
        elif self.settings['faults_mode'] == 'image':
            self.geology.set_data_from_image('faults', self.settings['faults_datafile'])
        else:
            print('/!\\ Error : unrecognized faults mode.')
            sys.exit()
        self.geology.compute_stats_on_data()
        if self.settings['data_has_polygon']:
            self.geology_masked['faults'] = ma.MaskedArray(self.geology.data['faults']['data'], mask=self.mask)
        return None

    def update_fractures(self):
        """
        Update the fractures settings.
        """
        if self.settings['fractures_mode'] == 'null':
            self.geology.set_data_null('fractures')
        elif self.settings['fractures_mode'] == 'import':
            self.geology.set_data('fractures', self.settings['fractures_datafile'])
        elif self.settings['fractures_mode'] == 'image':
            self.geology.set_data_from_image('fractures', self.settings['fractures_datafile'])
        elif self.settings['fractures_mode'] == 'random':
            self.geology.generate_fractures(self.settings['fractures_numbers'], self.settings['fractures_min_orientation'], self.settings['fractures_max_orientation'], self.settings['alpha'], self.settings['fractures_min_length'], self.settings['fractures_max_length'])
        else:
            print('/!\\ Error : unrecognized fractures mode.')
            sys.exit()
        self.geology.compute_stats_on_data()
        if self.settings['data_has_polygon']:
            self.geology_masked['fractures'] = ma.MaskedArray(self.geology.data['fractures']['data'], mask=self.mask)
        return None

    def update_all(self):
        """
        Update all the parameters.
        """
        self.update_polygon()
        self.update_inlets()
        self.update_outlets()
        self.update_geology()
        self.update_faults()
        self.update_fractures()
        return None

    def shuffle_inlets(self):
        """
        Shuffle the inlets order.
        """
        np.random.shuffle(self.inlets)
        return None

    ############################
    # initialization functions #
    ############################

    # __init__()
    def _get_reference_statistics(self):
        """
        Get the reference statistics for comparing it with karstnet outputs.
        """
        path = os.path.dirname(os.path.abspath(__file__)) + '/' + 'default_files/statistics.xlsx'
        dataframe = pd.read_excel(path, usecols = "B:H,K")
        self.reference_statistics = {}
        for key in dataframe:
            min_val = np.nanmin(dataframe[key])
            max_val = np.nanmax(dataframe[key])
            self.reference_statistics[key] = (min_val,max_val)
        return None

    # __init__()
    def _load_data(self):
        """
        Initialize all the data.
        """
        # The grid
        self.grid    = self._set_grid(
            self.settings['x0'], self.settings['y0'], self.settings['z0'],
            self.settings['xnum'], self.settings['ynum'], self.settings['znum'],
            self.settings['dx'], self.settings['dy'], self.settings['dz']
        )
        # The polygon and his mask
        self.polygon = self._set_polygon(self.grid)
        self.mask    = self.polygon.mask
        # The points
        self.points  = self._set_points(self.grid, self.polygon)
        self.inlets  = self.points.points['inlets'][:]
        self.outlets = self.points.points['outlets'][:]
        np.random.shuffle(self.inlets)
        # The geologic data
        self.geology = self._set_geology(self.grid)
        self.geology.compute_stats_on_data()

        self.geology_masked = {}
        if self.mask is not None:
            for key in self.geology.data:
                self.geology_masked[key] = ma.MaskedArray(self.geology.data[key]['data'], mask=self.mask)

        return None

    def _set_grid(self, x0, y0, z0, xnum, ynum, znum, dx, dy, dz):
        """
        Set the grid object.
        """
        grid = Grid(x0, y0, z0, xnum, ynum, znum, dx, dy, dz)
        return grid

    def _set_polygon(self, grid):
        """
        Set the polygon object.
        """
        polygon = Polygon(grid)
        if self.settings['data_has_polygon']:
            polygon.set_polygon(self.settings['polygon_data'])
            polygon.inspect_polygon()
        return polygon

    def _set_points(self, grid, polygon):
        """
        Set the point manager object.
        """
        points = PointManager(grid, polygon)
        ###Inlets###
        if self.settings['inlets_mode'] == 'random':
            points.generate_points('inlets', self.settings['inlets_number'])
        elif self.settings['inlets_mode'] == 'import':
            points.set_points('inlets', self.settings['inlets_data'])
        elif self.settings['inlets_mode'] == 'composite':
            points.composite_points('inlets', self.settings['inlets_data'], self.settings['inlets_number'])
        else:
            print('/!\\ Error : unrecognized inlets setting.')
            sys.exit()
        ###Outlets###
        if self.settings['outlets_mode'] == 'random':
            points.generate_points('outlets', self.settings['outlets_number'])
        elif self.settings['outlets_mode'] == 'import':
            points.set_points('outlets', self.settings['outlets_data'])
        elif self.settings['outlets_mode'] == 'composite':
            points.composite_points('outlets', self.settings['outlets_data'], self.settings['outlets_number'])
        else:
            print('/!\\ Error : unrecognized outlets setting.')
            sys.exit()
        points.inspect_points()
        return points

    def _set_geology(self, grid):
        """
        Set the geology manager object.
        """
        geology = GeologyManager(grid)
        # Geology
        if self.settings['geological_mode'] == 'null':
            geology.set_data_null('geology')
        elif self.settings['geological_mode'] == 'import':
            geology.set_data('geology', self.settings['geological_datafile'])
        elif self.settings['geological_mode'] == 'image':
            geology.set_data_from_image('geology', self.settings['geological_datafile'])
        else:
            print('/!\\ Error : unrecognized geological mode.')
            sys.exit()
        # Faults
        if self.settings['faults_mode'] == 'null':
            geology.set_data_null('faults')
        elif self.settings['faults_mode'] == 'import':
            geology.set_data('faults', self.settings['faults_datafile'])
        elif self.settings['faults_mode'] == 'image':
            geology.set_data_from_image('faults', self.settings['faults_datafile'])
        else:
            print('/!\\ Error : unrecognized faults mode.')
            sys.exit()
        # Fractures
        if self.settings['fractures_mode'] == 'null':
            geology.set_data_null('fractures')
        elif self.settings['fractures_mode'] == 'import':
            geology.set_data('fractures', self.settings['fractures_datafile'])
        elif self.settings['fractures_mode'] == 'image':
            geology.set_data_from_image('fractures', self.settings['fractures_datafile'])
        elif self.settings['fractures_mode'] == 'random':
            geology.generate_fractures(self.settings['fractures_numbers'], self.settings['fractures_min_orientation'], self.settings['fractures_max_orientation'], self.settings['alpha'], self.settings['fractures_min_length'], self.settings['fractures_max_length'])
        else:
            print('/!\\ Error : unrecognized fractures mode.')
            sys.exit()
        return geology

    ##############################
    # karst simulation functions #
    ##############################

    def compute_karst_network(self, verbose = False):
        """
        Compute the karst network according to the parameters.

        Save the results in the `karst_simulations` list attribute.
        """
        # 1 - Initialize the parameters
        self._initialize_karst_network_parameters()

        # 2 - Compute conduits for each generation
        self._compute_iterations_karst_network()

        # 3 - Gather the conduits points to construct a node network
        edges, nodes = self._compute_nodes_network()

        # 4 - Calculate the karst network statistics indicators with karstnet and save karst network
        k = kn.KGraph(edges, nodes)
        stats = k.characterize_graph(verbose)

        maps = self.maps.copy()

        points = {}
        points['inlets']  = self.inlets
        points['outlets'] = self.outlets

        network = {}
        network['edges'] = edges
        network['nodes'] = nodes

        config = self.settings

        self.karst_simulations.append(KarstNetwork(maps, points, network, stats, config))
        return None

    # 1
    def _initialize_karst_network_parameters(self):
        """
        Initialize the karst network parameters.
        """
        self.nbr_iteration = len(self.settings['importance_factor'])
        self.inlets        = self._set_inlets_repartition()

        # Raster maps
        self.maps = {}
        self.maps['phi']      = np.ones((self.grid.ynum, self.grid.xnum))
        self.maps['velocity'] = np.zeros((self.nbr_iteration, self.grid.ynum, self.grid.xnum))
        self.maps['time']     = np.zeros((self.nbr_iteration, self.grid.ynum, self.grid.xnum))
        self.maps['karst']    = np.zeros((self.nbr_iteration, self.grid.ynum, self.grid.xnum))

        # Vector maps
        self.nodeID      = 0
        self.conduits    = []
        self.outletsNode = []
        return None

    # 1
    def _set_inlets_repartition(self):
        """
        Distribute the inlets points between the different simulating generation according to the `importance_factor` setting.
        """
        inlets_repartition = []
        inlets_proportion  = []
        total = float(sum(self.settings['importance_factor']))

        for iteration in range(self.nbr_iteration):
            inlets_proportion.append(float(self.settings['importance_factor'][iteration])/total)
            inlets_repartition.append(round(inlets_proportion[iteration]*len(self.inlets)))
        inlets_repartition[-1] = len(self.inlets)-sum(inlets_repartition[0:-1])

        inlets = []
        i = 0
        for k,repartition_nbr in enumerate(inlets_repartition):
            for num in range(repartition_nbr):
                inlets.append((self.inlets[i][0],self.inlets[i][1],k))
                i += 1
        return inlets

    # 2
    def _compute_iterations_karst_network(self):
        """
        Compute each generation of karst conduits.
        """
        # Define phi map according to outlets emplacements
        #print('-START-')
        self._compute_phi_map()

        # Compute velocity map and travel time for each iteration and draw network
        for iteration in range(self.nbr_iteration):
            self._compute_velocity_map(iteration)
            self._compute_time_map(iteration)
            self._compute_karst_map(iteration)
#            print('iteration:{}/{}'.format(iteration+1,self.nbr_iteration))
#        print('- END -')
        return None

    # 2
    def _compute_phi_map(self):
        """
        Compute the phi map.
        """
        for (x,y) in self.outlets:
            X = int(math.ceil((x - self.grid.x0 - self.grid.dx/2) / self.grid.dx))
            Y = int(math.ceil((y - self.grid.y0 - self.grid.dy/2) / self.grid.dy))
            self.maps['phi'][Y][X] = -1
        return None

    # 2
    def _compute_velocity_map(self, iteration):
        """
        Compute the velocity map.
        """
        # If it's the first iteration, iniatialize the velocity map according to the geological settings.
        if iteration == 0:
            # Geology
            if self.geology.data['geology']['mode'] == 'null':
                self.maps['velocity'][0] = np.full((self.grid.ynum, self.grid.xnum), self.settings['code_aquifere'])
            elif self.geology.data['geology']['mode'] == 'image':
                self.maps['velocity'][0] = np.where(self.geology.data['geology']['data']==1, self.settings['code_aquiclude'], self.settings['code_aquifere'])
            elif self.geology.data['geology']['mode'] == 'import':

                tableFMM = {}
                if len(self.settings['geology_id']) != len(self.settings['geology_velocity']):
                    print("- _compute_velocity_map() - Error : number of lithologies does not match with number of FMM code.")
                    sys.exit()

                for geology, codeFMM in zip(self.settings['geology_id'], self.settings['geology_velocity']):
                    if geology in self.geology.data['geology']['stats']['ID']:
                        tableFMM[geology] = codeFMM
                    else:
                        print("- initialize_velocityMap() - Warning : no geology n {} found.".format(geology))
                        tableFMM[geology] = self.settings['code_out']

                for y in range(self.grid.ynum):
                    for x in range(self.grid.xnum):
                        geology = self.geology.data['geology']['data'][y][x]
                        self.maps['velocity'][0][y][x] = tableFMM[geology]

            # Faults
            self.maps['velocity'][0] = np.where(self.geology.data['faults']['data'] > 0, self.settings['code_faults'], self.maps['velocity'][0])

            # Fractures
            self.maps['velocity'][0] = np.where(self.geology.data['fractures']['data'] > 0, self.settings['code_fractures'], self.maps['velocity'][0])

            # If out of polygon
            if self.mask is not None:
                self.maps['velocity'][0] = np.where(self.mask==1, self.settings['code_out'], self.maps['velocity'][0])

        else:
            self.maps['velocity'][iteration] = self.maps['velocity'][iteration-1]
            self.maps['velocity'][iteration] = np.where(self.maps['karst'][iteration-1] > 0, self.settings['code_conduits'], self.maps['velocity'][iteration])
        return None

    # 2
    def _compute_time_map(self, iteration):
        """
        Compute the time map.
        """
        try:
            self.maps['time'][iteration] = skfmm.travel_time(self.maps['phi'], self.maps['velocity'][iteration], dx=self.grid.dx, order=2)
        except:
            try:
                self.maps['time'][iteration] = skfmm.travel_time(self.maps['phi'], self.maps['velocity'][iteration], dx=self.grid.dx, order=1)
            except:
                raise
        return None

    # 2
    def _compute_karst_map(self, iteration):
        """
        Compute the karst map.
        """
        factor = 2
        self.step = self.grid.dx/factor

        grad_y, grad_x = np.gradient(self.maps['time'][iteration], self.grid.dx, self.grid.dy)

        get_X = lambda x : int(math.ceil((x - self.grid.x0 - self.grid.dx/2) / self.grid.dx))
        get_Y = lambda y : int(math.ceil((y - self.grid.y0 - self.grid.dy/2) / self.grid.dy))

        # get karst map from previous iteration
        if iteration > 0:
            self.maps['karst'][iteration] = self.maps['karst'][iteration-1]

        for (x,y,i) in self.inlets:
            # check if inlet iteration match with actual iteration
            if i == iteration:
                conduit = Conduit(iteration)
                X = get_X(x)
                Y = get_Y(y)
                while self.maps['time'][iteration][Y][X] != np.amin(self.maps['time'][iteration]):

                    # If X,Y is on an old karstic channel, then stop
                    if self.maps['karst'][iteration-1][Y][X] == 1:
                        conduit.add_node(self.nodeID,x,y,1) # node type 1 (conjonction)
                        self.nodeID += 1
                        break
                    else:
                        self.maps['karst'][iteration][Y][X] = 1

                    # If X,Y is on an outlet, then stop
                    #if iteration == 0 and self.phiMap[Y][X] == -1:
                    if self.maps['phi'][Y][X] == -1:
                        self.maps['karst'][iteration][Y][X] = 1
                        conduit.add_node(self.nodeID,x,y,2) # node type 2 (end)
                        self.nodeID += 1
                        break

                    # else conduit is here
                    conduit.add_node(self.nodeID,x,y,0)
                    self.nodeID += 1

                    # new move
                    alpha = math.sqrt((self.step**2)*(1/(grad_x[Y][X]**2+grad_y[Y][X]**2)))
                    dx = grad_x[Y][X] * alpha
                    dy = grad_y[Y][X] * alpha

                    # check if we are going out boundaries
                    X_ = get_X(x - dx)
                    Y_ = get_Y(y - dy)
                    if (X_ < 0) or (Y_ < 0) or (X_ > self.grid.xnum - 1) or (Y_ > self.grid.ynum - 1):
                        dx_,dy_ = self._check_boundary_conditions(iteration,X,Y)
                        x = x + dx_
                        y = y + dy_
                    else: # otherwise acts normal
                        x = x - dx
                        y = y - dy

                    X = get_X(x)
                    Y = get_Y(y)

                self.conduits.append(conduit)
        return None

    # 2
    def _check_boundary_conditions(self,iteration,X,Y):
        borders_conditions = [     #(X,Y)
            [(0,-1),(0,1),(1,0) ], # X = 0
            [(-1,0),(1,0),(0,1) ], # Y = 0
            [(0,-1),(0,1),(-1,0)], # X = xnum
            [(-1,0),(1,0),(0,-1)]] # Y = ynum

        corners_conditions = [ #(X,Y)
            [(1,0) ,(0,1) ],   # X = 0    && Y = 0
            [(1,0) ,(0,-1)],   # X = 0    && Y = ynum
            [(-1,0),(0,1) ],   # X = xnum && Y = ynum
            [(-1,0),(0,-1)]]   # X = xnum && Y = 0

        time_values = []
        rank = 0
        if (X <= 0) and (Y <= 0):
            for row in corners_conditions[0]:
                time_values.append(self.maps['time'][iteration][Y+row[1]][X+row[0]])
            rank  = time_values.index(min(time_values))
            moove = corners_conditions[0][rank]
        elif (X <= 0) and (Y >= self.grid.ynum-1):
            for row in corners_conditions[1]:
                time_values.append(self.maps['time'][iteration][Y+row[1]][X+row[0]])
            rank  = time_values.index(min(time_values))
            moove = corners_conditions[1][rank]
        elif (X >= self.grid.xnum-1) and (Y >= self.grid.ynum-1):
            for row in corners_conditions[2]:
                time_values.append(self.maps['time'][iteration][Y+row[1]][X+row[0]])
            rank  = time_values.index(min(time_values))
            moove = corners_conditions[2][rank]
        elif (X >= self.grid.xnum-1) and (Y <= 0):
            for row in corners_conditions[3]:
                time_values.append(self.maps['time'][iteration][Y+row[1]][X+row[0]])
            rank  = time_values.index(min(time_values))
            moove = corners_conditions[3][rank]
        elif X <= 0:
            for row in borders_conditions[0]:
                time_values.append(self.maps['time'][iteration][Y+row[1]][X+row[0]])
            rank  = time_values.index(min(time_values))
            moove = borders_conditions[0][rank]
        elif Y <= 0:
            for row in borders_conditions[1]:
                time_values.append(self.maps['time'][iteration][Y+row[1]][X+row[0]])
            rank  = time_values.index(min(time_values))
            moove = borders_conditions[1][rank]
        elif X >= self.grid.xnum-1:
            for row in borders_conditions[2]:
                time_values.append(self.maps['time'][iteration][Y+row[1]][X+row[0]])
            rank  = time_values.index(min(time_values))
            moove = borders_conditions[2][rank]
        elif Y >= self.grid.ynum-1:
            for row in borders_conditions[3]:
                time_values.append(self.maps['time'][iteration][Y+row[1]][X+row[0]])
            rank  = time_values.index(min(time_values))
            moove = borders_conditions[3][rank]

        else:
            print(X,Y)
        dx = moove[0] * self.step
        dy = moove[1] * self.step

        return (dx,dy)

    # 3
    def _compute_nodes_network(self):
        """
        This script computes the nodes network with all the previous calculated conduits.
        """
        # Add outlets
        outlets = self.points.points['outlets']
        for outlet in outlets:
            self.outletsNode.append((self.nodeID,outlet[0],outlet[1]))
            self.nodeID += 1

        # Connect regular nodes
        [conduit.compute_edges() for conduit in self.conduits]

        for conduit in self.conduits:
            dist = []
            last_node = conduit.nodes[-1]
            x = last_node[1]
            y = last_node[2]
            # connect last node on old karstic conduit
            if last_node[-1] == 1:
                older_conduits = [conduit_ for conduit_ in self.conduits if conduit_.iteration < conduit.iteration]
                all_nodes = []
                for older_conduit in older_conduits:
                    for node in older_conduit.nodes:
                        all_nodes.append(node)
                for node in all_nodes:
                    dx   = abs(x - node[1])
                    dy   = abs(y - node[2])
                    dist.append((math.sqrt(dx**2 + dy**2),node[0]))
                min_node = min(dist)
                conduit.edges.append((last_node[0],min_node[1]))

            # connect last node on outlet coordinates
            if last_node[-1] == 2:
                for outlet in self.outletsNode:
                    dx = abs(x - outlet[1])
                    dy = abs(y - outlet[2])
                    dist.append((math.sqrt(dx**2 + dy**2),outlet[0]))
                min_outlet = min(dist)
                conduit.edges.append((last_node[0],min_outlet[1]))

        # Get data for karstnet
        EDGES = []
        NODES = {}
        for conduit in self.conduits:
            for edge in conduit.edges:
                EDGES.append(edge)
            for node in conduit.nodes:
                NODES[node[0]] = (node[1],node[2])
        for outlet in self.outletsNode:
            NODES[outlet[0]] = (outlet[1],outlet[2])
        return (EDGES,NODES)


    ###########################
    # visualization functions #
    ###########################

    def show_catchment(self, data='geology', title=None, mask=False, cmap='binary'):
        """
        Show the entire study domain.
        """
        fig, ax1 = plt.subplots()
        if title is None:
            title = data
        fig.suptitle(title, fontsize=16)

        if   data == 'geology':
            d = self.geology.data['geology']['data']
        elif data == 'faults':
            d = self.geology.data['faults']['data']
        elif data == 'fractures':
            d = self.geology.data['fractures']['data']

        if mask==True:
            if   data == 'geology':
                d = self.geology_masked['geology']
            elif data == 'faults':
                d = self.geology_masked['faults']
            elif data == 'fractures':
                d = self.geology_masked['fractures']

        im1 = ax1.imshow(d, extent=_extents(self.grid.x) + _extents(self.grid.y), origin='lower', cmap=cmap)
        fig.colorbar(im1, ax=ax1)
        if self.settings['data_has_polygon']:
            closed_polygon = self.polygon.polygon[:]
            closed_polygon.append(closed_polygon[0])
            x,y = zip(*closed_polygon)
            ax1.plot(x,y, color='red', label='polygon')
        for key in self.points.points:
            x,y = zip(*self.points.points[key])
            ax1.plot(x,y,'o',label=key)
        ax1.set_aspect('equal', 'box')
        plt.legend(loc='upper right')
        plt.show()
        return None

    def _show_karst_network(self, sim=-1, iteration=-1, cmap='binary'):
        """
        Show the simulated karst network.
        """
        karst_network = self.karst_simulations[sim]

        fig, ([ax1,ax2],[ax3,ax4]) = plt.subplots(2, 2, sharex=True, sharey=True)
        fig.suptitle('Karst Network', fontsize=16)

        ax1.imshow(karst_network.maps['phi'], extent=_extents(self.grid.x) + _extents(self.grid.y), origin='lower', cmap=cmap)
        ax1.set_title('Phi')

        ax2.imshow(karst_network.maps['velocity'][iteration], extent=_extents(self.grid.x) + _extents(self.grid.y), origin='lower', cmap=cmap)
        ax2.set_title('Velocity')

        ax3.imshow(karst_network.maps['time'][iteration], extent=_extents(self.grid.x) + _extents(self.grid.y), origin='lower', cmap=cmap)
        ax3.set_title('Time')

        ax4.imshow(karst_network.maps['karst'][iteration], extent=_extents(self.grid.x) + _extents(self.grid.y), origin='lower', cmap=cmap)
        ax4.set_title('Karst')

        fig.subplots_adjust(hspace=0.5)
        plt.show()
        return None

    def show(self, data=None, title=None, cmap='binary', probability=False):
        """
        Show the entire study domain.
        """
        if data is None:
            data = self.karst_simulations[-1].maps['karst'][-1]

        if probability == True:
            data = self._compute_average_paths()

        fig, ax1 = plt.subplots()
        if title is not None:
            fig.suptitle(title, fontsize=16)

        im1 = ax1.imshow(data, extent=_extents(self.grid.x) + _extents(self.grid.y), origin='lower', cmap=cmap)

        fig.colorbar(im1, ax=ax1)

        if self.settings['data_has_polygon']:
            closed_polygon = self.polygon.polygon[:]
            closed_polygon.append(closed_polygon[0])
            x,y = zip(*closed_polygon)
            ax1.plot(x,y, color='red', label='polygon')
        for key in self.points.points:
            x,y = zip(*self.points.points[key])
            ax1.plot(x,y,'o',label=key)
        ax1.set_aspect('equal', 'box')
        plt.legend(loc='upper right')
        plt.show()
        return None

    def _compute_average_paths(self, show=False):
        """
        Compute the mean of all the simulations.
        """
        karst_maps = []
        for karst_simulation in self.karst_simulations:
            karst_maps.append(karst_simulation.maps['karst'][-1])

        karst_prob = sum(karst_maps)/len(karst_maps)

        if self.mask is not None:
            karst_prob = np.ma.MaskedArray(karst_prob, self.mask)

        return karst_prob

    def compare_stats(self, data=None, mean=False):
        """
        Compare statistics between reference indicators and calculated networks.
        """

        cpd  = []
        cv_d = []
        cv_l = []
        o_e  = []
        l_e  = []
        spl  = []
        m_d  = []
        cvd  = []
        vars = [cpd,cv_d,cv_l,o_e,l_e,spl,m_d,cvd]

        if data == None:
            i = -1
        else:
            i = 0

        for karst_network in self.karst_simulations[i:]:
            stats_list = []
            stats_list.append(karst_network.stats['cpd'])
            stats_list.append(karst_network.stats['cv degree'])
            stats_list.append(karst_network.stats['cv length'])
            stats_list.append(karst_network.stats['orientation entropy'])
            stats_list.append(karst_network.stats['length entropy'])
            stats_list.append(karst_network.stats['aspl'])
            stats_list.append(karst_network.stats['mean degree'])
            stats_list.append(karst_network.stats['correlation vertex degree'])

            cpd.append(karst_network.stats['cpd'])
            cv_d.append(karst_network.stats['cv degree'])
            cv_l.append(karst_network.stats['cv length'])
            o_e.append(karst_network.stats['orientation entropy'])
            l_e.append(karst_network.stats['length entropy'])
            spl.append(karst_network.stats['aspl'])
            m_d.append(karst_network.stats['mean degree'])
            cvd.append(karst_network.stats['correlation vertex degree'])

        if mean == False:
            print('\n')
            print('STATS for modelisation')
            print('%-10s%-12s%-12s%-12s%-12s' % ('Var', 'Value', 'Min ref', 'Max ref', 'Result'))
            print(54*'-')
            for stat_calc, key in zip(stats_list, self.reference_statistics):
                if stat_calc > self.reference_statistics[key][1]:
                    result = 'out +'
                elif stat_calc < self.reference_statistics[key][0]:
                    result = 'out -'
                else:
                    result = 'IN'
                print('%-10s%-12s%-12s%-12s%-12s' % (key, round(stat_calc,3), round(self.reference_statistics[key][0],3), round(self.reference_statistics[key][1],3), result))

        else:
            print('\n')
            print('MEAN STATS')
            print('%-10s%-12s%-12s%-12s%-12s' % ('Var', 'Value', 'Min ref', 'Max ref', 'Result'))
            print(54*'-')
            for var, key in zip(vars, self.reference_statistics):
                mean = sum(var)/len(var)
                if mean > self.reference_statistics[key][1]:
                    result = 'out +'
                elif mean < self.reference_statistics[key][0]:
                    result = 'out -'
                else:
                    result = 'IN'
                print('%-10s%-12s%-12s%-12s%-12s' % (key, round(mean,3), round(self.reference_statistics[key][0],3), round(self.reference_statistics[key][1],3), result))
        return None


################
class Conduit():
    """
    A class for storing simulated conduits.
    """
    def __init__(self,iteration):
        self.iteration = iteration
        self.nodes     = []
        self.edges     = []

    def add_node(self,nodeID,x,y,nodeType):
        """
        Add a new node in the current constructed conduit.

        Parameters
        ----------
        nodeID:
            ID of the node
        x:
            x position
        y:
            y position
        nodeType:
            0 - normal
            1 - old karstic channel
            2 - outlet
        """
        self.nodes.append((nodeID,x,y,nodeType))
        return None

    def compute_edges(self):
        for i in range(len(self.nodes)):
            if i == 0:
                continue
            node_1 = self.nodes[i-1][0]
            node_2 = self.nodes[i][0]
            self.edges.append((node_1,node_2))
        return None

#####################
class KarstNetwork():
    """
    A class for stroring a calculated karst network.
    """
    def __init__(self, maps, points, network, stats, settings):
        self.maps     = maps
        self.points   = points
        self.network  = network
        self.stats    = stats
        self.settings = settings
