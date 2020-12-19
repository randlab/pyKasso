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
import copy
import yaml # dependance

import math
import mpmath # dependance
import numpy    as np # dependance
import numpy.ma as ma
import pandas   as pd # dependance

import matplotlib
import matplotlib.pyplot as plt
from   matplotlib.path import Path

import karstnet as kn # dependance
import agd # dependance
from agd import Eikonal
from agd.Metrics import Riemann



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
		x-coordinate origin (centerpoint of bottom left cell, NOT bottom left corner of bottom left cell).
	y0 : float
		y-coordinate origin (centerpoint of bottom left cell, NOT bottom left corner of bottom left cell).
	xnum : integer
		Number of cells on x dimension.
	ynum : integer
		Number of cells on y dimension.
	dx : float
		Width of a cell.
	dy : float
		Length of a cell.

	Notes
	-----
	- On this version, dx and dy must be similar. pyKasso does not support unstructured grid yet.  
	- x0 and y0 are considered to be in the center of the first cell.
	"""        
    
    def __init__(self, x0, y0, xnum, ynum, dx, dy):
        self.x0   = x0        #x coord of centerpoint of bottom left cell 
        self.y0   = y0        #y coord of centerpoint of bottom left cell
        self.xnum = xnum      #x resolution: number of cells in x direction (aka columns)
        self.ynum = ynum      #y resolution: number of cells in y direction (aka rows)
        self.dx   = dx        #cell width in x direction
        self.dy   = dy        #cell width in y direction

        self.x        = np.arange(self.x0,self.x0+(self.xnum)*self.dx,self.dx,dtype=np.float_)  #1D array of centerpoints of each cell along x axis
        self.y        = np.arange(self.y0,self.y0+(self.ynum)*self.dy,self.dy,dtype=np.float_)  #1D array of centerpoints of each cell along y axis
        self.X,self.Y = np.meshgrid(self.x,self.y)            #2D array of dim (xnum, ynum) with xy coord of each cell's centerpoint - useful for plotting
        self.xlimits  = [self.x0-self.dx/2,self.x0-self.dx/2,self.x[-1]+self.dx/2,self.x[-1]+self.dx/2,self.x0-self.dx/2] #x coord of outermost cell edges [bottom left, top left, top right, bottom right, bottom left]
        self.ylimits  = [self.y0-self.dy/2,self.y[-1]+self.dy/2,self.y[-1]+self.dy/2,self.y0-self.dy/2,self.y0-self.dy/2] #y coord of outermost cell edges [bottom left, top left, top right, bottom right, bottom left]
        self.limits   = list(zip(self.xlimits,self.ylimits)) #array of tuples with coordinates of corners [(bottom left/origin), (top left), (top right), (bottom right), (bottom left/origin)]
        self.xmin     = self.x0    - self.dx/2      #coordinate of leftmost edge of leftmost cells
        self.xmax     = self.x[-1] + self.dx/2      #coordinate of rightmost edge of rightmost cells
        self.ymin     = self.y0    - self.dy/2      #coordinate of bottom edge of bottom cells
        self.ymax     = self.y[-1] + self.dy/2      #coordinate of top edge of top cells
        self.extent   = [self.xmin, self.xmax, self.ymin, self.ymax]  #coordinates of extent for use in plt.imshow()

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
            text         = _opendatafile(vertices)
            self.polygon = _loadpoints(text)
        else:
            self.polygon = vertices
            
        self.mask = self._set_mask()
        return None
    
    def _set_mask(self):
        """
        Set the mask.
        """
        mask = np.zeros((self.grid.ynum, self.grid.xnum))
        for y in range(self.grid.ynum):
            for x in range(self.grid.xnum):
                mask[y][x] = not int(Path(self.polygon).contains_point((self.grid.X[y][x],self.grid.Y[y][x])))
        return mask

    def inspect_polygon(self):
        """
        Check if the vertices of the polygon are well located inside the grid.
        """
        if self.polygon is not None:
            unvalidated_vertices = []
            for k,(x,y) in enumerate(self.polygon):
                if not int(Path(self.grid.limits).contains_point((x,y))):
                    unvalidated_vertices.append(k+1)
            validated_vertices = len(self.polygon) - len(unvalidated_vertices)
            if validated_vertices < 3:
                print(' inspect_polygon() - Warning : not enough vertices inside the grid limits to proceed further ({} exactly).'.format(validated_vertices))
                print("/!\\ No polygon set /!\\")
                self.polygon = None
                return None
            if len(unvalidated_vertices) > 0:
                print('- inspect_polygon() - Warning : {} vertices not inside the grid limits on {} vertices.'.format(len(unvalidated_vertices),len(self.polygon)))
                for vertex in unvalidated_vertices:
                    print('- vertice {}'.format(vertex))
                return unvalidated_vertices
        else:
            print('- inspect_polygon() - Error : no polygon to inspect. Please set a polygon.')
        return None

    def show(self):
        """
        Show the delimitation of the grid and, if a polygon is present, display its limits.
        """
        if self.polygon is not None:
            closed_polygon = self.polygon[:]
            closed_polygon.append(closed_polygon[0])
            x,y = zip(*closed_polygon)

            fig, ax = plt.subplots()
            fig.suptitle('Polygon', fontsize=16)
            ax.plot(x,y, color='red')

            ax.plot(self.grid.xlimits,self.grid.ylimits, color='black')
            ax.set_aspect('equal', 'box')
            plt.legend(('polygon','grid limits'), loc='center left', bbox_to_anchor=(1, 0.5))
            plt.show()
            return None
        else:
            print('- show_Polygon() - Error : no polygon to show.')
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
    
    def __init__(self, grid, polygon=None):
        self.points  = {}
        self.grid    = grid
        self.polygon = polygon

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

    def generate_points(self, points_key, points_number):
        """
        Generate random points on the grid, according to the parameters.

        Parameters
        ----------
        points_key : string
            Type of points : 'inlets' or 'outlets'.
        points_number : integer
            Number of points to generate.
        """
        if self.polygon.polygon is None:
            rand_x = [self.grid.x0 - self.grid.dx/2 + self.grid.xnum*np.random.random() * self.grid.dx for x in range(points_number)]
            rand_y = [self.grid.y0 - self.grid.dy/2 + self.grid.ynum*np.random.random() * self.grid.dy for y in range(points_number)]
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

    def show(self):
        """
        Show the delimitation of the grid, of the polygon (if present) and of the locations of the points (if present).
        """
        fig, ax = plt.subplots()
        fig.suptitle('Show points', fontsize=16)
        ax.plot(self.grid.xlimits,self.grid.ylimits,color='black',label='grid limits')

        if self.polygon.polygon is not None:
            closed_polygon = self.polygon.polygon[:]
            closed_polygon.append(closed_polygon[0])
            x,y = zip(*closed_polygon)
            ax.plot(x,y,color='red',label='basin')

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
            Type of data : 'geology', 'topography', 'orientation', 'faults' or 'fractures'.
        """
        self.data[data_key] = {}
        self.data[data_key]['data'] = np.zeros((self.grid.ynum, self.grid.xnum)) 
        self.data[data_key]['mode'] = 'null'
        return None

    def set_data(self, data_key, datafile_location, selected_var=0):
        """
        Set data from a datafile for an indicated data key.

        Parameters
        ----------
        data_key : string
            Type of data : 'geology', 'topography', 'orientation', 'faults' or 'fractures'.
        datafile_location : string
            Path of the datafile.
        selected_var : integer, optionnal
            A way to select the column which will be extracted from the datafile.
        """
        self.data[data_key]         = {}
        self.data[data_key]['data'] = self._fill(datafile_location, selected_var)
        self.data[data_key]['mode'] = 'import'
        return None

    def set_data_from_image(self, data_key, datafile_location):
        """
        Set data from an image for an indicated data key.

        Parameters
        ----------
        data_key : string
            Type of data : 'geology', 'topography', 'orientation', 'faults' or 'fractures'.
        datafile_location : string
            Path of the datafile.
        """
        try:
           image_data = np.flipud(plt.imread(datafile_location))
        except:
            print('- set_data_from_image() - Error : unable to read datafile.')
            raise
        self.data[data_key] = {}
        self.data[data_key]['data'] = (image_data[:,:,0] == 0)*1
        self.data[data_key]['mode'] = 'image'
        return None

    def set_data_from_csv(self, data_key, datafile_location):
        """
        Set data from a csv file for indicated data key.

        Parameters
        ----------
        data_key : string
            Type of data : 'geology', 'topography', 'orientation', 'faults' or 'fractures'.
        datafile_location : string
            Path of the datafile.
        """
        self.data[data_key]         = {}
        self.data[data_key]['data'] = np.genfromtxt(datafile_location, delimiter=',')
        self.data[data_key]['mode'] = 'csv'
        return None

    def set_data_from_gslib(self, data_key, datafile_location):
        """
        Set data from a gslib file for indicated data key in the format that work with the fast marching algorithm.

        Parameters
        ----------
        data_key : string
            Type of data : 'geology', 'topography', 'orientation', 'faults' or 'fractures'.
        datafile_location : string
            Path to the datafile.
        """
        self.data[data_key]         = {}
        a = np.genfromtxt(datafile_location, dtype=float, skip_header=3) #read in gslib file as float numpy array without header rows
        a[a==0] = np.nan                                                 #replace zeros with nans (array must be float first)
        a = np.reshape(a, (self.grid.ynum,self.grid.xnum), order='F')    #reshape to xy grid using Fortran ordering
        self.data[data_key]['data'] = a                                  #store  
        self.data[data_key]['mode'] = 'gslib'
        return None
    
    def generate_orientations(self, surface):
        """
        Generate maps of x and y components of orientation.

        Parameters
        ------------
        surface : array
            A 2D array of the elevation or potential in cell, used to calculate the orientation (i.e. slope or dip) of that cell.
            Use either the land surface ('topography'), the surface of the bottom of the karst unit ('contact'), or the potential returned by the geologic model 
        """

        self.data['orientationx'] = {}
        self.data['orientationx']['data'] = np.zeros((self.grid.ynum, self.grid.xnum))
        self.data['orientationy'] = {}
        self.data['orientationy']['data'] = np.zeros((self.grid.ynum, self.grid.xnum))
        self.data['orientationx']['data'], self.data['orientationy']['data'] = np.gradient(surface, self.grid.dx, self.grid.dy, axis=(0,1))   #x and y components of gradient in each cell of array 
        self.data['orientationx']['mode'] = 'topo'
        self.data['orientationy']['mode'] = 'topo'
        return None
    
    def generate_fractures(self, fractures_numbers, fractures_min_orientation, fractures_max_orientation, fractures_alpha, fractures_min_length, fractures_max_length):
        """
        Generate fractures.

        Parameters
        ----------
        fractures_numbers : list
            Fractures number for each fracture family.
        fractures_min_orientation : list
            Fractures minimum orientation for each fracture family.
        fractures_max_orientation : list
            Fractures maximum orientation for each fracture family.
        fractures_alpha : float
            Degree of power law.
        fractures_min_length : float
            The minimum lenght of the fractures. For all the families.
        fractures_max_length : float
            The maximum lenght of the fractures. For all the families.
        """
        self.data['fractures'] = {}
        self.data['fractures']['data'] = np.zeros((self.grid.ynum, self.grid.xnum)) #need to flip

        self.fractures = {}
        fracture_id = 0

        for frac_family in range(len(fractures_numbers)):
            fracs = []
            min_angle = fractures_min_orientation[frac_family]
            max_angle = fractures_max_orientation[frac_family]

            if min_angle > max_angle:
                min_angle = min_angle - 360

            mean_angle = (max_angle + min_angle)  / 2
            std_angle  = (max_angle - mean_angle) / 3
            mean_angle = math.radians(mean_angle)
            std_angle  = math.radians(std_angle)

            func  = lambda k : std_angle**2 - 1 + mpmath.besseli(1,k) / mpmath.besseli(0,k)
            kappa = mpmath.findroot(func,1)

            # Generate poisson number for each fracture family
            real_frac_number = np.random.poisson(fractures_numbers[frac_family])

            for i in range(1, real_frac_number+1):
                # FRACTURE CENTER LOCATION
                # -> from double uniform distribution
                frac_x_start = self.grid.x0 - self.grid.dx/2 + np.random.random() * self.grid.xnum * self.grid.dx
                frac_y_start = self.grid.y0 - self.grid.dy/2 + np.random.random() * self.grid.ynum * self.grid.dy

                # FRACTURE LENGHT
                # -> from power law distribution
                # x   = np.arange(min_fracture_length,max_fracture_length,0.01)
                # C = 1 - np.float_power(max_fracture_length, -alpha+1) / np.float_power(min_fracture_length, -alpha+1)
                # pdf = ((alpha - 1)/min_fracture_length) * np.float_power(x / min_fracture_length, -alpha) / C
                # cdf = (np.float_power(min_fracture_length, -alpha+1) - np.float_power(x, -alpha+1)) / (np.float_power(min_fracture_length, -alpha+1) - np.float_power(max_fracture_length, -alpha+1))
                # cdf_reverse = np.float_power( np.float_power(min_fracture_length, -alpha+1) - P(x) * (np.float_power(min_fracture_length, -alpha+1) - np.float_power(max_fracture_length, -alpha+1)) , (1/(-alpha+1)))
                frac_length = np.float_power(np.float_power(fractures_min_length, -fractures_alpha+1) - np.random.rand() * (np.float_power(fractures_min_length, -fractures_alpha+1) - np.float_power(fractures_max_length, -fractures_alpha+1)) , (1/(-fractures_alpha+1)))

                # FRACTURE ORIENTATION
                # -> from von Mises distribution
                mu = mean_angle
                frac_orientation = math.degrees(np.random.vonmises(mu, kappa))
                if frac_orientation < 0:
                    frac_orientation += 360

                # Create fracture object
                fracs.append(Fracture(fracture_id, frac_family, frac_x_start, frac_y_start, frac_length, frac_orientation))
                fracture_id += 1

            # Complete fractures dictionary
            self.fractures[frac_family] = {'frac_nbr' : real_frac_number, 'fractures' : fracs}

        # Draw fractures maps
        self._generate_fractures_maps()
        return None

    def _generate_fractures_maps(self):
        fractures_maps = np.zeros((len(self.fractures)+1,self.grid.ynum,self.grid.xnum))

        # For each family, compute segment points and draw a map
        for frac_family in self.fractures:
            fractures = self.fractures[frac_family]['fractures']
            [fracture.compute_segment_points(self.grid.dx) for fracture in fractures]
        
            for fracture in fractures:
                x1,y1 = fracture.coordinates[1]
                x2,y2 = fracture.coordinates[2]
                x,y   = (x1,y1)

                if x1 - x2 > 0:
                    fracture.growth_rate_x *= -1
                if y1 - y2 > 0:
                    fracture.growth_rate_y *= -1

                X  = int(math.ceil((x1 - self.grid.x0 - self.grid.dx/2) / self.grid.dx))
                Y  = int(math.ceil((y1 - self.grid.y0 - self.grid.dy/2) / self.grid.dy))
                X2 = int(math.ceil((x2 - self.grid.x0 - self.grid.dx/2) / self.grid.dx))
                Y2 = int(math.ceil((y2 - self.grid.y0 - self.grid.dy/2) / self.grid.dy))

                while X != X2 and Y != Y2:
                    X = int(math.ceil((x - self.grid.x0 - self.grid.dx/2) / self.grid.dx))
                    Y = int(math.ceil((y - self.grid.y0 - self.grid.dy/2) / self.grid.dy))
                    if not (X < 0) and not (Y < 0) and not (X > self.grid.xnum-1) and not (Y > self.grid.ynum-1):
                        fractures_maps[frac_family][Y][X] = 1
                    x = x + fracture.growth_rate_x
                    y = y + fracture.growth_rate_y

            self.fractures[frac_family]['frac_map'] = fractures_maps[frac_family]

        self.data['fractures']['data'] = (sum([fractures_maps[i] for i in range(len(self.fractures))]) > 0).astype(int)
        self.data['fractures']['mode'] = 'generate'
        return None

    def _fill(self, datafile_location, selected_var=0):
        # Try to open datafile
        text = _opendatafile(datafile_location)

        # Check if second line is an integer
        try:
            nvar = int(text[1])
        except:
            print("- set_data() - Error : Second line of gslib's file must be an integer.")
            raise

        # Control if size declared match with gslib's size file
        data_lines = text[2+nvar:]
        if len(data_lines) != (self.grid.xnum*self.grid.ynum):
            print("- set_data() - Error : Dimensions declared does not match with gslib file's dimensions.")
            return None

        # Get variables names
        dataTitles = [var_name.strip() for var_name in text[2:2+nvar]]
        #print("Data variables's name : ", self.dataTitles)
        if selected_var > len(dataTitles):
            print("- set_data() - Error : Selected variable doesn't exist.")
            return None

        # Get values from text files
        data = np.zeros((len(data_lines), nvar))
        try:
            for data_line, k in zip(data_lines, range(len(data_lines))):
                data[k] = data_line.strip().split()
        except ValueError:
            print("- set_data() - Error : Dimensions of data frame does not match with number of variables. Missing values or values in excess.")
            raise
        except:
            print("- set_data() - Unexpected error:", sys.exc_info()[0])
            raise

        # Put values in x,y matrix
        maps = np.zeros((nvar,self.grid.ynum,self.grid.xnum))
        for var in range(nvar):
            k = 0
            for y in range(self.grid.ynum):
                for x in range(self.grid.xnum):
                    maps[var,y,x] = data[k,var]
                    k += 1

        return maps[selected_var]
    
    def compute_stats_on_data(self):
        """
        Compute statistics on the geologic data (not including the orientations).
        """
        #for key in self.data:
        for key in ['geology','faults','fractures']:
            stats = {}
            for y in range(self.grid.ynum):
                for x in range(self.grid.xnum):
                    value = self.data[key]['data'][y][x] 
                    try:
                        if stats.get(value) == None:
                            stats[value] = 1
                        else:
                            stats[value] += 1
                    except:
                        print('Unable to comput stats for value', value)

            #if self.settings['verbosity'] > 1:
                #print('\n')
                #print('STATS for {}'.format(key))
                #print('%-12s%-12s%-12s%-12s' % ('ID', 'Number', '%', 'Superficy'))
                #print(48*'-')
            ID        = []
            occurence = []
            frequency = []
            superficy = []
            for k in stats:
                #if self.settings['verbosity'] > 1:
                   #print('%-12i%-12i%-12f%-12i' % (k, stats[k], 100*stats[k]/(self.grid.xnum*self.grid.ynum), stats[k]*self.grid.dx*self.grid.dy))
                ID.append(k)
                occurence.append(stats[k])
                frequency.append(100*stats[k]/(self.grid.xnum*self.grid.ynum))
                superficy.append(stats[k]*self.grid.dx*self.grid.dy)
            self.data[key]['stats'] = {'ID':ID, 'occurence':occurence, 'frequency':frequency, 'superficy':superficy}
        return None

    def show(self, frac_family=None, cmap='gray_r'):
        """
        Show the geology, faults and fractures maps.
        """
        fig, (ax1,ax2,ax3) = plt.subplots(1,3,sharex=True,sharey=True, figsize=(15,5))
        fig.suptitle('Geologic data', fontsize=16)

        try:
            im1 = ax1.imshow(self.data['geology']['data'], extent=self.grid.extent, origin='lower', cmap=cmap)
        except:
            im1 = []
        ax1.set_title('Geology')

        try:
            im2 = ax2.imshow(self.data['faults']['data'], extent=self.grid.extent, origin='lower', cmap=cmap)
        except:
            im2 = []
        ax2.set_title('Faults')

        try:
            if frac_family is None:
                data = self.data['fractures']['data']
            else:
                data = self.fractures[frac_family]['frac_map']
            im3 = ax3.imshow(data, extent=self.grid.extent, origin='lower', cmap=cmap)
        except:
            im3 = []
        ax3.set_title('Fractures')

        plt.show()
        return None

#def _extents(f):
    """
    Little private function to correctly locate data on plots.
    """
    #delta = f[1] - f[0]
    #return [f[0] - delta/2, f[-1] + delta/2]


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
    x_start : float
        x-coordinate of the center of the fracture
    y_start : flaot
        y-coordinate of the center of the fracture
    length : float
        Length of the fracture
    orientation : float
        Orientation of the fracture
    """
    
    def __init__(self, ID, family, x_start, y_start, length, orientation):

        self.ID           = ID
        self.family       = family
        self.x_start      = x_start
        self.y_start      = y_start
        self.length       = length
        self.orientation  = orientation
        self.coordinates  = []
        self.coordinates.append(np.array((x_start,y_start)))

    def __repr__(self):
        return '[id:{}, x:{}, y:{}, len:{}, or.:{}] \n'.format(self.ID, round(self.x_start,2), round(self.y_start,2), round(self.length,2), round(self.orientation,2))

    def compute_segment_points(self, dx):
        """
        Compute the coordinates of the points of the segment.
        """
        mid_length = self.length / 2
        self.growth_rate_x = abs(math.sin(self.orientation*math.pi/180) * dx) #/ 2
        self.growth_rate_y = abs(math.cos(self.orientation*math.pi/180) * dx) #/ 2
        dx_point1 = mid_length * math.sin(self.orientation*math.pi/180)
        dy_point1 = mid_length * math.cos(self.orientation*math.pi/180)
        dx_point2 = -mid_length * math.sin(self.orientation*math.pi/180)
        dy_point2 = -mid_length * math.cos(self.orientation*math.pi/180)
        self.coordinates.append(self.coordinates[0] + np.array((dx_point1,dy_point1)))
        self.coordinates.append(self.coordinates[0] + np.array((dx_point2,dy_point2)))
        return None


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
        print('CAUTION: You are using the development version of this package.') #uncomment this to tell apart the development version and the main version

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
        
        #Set travel cost through each geologic formation
        geology_cost = []
        for elem in self.settings['geology_cost']:
            geology_cost.append(self.settings[elem])
        self.settings['geology_cost'] = geology_cost
        
        if self.settings['rand_seed'] > 0:
            np.random.seed(self.settings['rand_seed'])
            
        if rand_seed != None:
            np.random.seed(rand_seed)

        self._get_reference_statistics()
        self._load_data()
        self.karst_simulations = []

    ################
    # get_ methods #
    ################

    def get_x0(self):
        print('test')
        return self.settings['x0']

    def get_y0(self):
        return self.settings['y0']

    def get_xnum(self):
        return self.settings['xnum']

    def get_ynum(self):
        return self.settings['ynum']

    def get_dx(self):
        return self.settings['dx']

    def get_dy(self):
        return self.settings['dy']

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

    def get_inlets_shuffle(self):
        return self.settings['inlets_shuffle']

    def get_outlets_mode(self):
        return self.settings['outlets_mode']

    def get_outlets_data(self):
        return self.settings['outlets_data']

    def get_outlets_number(self):
        return self.settings['outlets_number']

    def get_outlets_shuffle(self):
        return self.settings['outlets_shuffle']

    def get_geological_mode(self):
        return self.settings['geological_mode']

    def get_geological_datafile(self):
        return self.settings['geological_datafile']

    def get_topography_mode(self):
        return self.settings['topography_mode']

    def get_topography_datafile(self):
        return self.settings['topography_datafile']

    def get_orientation_mode(self):
        return self.settings['orientation_mode']

    def get_orientation_datafile(self):
        return self.settings['orientation_datafile']

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

    def get_fractures_alpha(self):
        return self.settings['fractures_alpha']

    def get_fractures_min_length(self):
        return self.settings['fractures_min_length']

    def get_fractures_max_length(self):
        return self.settings['fractures_max_length']

    def get_algorithm(self):
        return self.settings['algorithm']

    def get_cost_out(self):
        return self.settings['cost_out']

    def get_cost_aquifer(self):
        return self.settings['cost_aquifer']

    def get_cost_aquiclude(self):
        return self.settings['cost_aquiclude']

    def get_cost_faults(self):
        return self.settings['cost_faults']

    def get_cost_fractures(self):
        return self.settings['cost_fractures']

    def get_cost_conduits(self):
        return self.settings['cost_conduits']

    def get_cost_ratio(self):
        return self.settings['cost_ratio']

    def get_geology_id(self):
        return self.settings['geology_id']
    
    def get_geology_cost(self):
        return self.settings['geology_cost']

    def get_importance_factor(self):
        return self.settings['importance_factor']

    def get_rand_seed(self):
        return self.settings['rand_seed']

    def get_verbosity(self):
        return self.settings['verbosity']

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
        Get the inlets as an array.
        """
        return self.inlets

    def get_outlets(self):
        """
        Get the outlets as an array.
        """
        return self.outlets

    def get_geology(self):
        """
        Get the geological data as a numpy-array.
        """
        return self.geology.data['geology']['data']

    def get_topography(self):
        """
        Get the topography data as a numpy-array.
        """
        return self.geology.data['topography']['data']

    def get_contact_surface(self):
        """
        Get the lower contact surface of the karst unit as a numpy-array.
        """
        return self.geology.data['contact']['data']

    def get_orientation(self):
        """
        Get the orientation data as two numpy arrays (for x and y components).
        """
        return [self.geology.data['orientationx']['data'], self.geology.data['orientationy']['data']]

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
        Useful only when inlets mode is on 'import' or 'composite'.

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
        Useful only when inlets mode is on 'random' or 'composite'.

        Parameter
        ---------
        inlets_number : string
            Number of inlets to generate. 
        """
        self.settings['inlets_number'] = inlets_number
        return None

    def set_inlets_shuffle(self, inlets_shuffle):
        """
        Define whether to shuffle the order of the inlets randomly.
        Useful only when iterating over inlets.

        Parameter
        ---------
        inlets_shuffle : string
            0: don't shuffle, 1: shuffle randomly. 
        """
        self.settings['inlets_shuffle'] = inlets_shuffle
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

    def set_outlets_shuffle(self, outlets_shuffle):
        """
        Define whether to shuffle the order of the outlets randomly.
        Useful only when iterating over outlets.

        Parameter
        ---------
        outlets_shuffle : string
            0: don't shuffle, 1: shuffle randomly. 
        """
        self.settings['outlets_shuffle'] = outlets_shuffle
        return None

    def set_geological_mode(self, geological_mode):
        """
        Define the geological mode.

        Parameter
        ---------
        geological_mode : string
            'null'   - No geology
            'import' - Import from old gslib file format
            'gslib'  - Import from new gslib format
            'csv'    - Import from csv
            'image'  - Import via image
        """
        self.settings['geological_mode'] = geological_mode
        return None

    def set_geological_datafile(self, geological_datafile):
        """
        Define the geological datafile path.
        Usefull only when geological mode is not 'null'.

        Parameter
        ---------
        geological_datafile : string
            Geological datafile path. 
        """
        self.settings['geological_datafile'] = geological_datafile
        return None

    def set_topography_mode(self, topography_mode):
        """
        Define the topography mode.

        Parameter
        ---------
        topography_mode : string
            'null'   - No topography
            'import' - Import from old gslib file format
            'gslib'  - Import from new gslib format
            'csv'    - Import from csv
            'image'  - Import via image
        """
        self.settings['topography_mode'] = topography_mode
        return None

    def set_topography_datafile(self, topography_datafile):
        """
        Define the topography datafile path.
        Usefull only when topography mode is not 'null'.

        Parameter
        ---------
        topography_datafile : string
            topography datafile path. 
        """
        self.settings['topography_datafile'] = topography_datafile
        return None

    def set_orientation_mode(self, orientation_mode):
        """
        Define the orientation mode.

        Parameter
        ---------
        orientation_mode : string
            'null'   - No orientation
            'import' - Import from old gslib file format
            'gslib'  - Import from new gslib format
            'csv'    - Import from csv
            'image'  - Import via image
        """
        self.settings['orientation_mode'] = orientation_mode
        return None

    def set_orientation_datafile(self, orientation_datafile):
        """
        Define the orientation datafile path.
        Usefull only when orientation mode is not 'null'.

        Parameter
        ---------
        orientation_datafile : string
            orientation datafile path. 
        """
        self.settings['orientation_datafile'] = orientation_datafile
        return None

    def set_faults_mode(self, faults_mode):
        """
        Define the mode for the faults.

        Parameter
        ---------
        faults_mode : string
            'null'   - No faults
            'import' - Import from old gslib file format
            'gslib'  - Import from new gslib format
            'csv'    - Import from csv
            'image'  - Import via image
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
            'null'   - No fractures
            'import' - Import from old gslib file format
            'gslib'  - Import from new gslib format
            'csv'    - Import from csv
            'image'  - Import via image
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

    def set_fractures_alpha(self, fractures_alpha):
        """
        Define the ???.
        Usefull only when the mode for fractures is on 'random'.

        Parameter
        ---------
        fractures_alpha : float
        """
        self.settings['fractures_alpha'] = fractures_alpha
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

    def set_cost_out(self, cost_out):
        """
        Define the fast-marching value for the outside of the study area.
        The value must be between 0 and 1 and should be high to avoid unrealistic conduits.

        Parameter
        ---------
        cost_out : float  (default: 0.999)
        """
        self.settings['cost_out'] = cost_out
        return None

    def set_cost_aquifer(self, cost_aquifer):
        """
        Define the fast-marching value for the aquifer cells.
        Should be between 0 and 1 and lower than aquiclude but higher than conduits.

        Parameter
        ---------
        cost_aquifer : float  (default: 0.3)
        """
        self.settings['cost_aquifer'] = cost_aquifer
        return None

    def set_cost_aquiclude(self, cost_aquiclude):
        """
        Define the fast-marching value for the aquiclude cells.
        Should be between 0 and 1 and higher than aquiclude but lower than cost_out

        Parameter
        ---------
        cost_aquiclude : float  (default: 0.8)
        """
        self.settings['cost_aquiclude'] = cost_aquiclude
        return None

    def set_cost_faults(self, cost_faults):
        """
        Define the fast-marching value for the faults cells.
        Should be between 0 and 1 and between conduits and cost_out. Higher = conduit will avoid faults, lower = conduits will follow faults

        Parameter
        ---------
        cost_faults : float  (default: 0.2)
        """
        self.settings['cost_faults'] = cost_faults
        return None

    def set_cost_fractures(self, cost_fractures):
        """
        Define the fast-marching value for the fractures cells.
        Should be between 0 and 1 and between conduits and cost_out. Higher = conduit will avoid fractures, lower = conduits will follow fractures

        Parameter
        ---------
        cost_fractures : float (default: 0.2)
        """
        self.settings['cost_fractures'] = cost_fractures
        return None

    def set_cost_conduits(self, cost_conduits):
        """
        Define the fast-marching value for the conduits cells.
        Should be between 0 and 1 but lower than aquifer (for conduits to preferentially follow each other)

        Parameter
        ---------
        cost_conduits : float  (default: 0.01)
        """
        self.settings['cost_conduits'] = cost_conduits
        return None

    def set_cost_ratio(self, cost_ratio):
        """
        Define the fast-marching ratio of travel cost parallel to gradient / travel cost prependicular to gradient.
        Should be between 0 and 0.5. 

        Parameter
        ---------
        cost_ratio : float  (default: 0.25)
        """
        self.settings['cost_ratio'] = cost_ratio
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
    
    def set_geology_cost(self, geology_cost):
        """
        Define the travel cost to apply to each geology id.
        
        Parameter
        ---------
        geology_cost : list
        """
        self.settings['geology_cost'] = geology_cost
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

    def set_verbosity(self, verbosity):
        """
        Define the verbosity (how much output to print during runs).
        
        Parameter
        ---------
        verbosity: integer 
            0: print minimal output,  1: print medium output,  2: print max output
        """
        self.settings['verbosity'] = verbosity
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
        elif self.settings['geological_mode'] == 'gslib':
            self.geology.set_data_from_gslib('geology', self.settings['geological_datafile'])
        elif self.settings['geological_mode'] == 'csv':
            self.geology.set_data_from_csv('geology', self.settings['geological_datafile'])
        else:
            print('/!\\ Error : unrecognized geological mode', self.settings['geological_mode'])
            sys.exit() 
        self.geology.compute_stats_on_data()
        if self.settings['polygon_data']:
            self.geology_masked['geology'] = ma.MaskedArray(self.geology.data['geology']['data'], mask=self.mask)
        return None

    def update_topography(self):
        """
        Update the topography settings.
        """
        if self.settings['topography_mode'] == 'null':
            self.geology.set_data_null('topography')
        elif self.settings['topography_mode'] == 'import':
            self.geology.set_data('topography', self.settings['topography_datafile'])
        elif self.settings['topography_mode'] == 'image':
            self.geology.set_data_from_image('topography', self.settings['topography_datafile'])
        elif self.settings['topography_mode'] == 'gslib':
            self.geology.set_data_from_gslib('topography', self.settings['topography_datafile'])
        elif self.settings['topography_mode'] == 'csv':
            self.geology.set_data_from_csv('topography', self.settings['topography_datafile'])
        else:
            print('/!\\ Error : unrecognized topography mode', self.settings['topography_mode'])
            sys.exit() 
        self.geology.compute_stats_on_data()
        if self.settings['polygon_data']:
            self.geology_masked['topography'] = ma.MaskedArray(self.geology.data['topography']['data'], mask=self.mask)
        return None

    def update_orientation(self):
        """
        Update the orientation settings.
        """
        if self.settings['orientation_mode'] == 'null':
            self.geology.set_data_null('orientation')
        elif self.settings['orientation_mode'] == 'import':
            self.geology.set_data('orientation', self.settings['orientation_datafile'])
        elif self.settings['orientation_mode'] == 'gslib':
            self.geology.set_data_from_gslib('orientation', self.settings['orientation_datafile'])
        elif self.settings['orientation_mode'] == 'csv':
            self.geology.set_data_from_csv('orientation', self.settings['orientation_datafile'])
        else:
            print('/!\\ Error : unrecognized orientation mode', self.settings['orientation_mode'])
            sys.exit() 
        self.geology.compute_stats_on_data()
        if self.settings['polygon_data']:
            self.geology_masked['orientation'] = ma.MaskedArray(self.geology.data['orientation']['data'], mask=self.mask)
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
        if self.settings['polygon_data']:
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
            self.geology.generate_fractures(self.settings['fractures_numbers'], self.settings['fractures_min_orientation'], self.settings['fractures_max_orientation'], self.settings['fractures_alpha'], self.settings['fractures_min_length'], self.settings['fractures_max_length'])
        else:                                                                                  
            print('/!\\ Error : unrecognized fractures mode.')
            sys.exit()
        self.geology.compute_stats_on_data()
        if self.settings['polygon_data']:
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
        self.update_topography()
        self.update_orientation()
        self.update_faults()
        self.update_fractures()
        return None

    def shuffle_inlets(self):
        """
        Shuffle the inlets order.
        """
        np.random.shuffle(self.inlets)
        return None

    def shuffle_outlets(self):
        """
        Shuffle the outlets order.
        """
        np.random.shuffle(self.outlets)
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
        self.grid    = self._set_grid(self.settings['x0'],self.settings['y0'],self.settings['xnum'],self.settings['ynum'],self.settings['dx'],self.settings['dy'])
        # The polygon and its mask
        self.polygon = self._set_polygon(self.grid)
        self.mask    = self.polygon.mask
        # The points
        self.points  = self._set_points(self.grid, self.polygon)
        self.inlets  = self.points.points['inlets'][:]
        self.outlets = self.points.points['outlets'][:]
        if self.settings['inlets_shuffle'] == True:
            np.random.shuffle(self.inlets)
        if self.settings['outlets_shuffle'] == True:
            np.random.shuffle(self.outlets)
        # The geologic data
        self.geology = self._set_geology(self.grid)
        self.geology.compute_stats_on_data()
        
        self.geology_masked = {}
        if self.mask is not None:
            for key in self.geology.data:
                self.geology_masked[key] = ma.MaskedArray(self.geology.data[key]['data'], mask=self.mask)
        
        return None

    def _set_grid(self, x0, y0, xnum, ynum, dx, dy):
        """
        Set the grid object.
        """
        grid = Grid(x0, y0, xnum, ynum, dx, dy)
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
        elif self.settings['geological_mode'] == 'csv':
            geology.set_data_from_csv('geology', self.settings['geological_datafile'])
        elif self.settings['geological_mode'] == 'gslib':
            geology.set_data_from_gslib('geology', self.settings['geological_datafile'])
        else:
            print('/!\\ Error : unrecognized geological mode', self.settings['geological_mode'])
            sys.exit()
        # Topography
        if self.settings['topography_mode'] == 'null':
            geology.set_data_null('topography')
        elif self.settings['topography_mode'] == 'csv':
            geology.set_data_from_csv('topography', self.settings['topography_datafile'])
        else:
            print('/!\\ Error : unrecognized topography mode', self.settings['topography_mode'])
            sys.exit()
        # Orientation
        if self.settings['orientation_mode'] == 'null':
            geology.set_data_null('orientationx')
            geology.set_data_null('orientationy')
        elif self.settings['orientation_mode'] == 'topo':
            geology.generate_orientations(geology.data['topography']['data'])
        elif self.settings['orientation_mode'] == 'contact':
            geology.set_data_from_csv('contact', self.settings['orientation_datafile'])
            geology.generate_orientations(geology.data['contact']['data'])
        else:
            print('/!\\ Error : unrecognized orientation mode', self.settings['orientation_mode'])
            sys.exit()
        # Faults
        if self.settings['faults_mode'] == 'null':
            geology.set_data_null('faults')
        elif self.settings['faults_mode'] == 'import':
            geology.set_data('faults', self.settings['faults_datafile'])
        elif self.settings['faults_mode'] == 'image':
            geology.set_data_from_image('faults', self.settings['faults_datafile'])
        else:
            print('/!\\ Error : unrecognized faults mode', self.settings['faults_mode'])
            sys.exit()
        # Fractures
        if self.settings['fractures_mode'] == 'null':
            geology.set_data_null('fractures')
        elif self.settings['fractures_mode'] == 'import':
            geology.set_data('fractures', self.settings['fractures_datafile'])
        elif self.settings['fractures_mode'] == 'image':
            geology.set_data_from_image('fractures', self.settings['fractures_datafile'])
        elif self.settings['fractures_mode'] == 'random':
            geology.generate_fractures(self.settings['fractures_numbers'], self.settings['fractures_min_orientation'], self.settings['fractures_max_orientation'], self.settings['fractures_alpha'], self.settings['fractures_min_length'], self.settings['fractures_max_length'])
        else:
            print('/!\\ Error : unrecognized fractures mode', self.settings['fractures_mode'])
            sys.exit()
        return geology
    
    ##############################
    # Karst simulation functions #
    ##############################
    
    def compute_karst_network(self, verbose = False):
        """
        Compute the karst network according to the parameters.
        
        Save the results in the `karst_simulations` list attribute.
        """
        
        # 1 - Initialize the parameters
        self._initialize_karst_network_parameters()
        
        # 2 - Compute conduits for each generation & store nodes and edges for network
        self._compute_iterations_karst_network()
        
        # 3 - Calculate the karst network statistics indicators with karstnet and save karst network
        karstnet_edges = list(self.edges.values()) #convert to format needed by karstnet (list)
        karstnet_nodes = copy.deepcopy(self.nodes)       #convert to format needed by karstnet (dic with only coordinates) - make sure to only modify a copy!
        for key, value in karstnet_nodes.items(): #drop last item in list (the node type) for each dictionary entry
            value.pop()
        k = kn.KGraph(karstnet_edges, karstnet_nodes)  #make graph - edges must be a list, and nodes must be a dic of format {nodeindex: [x,y]}
        stats = k.characterize_graph(verbose)
        
        # 4 - Store all the relevant data for this network in dictionaries:
        maps = self.maps.copy() 
        points = {}
        points['inlets']  = self.inlets
        points['outlets'] = self.outlets
        network = {}
        network['edges'] = self.edges   #store edges list
        network['nodes'] = self.nodes   #store nodes list
        network['karstnet'] = k    #store karstnet network object (including graph)
        config = self.settings
        
        self.karst_simulations.append(KarstNetwork(maps, points, network, stats, config))
        return None
    
    # 1
    def _initialize_karst_network_parameters(self):
        """
        Initialize the karst network parameters.
        """
        #Iterations
        self.nbr_iteration   = len(self.settings['outlets_importance'])*len(self.settings['inlets_importance']) #total number of iterations that will occur
        outlets_repartition  = self._repartition_points(self.outlets, self.settings['outlets_importance']) #correct for outlets_importance not summing to correct number of actual outlets
        self.outlets         = self._distribute_outlets(outlets_repartition)  # distribute outlets to iterations
        inlets_repartition   = self._repartition_points(self.inlets, self.settings['inlets_per_outlet']) #correct for inlets_per_outlet not summing to correct number of actual inlets
        self.inlets          = self._distribute_inlets(inlets_repartition)  # distribute inlets to outlets
        
        inlets  = []
        for o,outlet in enumerate(self.outlets):                 #loop over outlets
            inlets_current = self.inlets[self.inlets[:,2]==o]    #select only inlets assigned to current outlet
            repartition    = self._repartition_points(inlets_current, self.settings['inlets_importance']) #correct for inlets_importance not summing to correct number of actual inlets available for current outlet
            i = 0            #total inlet counter for current outlet
            for k,n in enumerate(repartition):   #get iteration index and number of inlets assigned to that iteration
                for j in range(n):               #loop over each inlet in current iteration
                    inlets.append((inlets_current[i,0], inlets_current[i,1], inlets_current[i,2], inlets_current[i,3], inlets_current[i,4], k))
                    i += 1
        self.inlets  = pd.DataFrame(inlets,       columns = ['x','y', 'outlet', 'outletx', 'outlety', 'inlet_iteration']) #store as pandas df for easier handling
        self.outlets = pd.DataFrame(self.outlets, columns = ['x','y', 'outlet_iteration']) #store as df for easier indexing

        # Raster maps
        self.maps              = {}
        self.maps['outlets']   = np.full((self.grid.ynum, self.grid.xnum), np.nan) #map of null values where each cell with an outlet will have the index of that outlet
        self.maps['nodes']     = np.full((self.grid.ynum, self.grid.xnum), np.nan) #map of null values where each cell that has a node will be updated with that node index
        self.maps['cost']      = np.zeros((self.nbr_iteration, self.grid.ynum, self.grid.xnum)) #cost of travel through each cell
        self.maps['alpha']     = np.zeros((self.nbr_iteration, self.grid.ynum, self.grid.xnum)) #cost of travel along gradient through each cell
        self.maps['beta']      = np.zeros((self.nbr_iteration, self.grid.ynum, self.grid.xnum)) #cost of travel perpendicular to gradient through each cell
        self.maps['time']      = np.zeros((self.nbr_iteration, self.grid.ynum, self.grid.xnum)) #travel time to outlet from each cell
        self.maps['karst']     = np.zeros((self.nbr_iteration, self.grid.ynum, self.grid.xnum)) #presence/absence of karst conduit in each cell
        
        # Vector maps
        self.nodes     = {} #empty dic to store nodes (key: nodeID, val: [x, y, type])
        self.edges     = {} #empty dic to store edges (key: edgeID, val: [inNode, outNode])
        self.n         = 0  #start node counter at zero 
        self.e         = 0  #start edge counter at zero
        self.geodesics = [] #empty list to store raw fast-marching path output

        # Set up fast-marching:
        #Note: AGD-HFM library has different indexing, so model dimensions must be [ynum,xnum],
        # and model extent must be [ymin,ymax, xmin,xmax] (NOT x0,y0)
        self.riemannMetric = []                    #this changes at every iteration, but cannot be stored?
        self.fastMarching = agd.Eikonal.dictIn({
            'model'             : self.settings['algorithm'],      #set algorithm from settings file ('Isotropic2', 'Isotropic3', 'Riemann2', 'Riemann3')
            'order'             : 2,               #recommended setting: 2 (replace by variable)
            'exportValues'      : 1,               #export the travel time field
            'exportGeodesicFlow': 1                #export the walker paths (i.e. the conduits)
        })
        self.fastMarching.SetRect(                 #give the fast-marching algorithm the model grid 
            sides=[[self.grid.ymin, self.grid.ymax],    #leftmost edge, rightmost edge (NOT centerpoint)
                   [self.grid.xmin, self.grid.xmax]],   #bottom edge,   top edge (NOT centerpoint)
            dims=[self.grid.ynum, self.grid.xnum])      #number of cells, number of cells
        return None

    # 1
    def _repartition_points(self, points, importance):
        '''Correct for integers in importance factors list not summing correctly to total number of points'''
        repartition = []
        proportion  = []
        total       = float(sum(importance))     #total number of points as assigned (this may not be the actual total)
        for i,n in enumerate(importance):        #get iteration index and number of points being assigned to that iteration
            proportion.append(float(n)/total)    #percent of points to use this iteration
            repartition.append(round(proportion[i]*len(points))) #number of points to use this iteration (need to round bc percentage will not result in whole number)
        repartition[-1] = len(points)-sum(repartition[0:-1])     #leftover points are assignd to last iteration
        return repartition
    
    #1
    def _distribute_outlets(self, repartition):
        '''Distribute points into iteration groups based on the corrected repartition from the importance settings'''
        outlets = []
        i = 0              #point counter
        for k,n in enumerate(repartition):        #get iteration index and number of points being assigned to that iteration
            for j in range(n):                    #loop over each point in current iteration
                outlets.append((self.outlets[i][0], self.outlets[i][1], k))   #append point with its iteration number
                i += 1                                                                #increment point index by 1
        return np.asarray(outlets)

    #1
    def _distribute_inlets(self, repartition):
        '''Distribute points into iteration groups based on the corrected repartition from the importance settings'''
        inlets = []
        i = 0              #point counter
        for k,n in enumerate(repartition):        #get iteration index and number of points being assigned to that iteration
            for j in range(n):                    #loop over each point in current iteration
                inlets.append((self.inlets[i][0], self.inlets[i][1], k, self.outlets[k,0], self.outlets[k,1]))   #append point with its iteration number
                i += 1                                                                #increment point index by 1
        return np.asarray(inlets)

    # 2
    def _compute_iterations_karst_network(self):
        """
        Compute each generation of karst conduits.
        """
        if self.settings['verbosity'] > 0:
            print('-START-')

        # Define outlets map according to outlets emplacements
        self._compute_outlets_map()  #assign outlet indices

        #Plot for debugging:
        if self.settings['verbosity'] > 1:
            f,ax = plt.subplots(1,1)
            ax.imshow(self.geology.data['geology']['data'], extent=self.grid.extent, origin='lower', cmap='gray_r', alpha=0.5) #for debugging
        
        ## Set up iteration structure:
        iteration = 0                                   #initialize total iteration counter
        for outlet_iteration in range(len(self.settings['outlets_importance'])):  #loop over outlet iterations
            if self.settings['verbosity'] > 2:
                print('Total Iteration', iteration, 'Outlet iteration:', outlet_iteration)
            outlets_current = self.outlets[self.outlets.outlet_iteration==outlet_iteration]  #get the outlets assigned to the current outlet iteration
            if self.settings['verbosity'] > 1:
                ax.scatter(outlets_current.x,outlets_current.y, c='c')         #debugging
            for o,outlet in outlets_current.iterrows():                     #loop over outlets in current outlet iteration
                if self.settings['verbosity'] > 2:
                    print('\t Current outlet index:', outlet.name)             #debugging
                if self.settings['verbosity'] > 1:
                    ax.annotate(str(outlet_iteration), (outlet.x,outlet.y))
                inlets_outlet = self.inlets[self.inlets.outlet==outlet.name]          #get the inlets assigned to current outlet 
                if self.settings['verbosity'] > 1:
                    print('\t Inlets assigned to this outlet:\n', inlets_outlet)
                for inlet_iteration in range(len(self.settings['inlets_importance'])): #loop over inlet iterations
                    if self.settings['verbosity'] > 2:
                        print('\t\t Inlet iteration:', inlet_iteration)
                    inlets_current = inlets_outlet[inlets_outlet.inlet_iteration==inlet_iteration] #get the inlets assigned to the current inlet iteration
                    if self.settings['verbosity'] > 2:
                        print(inlets_current)                                                         #debugging
                    if self.settings['verbosity'] > 1:
                        ax.scatter(inlets_current.x, inlets_current.y)
                    for i,inlet in inlets_current.iterrows():                                 #loop over inlet in current inlet iteration
                        if self.settings['verbosity'] > 1:
                            ax.annotate(str(outlet_iteration)+'-'+str(inlet_iteration), (inlet.x,inlet.y))  #debugging
                        self.outlets.loc[self.outlets.index==outlet.name, 'iteration'] = iteration   #store total iteration counter
                        self.inlets.loc[ self.inlets.index ==inlet.name,  'iteration'] = iteration   #store total iteration counter
                    
                    #Compute travel time maps and conduit network:
                    if self.settings['algorithm'] == 'Isotropic2':
                        self._compute_cost_map(iteration)  
                        self._compute_time_map_isotropic(iteration) 
                        self._compute_karst_map(iteration) 
                    elif self.settings['algorithm'] == 'Riemann2':
                        self._compute_cost_map(iteration)     
                        self._compute_alpha_map(iteration) 
                        self._compute_beta_map(iteration) 
                        self._compute_riemann_metric(iteration) 
                        self._compute_time_map_riemann(iteration) 
                        self._compute_karst_map(iteration)  
                    else:
                        print('Unrecognized algorithm', self.settings['algorithm'])
                    iteration = iteration + 1                                              #increment total iteration number by 1 
                    
                    if self.settings['verbosity'] > 0:
                        print('iteration:{}/{}'.format(iteration+1,self.nbr_iteration))
        if self.settings['verbosity'] > 0:
            print('- END -')
        return None

    # 2
    def _compute_outlets_map(self):
        """
        Compute the outlets map (array indicating location of outlets as their index and everywhere else as nan).
        """
        for i,outlet in self.outlets.iterrows():
            X = int(math.ceil((outlet.x - self.grid.x0 - self.grid.dx/2) / self.grid.dx))
            Y = int(math.ceil((outlet.y - self.grid.y0 - self.grid.dy/2) / self.grid.dy))
            self.maps['outlets'][Y][X] = i 
        return None
    
    # 2
    def _compute_cost_map(self, iteration):
        """
        Compute the cost map (how difficult it is to traverse each cell).
        """
        # If it's the first iteration, iniatialize the cost map according to the geological settings.
        if iteration == 0:
            # Geology
            if self.geology.data['geology']['mode'] == 'null':
                self.maps['cost'][0] = np.full((self.grid.ynum, self.grid.xnum), self.settings['cost_aquifer']) #every cell has the same travel cost and is part of the aquifer
            elif self.geology.data['geology']['mode'] == 'image':
                self.maps['cost'][0] = np.where(self.geology.data['geology']['data']==1, self.settings['cost_aquiclude'], self.settings['cost_aquifer'])
            elif self.geology.data['geology']['mode'] == 'import' or self.geology.data['geology']['mode'] == 'csv' or self.geology.data['geology']['mode'] == 'gslib':

                tableFMM = {}
                if len(self.settings['geology_id']) != len(self.settings['geology_cost']):
                    print("- _compute_cost_map() - Error : number of lithologies does not match with number of FMM code.")
                    sys.exit()
    
                for geology, codeFMM in zip(self.settings['geology_id'], self.settings['geology_cost']):
                    if geology in self.geology.data['geology']['stats']['ID']:
                        tableFMM[geology] = codeFMM
                    else:
                        print("- initialize_costMap() - Warning : no geology n {} found.".format(geology))
                        tableFMM[geology] = self.settings['cost_out']
                    
                for y in range(self.grid.ynum):
                    for x in range(self.grid.xnum):
                        geology = self.geology.data['geology']['data'][y][x] 
                        self.maps['cost'][0][y][x] = tableFMM[geology]       
            else:
                print('geology mode', self.geology.data['geology']['mode'], 'not supported')
                sys.exit()

            # Faults
            self.maps['cost'][0] = np.where(self.geology.data['faults']['data'] > 0, self.settings['cost_faults'], self.maps['cost'][0])
    
            # Fractures
            self.maps['cost'][0] = np.where(self.geology.data['fractures']['data'] > 0, self.settings['cost_fractures'], self.maps['cost'][0])

            # If out of polygon
            if self.mask is not None:
                self.maps['cost'][0] = np.where(self.mask==1, self.settings['cost_out'], self.maps['cost'][0])

        else: #if not the first iteration
            self.maps['cost'][iteration] = self.maps['cost'][iteration-1] #set cost map to be same as previous iteration
            self.maps['cost'][iteration] = np.where(self.maps['karst'][iteration-1] > 0, self.settings['cost_conduits'], self.maps['cost'][iteration]) #where karst conduits are present from previous iteration, set cost to conduit cost, elsewhere, leave unchanged            
        return None
    
    # 2
    def _compute_alpha_map(self, iteration):
        """
        Compute the alpha map: travel cost in the same direction as the gradient.
        Cost map * topography map, so that the cost is higher at higher elevations, encouraging conduits to go downgradient.
        """
        self.maps['alpha'][iteration] = self.maps['cost'][iteration] * self.geology.data['topography']['data']
        return None
    
    # 2
    def _compute_beta_map(self, iteration):
        """
        Compute the beta map: travel cost perpendicular to the gradient. 
        If beta is higher than alpha, conduits will follow the steepest gradient.
        If beta is lower than alpha, conduits will follow contours.
        """
        self.maps['beta'][iteration] = self.maps['alpha'][iteration] / self.settings['cost_ratio']
        return None

    # 2
    def _compute_riemann_metric(self, iteration):
        """
        Compute the riemann metric: Define the Riemannian metric needed as input for the anisotropic fast marching.
        """
        self.riemannMetric = agd.Metrics.Riemann.needle([self.geology.data['orientationx']['data'], self.geology.data['orientationy']['data']], self.maps['alpha'][iteration], self.maps['beta'][iteration])
        return None
    
    # 2
    def _compute_time_map_isotropic(self, iteration):
        """
        Compute the travel time map (how long it takes to get to the outlet from each cell),
        using the isotropic agd-hfm fast-marching algorithm, and store travel time map.
        Note: the AGD-HFM library uses different indexing, so x and y indices are reversed for inlets and outlets.
        """
        self.fastMarching['seeds']  = np.rot90([self.outlets[self.outlets.iteration==iteration].x, self.outlets[self.outlets.iteration==iteration].y], k=3)         #set the outlets for this iteration
        self.fastMarching['tips']   = np.rot90([self.inlets[self.inlets.iteration==iteration].x, self.inlets[self.inlets.iteration==iteration].y], k=3) #select inlets for current iteration
        self.fastMarching['cost']   = self.maps['cost'][iteration]                       #set the isotropic travel cost through each cell
        self.fastMarching['verbosity'] = self.settings['verbosity']           #set verbosity of hfm run
        self.fastMarchingOutput     = self.fastMarching.Run()                 #run the fast marching algorithm and store the outputs
        self.maps['time'][iteration] = self.fastMarchingOutput['values']  #store travel time maps
        return None
    
    # 2
    def _compute_time_map_riemann(self, iteration):
        """
        Compute the travel time map (how long it takes to get to the outlet from each cell),
        using the anisotropic agd-hfm fast-marching algorithm, and store travel time map.
        Note: the AGD-HFM library uses different indexing, so x and y indices are reversed for inlets and outlets.
        """
        self.fastMarching['seeds']  = np.rot90([self.outlets[self.outlets.iteration==iteration].x, self.outlets[self.outlets.iteration==iteration].y], k=3)         #set the outlets for this iteration
        self.fastMarching['tips']   = np.rot90([self.inlets[self.inlets.iteration==iteration].x, self.inlets[self.inlets.iteration==iteration].y], k=3) #select inlets for current iteration
        self.fastMarching['metric'] = self.riemannMetric                      #set the travel cost through each cell
        self.fastMarching['verbosity'] = self.settings['verbosity']           #set verbosity of hfm run
        self.fastMarchingOutput     = self.fastMarching.Run()                 #run the fast marching algorithm and store the outputs
        self.maps['time'][iteration] = self.fastMarchingOutput['values']      #store travel time maps
        self.geodesics.append(self.fastMarchingOutput['geodesics'])         #store fastest travel paths
        return None

    # 2
    def _compute_karst_map(self, iteration):
        """
        Compute the karst map based on the paths from agd-hfm. 
        Array of all zeros, with ones in cells containing a karst conduit.
        """
        if iteration > 0:           # Get karst map from previous iteration (except for the very first iteration)
            self.maps['karst'][iteration] = self.maps['karst'][iteration-1] 
        
        #Debugging:
        #f1,ax1 = plt.subplots(1,1, figsize=(10,10))  #debugging
        #ax1.imshow(self.maps['karst'][iteration], origin='lower', extent=self.grid.extent, cmap='gray_r')
        #ax1.imshow(self.maps['nodes'], origin='lower', extent=self.grid.extent, cmap='gray_r')
        #ax1.scatter(self.outlets[self.outlets.iteration==iteration].x, self.outlets[self.outlets.iteration==iteration].y, c='c', s=100)
        #ax1.scatter(self.inlets[self.inlets.iteration==iteration].x, self.inlets[self.inlets.iteration==iteration].y, c='orange', s=100)   

        #Loop over conduit paths generated by fast marching:
        for path in self.fastMarchingOutput['geodesics']:   #loop over conduit paths in this iteration (there is one path from each inlet)
            merge = False                                   #reset indicator for whether this conduit has merged with an existing conduit
            for p in range(path.shape[1]):                  #loop over points making up this conduit path
                point = path[:,p]                           #get coordinates of current point
                [[ix,iy],error]  = self.fastMarching.IndexFromPoint(point) #convert to coordinates to indices
                #ax1.scatter(point[1],point[0], c='g',s=5)  #debugging

                #Place nodes and links:
                if np.isnan(self.maps['nodes'][ix,iy]):                                    #if there is no existing conduit node here 
                    if ~np.isnan(self.maps['outlets'][ix,iy]):                              #if there is an outlet here (cell value is not nan)
                        outlet = self.outlets.iloc[int(self.maps['outlets'][ix,iy])]         #get the outlet coordinates using the ID in the outlets map
                        self.nodes[self.n]             = [outlet.y, outlet.x, 'outfall']     #add a node at the outlet coordinates (with the node type for SWMM)
                        self.maps['nodes'][ix,iy] = self.n                                   #update node map with node index
                        #ax1.scatter(outlet.x,outlet.y, marker='o', c='b')                   #debugging
                        if p > 0:                                                           #if this is not the first point (i.e. the inlet) in the current path
                            if merge == False:                                               #if this conduit has not merged with an existing conduit
                                self.edges[self.e] = [self.n-1, self.n]                       #add an edge connecting the previous node to the current node
                                self.e = self.e+1                                             #increment edge counter up by one
                                #ax1.plot((self.nodes[self.n][0], self.nodes[self.n-1][0]),(self.nodes[self.n][1], self.nodes[self.n-1][1]))
                            else:                                                          #if this conduit HAS merged with an existing conduit
                                [[fromix,fromiy],error]  = self.fastMarching.IndexFromPoint(path[:,p-1]) #get xy indices of previous point in current conduit path
                                n_from = self.maps['nodes'][fromix,fromiy]                    #get node index of the node already in the cell where the previous point was
                                self.edges[self.e] = [n_from, self.n]                         #add an edge connecting existing conduit node to current node
                                self.e = self.e+1                                             #increment edge counter up by one
                                #ax1.plot((self.nodes[self.n].x, self.nodes[n_from].x),(self.nodes[self.n].y, self.nodes[n_from].y))
                        self.n = self.n+1                                                   #increment node counter up by one
                    else:                                                                  #if there is NOT an outlet here
                        if p > 0:                                                           #if this is not the first point in the current path
                            self.nodes[self.n] = [point[0], point[1], 'junction']            #add a junction node here (with the node type for SWMM)
                            self.maps['nodes'][ix,iy] = self.n                               #update node map with node index
                            #ax1.scatter(point[1],point[0], marker='o', edgecolor='g', facecolor='none')   #debugging
                            if merge == False:                                              #if this conduit has not merged with an existing conduit
                                self.edges[self.e] = [self.n-1, self.n]                      #add and edge connecting the previous node to the current node
                                self.e = self.e+1                                            #increment edge counter up by one
                                #ax1.plot((self.nodes[self.n][1], self.nodes[self.n-1][1]),(self.nodes[self.n][0], self.nodes[self.n-1][0]), c='gold', marker=None)
                            else:                                                           #if this conduit HAS merged with an existing conduit
                                [[fromix,fromiy],error]  = self.fastMarching.IndexFromPoint(path[:,p-1]) #get xy indices of previous point in current conduit path
                                n_from = self.maps['nodes'][fromix,fromiy]                   #get node index of the node already in the cell where the previous point was
                                self.edges[self.e] = [n_from, self.n]                        #add an edge connecting existing conduit node to current node
                                self.e = self.e+1                                            #increment edge counter up by one   
                                merge = False                                                #reset merge indicator to show that current conduit has left the path of the existing conduit
                                #ax1.plot((self.nodes[self.n][1], self.nodes[n_from][1]),(self.nodes[self.n][0], self.nodes[n_from][0]), c='m', marker=None)
                        else:                                                               #if this is the first point in current path
                            self.nodes[self.n] = [point[0], point[1], 'inlet']               #add an inlet node here (with the node type for SWMM)
                            self.maps['nodes'][ix,iy] = self.n                               #update node map with node index
                            #ax1.scatter(point[1],point[0], marker='o', edgecolor='sienna', facecolor='none')
                        self.n = self.n+1                                                   #increment node counter up by one
                elif ~np.isnan(self.maps['nodes'][ix,iy]):                                 #if there is already a node in this cell (either because there is a conduit here, or because there are two nodes in the same cell)
                    n_existing = self.maps['nodes'][ix,iy]                                  #get index of node already present in current cell
                    if merge == True:                                                       #if this conduit has already merged into an existing conduit
                        pass                                                                 #skip this node (there is already a node here)
                    elif n_existing == self.n-1:                                            #if existing index is only one less than next node to be added index, this is a duplicate node and can be skipped
                        pass                                                                 #skip this node (duplicate)
                    else:                                                                   #if existing node index is >1 less than next node to be added index
                        if p > 0:                                                           #if this is not the first point in the current path
                            self.edges[self.e] = [self.n-1, n_existing]                      #add an edge connecting most recently added node and existing node in cell
                            self.e = self.e+1                                                #increment edge counter up by one
                            merge = True                                                     #add a flag indicating that this conduit has merged into an existing one
                            #ax1.plot((self.nodes[self.n-1][1], self.nodes[n_existing][1]),(self.nodes[self.n-1][0], self.nodes[n_existing][0]), c='r', marker=None)
                        else:                                                                #if this is the first point in the current path (i.e. the inlet is on an exising conduit)
                            self.nodes[self.n] = [point[0], point[1], 'inlet']                #add a node here (with the node type for SWMM)- this will cause there to be two nodes in the same cell
                            self.maps['nodes'][ix,iy] = self.n                                #update node map with node index
                            self.n = self.n+1                                                 #increment node counter by 1
                            #ax1.scatter(point[1],point[0], marker='o', edgecolor='g', facecolor='none')  #debugging
                
                self.maps['karst'][iteration][ix,iy] = 1                               #update karst map to put a conduit in current cell
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

        if X == 0:
            for row in borders_conditions[0]:
                time_values.append(self.maps['time'][iteration][Y+row[1]][X+row[0]])
            rank  = time_values.index(min(time_values))
            moove = borders_conditions[0][rank]
        elif Y == 0:
            for row in borders_conditions[1]:
                time_values.append(self.maps['time'][iteration][Y+row[1]][X+row[0]])
            rank  = time_values.index(min(time_values))
            moove = borders_conditions[1][rank]
        elif X == self.grid.xnum:
            for row in borders_conditions[2]:
                time_values.append(self.maps['time'][iteration][Y+row[1]][X+row[0]])
            rank  = time_values.index(min(time_values))
            moove = borders_conditions[2][rank]
        elif Y == self.grid.ynum:
            for row in borders_conditions[3]:
                time_values.append(self.maps['time'][iteration][Y+row[1]][X+row[0]])
            rank  = time_values.index(min(time_values))
            moove = borders_conditions[3][rank]
        elif (X == 0) and (Y == 0):
            for row in corners_conditions[0]:
                time_values.append(self.maps['time'][iteration][Y+row[1]][X+row[0]])
            rank  = time_values.index(min(time_values))
            moove = corners_conditions[0][rank]
        elif (X == 0) and (Y == self.grid.ynum):
            for row in corners_conditions[1]:
                time_values.append(self.maps['time'][iteration][Y+row[1]][X+row[0]])
            rank  = time_values.index(min(time_values))
            moove = corners_conditions[1][rank]
        elif (X == self.grid.xnum) and (Y == self.grid.ynum):
            for row in corners_conditions[2]:
                time_values.append(self.maps['time'][iteration][Y+row[1]][X+row[0]])
            rank  = time_values.index(min(time_values))
            moove = corners_conditions[2][rank]
        elif (X == self.grid.xnum) and (Y == 0):
            for row in corners_conditions[3]:
                time_values.append(self.maps['time'][iteration][Y+row[1]][X+row[0]])
            rank  = time_values.index(min(time_values))
            moove = corners_conditions[3][rank]
        else:
            print('Error: none of the conditions are met.')

        dx = moove[0] * self.step
        dy = moove[1] * self.step
        
        return (dx,dy)


    ###########################
    # Visualization functions #
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
        elif data == 'topography':
            d = self.geology.data['topography']['data']

        if mask==True:
            if   data == 'geology':
                d = self.geology_masked['geology']
            elif data == 'faults':
                d = self.geology_masked['faults']
            elif data == 'fractures':
                d = self.geology_masked['fractures']
            elif data == 'topography':
                d = self.geology_masked['topography']

        im1 = ax1.imshow(d, extent=self.grid.extent, origin='lower', cmap=cmap) 
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
                
    def _show_maps(self, sim=-1, iteration=-1, cmap='binary'):
        """
        Show the simulated karst network as an image.
        """
        karst_network = self.karst_simulations[sim]
        
        fig, ([ax1,ax2],[ax3,ax4]) = plt.subplots(2, 2, sharex=True, sharey=True)
        fig.suptitle('Karst Network', fontsize=16)
        
        ax1.imshow(karst_network.maps['outlets'], extent=self.grid.extent, origin='lower', cmap=cmap)
        ax1.set_title('Outlets')
        
        ax2.imshow(karst_network.maps['cost'][iteration], extent=self.grid.extent, origin='lower', cmap=cmap)
        ax2.set_title('Cost')
        
        ax3.imshow(karst_network.maps['time'][iteration], extent=self.grid.extent, origin='lower', cmap=cmap)
        ax3.set_title('Time')
        
        ax4.imshow(karst_network.maps['karst'][iteration], extent=self.grid.extent, origin='lower', cmap=cmap)
        ax4.set_title('Karst')

        fig.subplots_adjust(hspace=0.5)
        plt.show()
        return None
        

    def show(self, data=None, title=None, probability=False):
        """
        Show the entire study domain (defaults to showing most recent simulation).
        """        
        if data is None:
            data = self.karst_simulations[-1]
        
        if probability == True:
            data = self._compute_average_paths()
            
        fig = plt.figure(figsize=(20,10))
        
        fig.add_subplot(131, aspect='equal')
        plt.xlabel('Cost array'+str(data.maps['cost'][-1].shape))
        plt.imshow(data.maps['cost'][-1], extent=self.grid.extent, origin='lower', cmap='gray') #darker=faster
        plt.colorbar(shrink=0.35)

        fig.add_subplot(132, aspect='equal')
        plt.xlabel('Travel time array'+str(data.maps['time'][-1].shape))
        plt.imshow(data.maps['time'][-1], extent=self.grid.extent, origin='lower', cmap='cividis') #darker=faster
        plt.colorbar(shrink=0.35)

        fig.add_subplot(133, aspect='equal')
        plt.xlabel('Karst array'+str(data.maps['time'][-1].shape))
        plt.imshow(data.maps['karst'][-1], extent=self.grid.extent, origin='lower', cmap='gray_r') #darker=conduits
        plt.colorbar(shrink=0.35)
        i = plt.scatter(data.points['inlets'].x,  data.points['inlets'].y,  c='orange')
        o = plt.scatter(data.points['outlets'].x, data.points['outlets'].y, c='steelblue')
        p = matplotlib.patches.Rectangle((0,0),0,0, ec='r', fc='none')
        if self.settings['data_has_polygon']:
            closed_polygon = self.polygon.polygon[:]
            closed_polygon.append(closed_polygon[0])
            x,y = zip(*closed_polygon)
            plt.plot(x,y, color='red', label='polygon')
        plt.legend([i,o,p], ['inlets', 'outlets', 'catchment'], loc='upper right')
        
        if title is not None:
            fig.suptitle(title, fontsize=16)
        plt.show()
        return None

    def show_network(self, data=None, simplify=False, ax=None, labels=['inlets', 'outlets'], title=None, cmap='cividis'):
        """
        Show the karst network as a graph with nodes and edges. Defaults to showing latest iteration.
        Inputs: 
        data: karst simulation object containing nodes, edges, points, etc. Can be obtained from self.karst_simulations[i]
        ax: axis to plot on
        label: None or list of strings ['nodes','edges','inlets','outlets'], indicating which components to label
        title: string, title of plot
        cmap: string, colormap to use when plotting
        """

        if ax == None:
            fig,ax = plt.subplots(figsize=(10,10))
            ax.set_aspect('equal')

        if data == None:
            data = self.karst_simulations[-1]

        if simplify == True:
            nodes = data.network['nodes']   #get all nodes
            nodes_simple = data.network['karstnet'].graph_simpl.nodes  #get indices of only the nodes in the simplified graph
            nodes_simple = {key: nodes[key] for key in nodes_simple}   #make df of only the nodes in the simplified graph, for plotting
            edges = data.network['edges']   #get all edges
            edges_simple = data.network['karstnet'].graph_simpl.edges  #get indices of only the nodes in the simplified graph
            edges_simple = {key: edges[key] for key in edges_simple}   #make df of only the nodes in the simplified graph, for p
        else:
            nodes = pd.DataFrame.from_dict(data.network['nodes'], orient='index', columns=['x','y','type']) #convert to pandas for easier plotting
            edges = pd.DataFrame.from_dict(data.network['edges'], orient='index', columns=['inNode','outNode']) 

        #Set up data for plotting:
        fromX = nodes.x.loc[edges.inNode]      #calculate coordinates for link start and end points
        fromY = nodes.y.loc[edges.inNode]
        toX   = nodes.x.loc[edges.outNode]
        toY   = nodes.y.loc[edges.outNode]

        #Plot nodes and edges:
        n = ax.scatter(nodes.y,                     nodes.x,                     c='k',         s=5)  #scatterplot nodes
        i = ax.scatter(data.points['inlets'].x,  data.points['inlets'].y,  c='orange',    s=30) #scatterplot inlets
        o = ax.scatter(data.points['outlets'].x, data.points['outlets'].y, c='steelblue', s=30) #scatterplot outlets
        e = matplotlib.lines.Line2D([0],[0])                                                  #line artist for legend 
        for ind in edges.index:                                                               #loop over edge indices
            ax.plot((fromY.iloc[ind], toY.iloc[ind]), (fromX.iloc[ind], toX.iloc[ind]), c=plt.cm.get_cmap(cmap)(ind/len(edges)))  #plot each edge, moving along color gradient to show order
        
        #Add labels:
        if 'nodes' in labels:                                         #label node indices
            for ind in nodes.index:                                   #loop over node indices
                ax.annotate(str(ind), xy=(nodes.y[ind]-10, nodes.x[ind]))  #annotate slightly to left of each node
        if 'edges' in labels:                                         
            for ind in edges.index:                                   
                ax.annotate(str(ind), xy=(edges.y[ind]-10, edges.x[ind]))  #annotate slightly to left of each edge
        if 'inlets' in labels:                                        
            for index,inlet in data.points['inlets'].iterrows(): 
                ax.annotate(str(int(inlet.outlet))+'-'+str(int(inlet.inlet_iteration)),  xy=(inlet.x-10,  inlet.y)) 
        if 'outlets' in labels:                                       
            for index,outlet in data.points['outlets'].iterrows():                     
                ax.annotate(str(int(outlet.name)), xy=(outlet.x-10, outlet.y)) 

        #Add legend & title:
        ax.legend([i,o,n,e],['inlets','outlets','nodes','edges'])
        if title is not None:
            ax.set_title(title, fontsize=16)

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

#####################   
class KarstNetwork():
    """
    A class for storing a calculated karst network.
    """
    def __init__(self, maps, points, network, stats, settings):
        self.maps     = maps
        self.points   = points
        self.network  = network
        self.stats    = stats
        self.settings = settings