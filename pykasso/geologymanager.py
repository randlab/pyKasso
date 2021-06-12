from .fracture  import Fracture
from tqdm       import tqdm

from skimage.io        import imread
from skimage.transform import resize

import sys
import math
import mpmath
import numpy    as np
import matplotlib.pyplot as plt

class GeologyManager():
    """
    Class modeling the geologic data of the studied domain.
    """

    def __init__(self, grid):
        """
        Creates a geology manager on a Grid instance.
        This class is designed to handle geology, faults and fractures data.
        Data is stored in a the 'data' attribute as a dictionnary.
        The data model is set as follow :
        {data_key :
            {
            data  : array,
            img   : array,
            stats : dict,
            mode  : str
            }
        }

        Parameters
        ----------
        grid : Grid instance
            The geology manager must be set on the studied grid.

        Examples
        --------
        >>> grid = pk.Grid(0, 0, 0, 10, 10, 10, 1, 1, 1)
        >>> geol = pk.GeologyManager(grid)
        """
        self.data    = {}
        self.grid    = grid

    def set_data_null(self, data_key):
        """
        Sets data to 'null' mode for the indicated data key in the 'data' dictionary attribute.

        For data keys 'topography', 'orientationx' and 'orientationy' sets : np.zeros((nx,ny))

        Otherwise, sets : np.zeros((nx,ny,nz))

        Parameters
        ----------
        data_key : str
            Type of data : 'geology', 'faults', 'fractures', 'topography', 'orientationx' or 'orientationy'.

        Examples
        --------
        >>> geol.set_data_null('geology')
        >>> print(geol.data['geology'])
        """
        self.data[data_key] = {}

        if data_key in ['topography', 'orientationx', 'orientationy']:
            self.data[data_key]['data'] = np.zeros((self.grid.nx, self.grid.ny), dtype=np.int_)
            self.data[data_key]['img']  = np.zeros((self.grid.nx, self.grid.ny), dtype=np.int_)
        else:
            self.data[data_key]['data'] = np.zeros((self.grid.nx, self.grid.ny, self.grid.nz), dtype=np.int_)
            self.data[data_key]['img']  = np.zeros((self.grid.nx, self.grid.ny, self.grid.nz), dtype=np.int_)

        self.data[data_key]['mode'] = 'null'
        return None

    def set_data_from_csv(self, data_key, datafile_location):
        """
        Sets data from a csv file for the indicated data key in the 'data' dictionary attribute.
        Delimiter is ','.

        If Grid.nz > 1, the layer is horizontally repeated.

        x-direction : from West to East
        y-direction : from North to South

        Parameters
        ----------
        data_key : str
            Type of data : 'geology', 'topography', 'surface' (orientation), 'faults' or 'fractures'.
        datafile_location : str
            Path of the datafile.

        Examples
        --------
        >>> geol.set_data_from_csv('geology', 'geology.csv')
        >>> print(geol.data['geology'])
        """
        self.data[data_key] = {}
        try:
            data = np.genfromtxt(datafile_location, delimiter=',', dtype=np.float_)
        except:
            print('- set_data_from_csv() - Error : unable to read datafile.')
            raise
        if data_key in ['topography', 'surface']:
            self.data[data_key]['data'] = np.transpose(data)
            self.data[data_key]['img']  = np.flipud(data)
        else:
            self.data[data_key]['data'] = np.zeros((self.grid.nx, self.grid.ny, self.grid.nz), dtype=np.float_)
            self.data[data_key]['img']  = np.zeros((self.grid.nx, self.grid.ny, self.grid.nz), dtype=np.float_)
            for z in range(self.grid.nz):
                self.data[data_key]['data'][:,:,z] = np.transpose(data)
                self.data[data_key]['img'] [:,:,z] = np.flipud(data)
        self.data[data_key]['mode'] = 'csv'
        return None

    def set_data_from_gslib(self, data_key, datafile_location):
        """
        Sets data from a gslib file for the indicated data key in the 'data' dictionary attribute.

        x-direction : from west to east
        y-direction : from south to north
        z-direction : from bottom to top

        Parameters
        ----------
        data_key : str
            Type of data : 'geology', 'faults' or 'fractures'.
        datafile_location : str
            Path to the datafile.

        Examples
        --------
        >>> geol.set_data_from_gslib('geology', 'geology.gslib')
        >>> print(geol.data['geology'])
        """
        self.data[data_key] = {}
        try:
            data = np.genfromtxt(datafile_location, skip_header=3, dtype=np.float_)    #read in gslib file as float numpy array without header rows
        except:
            print('- set_data_from_gslib() - Error : unable to read datafile.')
            raise
        #if data_key is "geology":
        #a[a==0] = np.nan                                                             #replace zeros with nans (array must be float first)
        data = np.reshape(data, (self.grid.nx, self.grid.ny, self.grid.nz), order='F')#reshape to xy grid using Fortran ordering
        self.data[data_key]['data'] = data                                            #store
        self.data[data_key]['img']  = np.flipud(np.transpose(data, (1,0,2)))
        self.data[data_key]['mode'] = 'gslib'
        return None

    def set_data_from_image(self, data_key, datafile_location):
        """
        Sets data from an image for the indicated data key in the 'data' dictionary attribute.
        The size of the image should be the same that the size of the grid, but this is optional.
        If nz > 1, the layer is horizontally repeated.
        This method usage is not recommended, it should be used only for quick testing.

        Parameters
        ----------
        data_key : str
            Type of data : 'geology', 'faults' or 'fractures'.
        datafile_location : str
            Path of the datafile.

        Examples
        --------
        >>> geol.set_data_from_image('geology', 'geology.jpg')
        >>> print(geol.data['geology'])
        """
        try:
            image = imread(datafile_location)
        except:
            print('- set_data_from_image() - Error : unable to read datafile.')
            raise
        self.data[data_key] = {}
        # Store image for show() methods
        image = (image[:,:,0] == 0)*1
        image = np.rint(resize(image, (self.grid.ny, self.grid.nx), anti_aliasing=False, preserve_range=True))
        self.data[data_key]['img'] = image
        # Store data
        self.data[data_key]['data'] = np.zeros((self.grid.nx, self.grid.ny, self.grid.nz), dtype=np.float_)
        image_data = np.flipud(image)                #we need to flip it since the data reading start from bottom
        image_data = np.transpose(image_data, (1,0)) #imread return image with mxn format so we also need to transpose it
        for z in range(self.grid.nz):
            self.data[data_key]['data'][:,:,z] = image_data
        self.data[data_key]['mode'] = 'image'
        return None

    def _set_data_from_pickle(self, data_key, datafile_location):
        """
        Sets data from a pickle data file.

        Parameters
        ----------
        data_key : str
            Type of data : 'geology', 'faults' or 'fractures'.
        datafile_location : str
            Path of the datafile.

        Examples
        --------
        >>> geol.set_data_from_pickle('geology', 'geology.pickle')
        >>> print(geol.data['geology'])
        """
        return None

    def generate_orientations(self, surface):
        """
        Generates maps of x and y components of orientation.

        Parameters
        ------------
        surface : 2D numpy-array
            A 2D array of the elevation or potential in cell, used to calculate the orientation (i.e. slope or dip) of that cell.
            Use either the land surface ('topography'), the surface of the bottom of the karst unit ('contact'), or the potential returned by the geologic model.

        Examples
        --------
        >>> geol.generate_orientations(surf)
        """
        self.data['orientationx'] = {}
        self.data['orientationx']['data'] = np.zeros((self.grid.nx, self.grid.ny))
        self.data['orientationy'] = {}
        self.data['orientationy']['data'] = np.zeros((self.grid.nx, self.grid.ny))

        self.data['orientationx']['data'], self.data['orientationy']['data'] = np.gradient(surface, self.grid.dx, self.grid.dy, axis=(0,1))   #x and y components of gradient in each cell of array

        self.data['orientationx']['img']  = np.flipud(np.transpose(self.data['orientationx']['data']))
        self.data['orientationy']['img']  = np.flipud(np.transpose(self.data['orientationy']['data']))
        self.data['orientationx']['mode'] = 'topo'
        self.data['orientationy']['mode'] = 'topo'
        return None

    def generate_fractures(self, fractures_densities, fractures_alpha, fractures_min_orientation, fractures_max_orientation, fractures_min_dip, fractures_max_dip, fractures_min_length, fractures_max_length):
        """
        Generates fractures as Fracture instances according to the parameters.

        Parameters
        ----------
        fractures_densities : array
            Fracture densities for each fracture family.
        fractures_alpha : array
            Degree of power law for each fracture family.
        fractures_min_orientation : array
            Fractures minimum orientation for each fracture family.
        fractures_max_orientation : array
            Fractures maximum orientation for each fracture family.
        fractures_min_dip : array
            Fractures minimum dip for each fracture family.
        fractures_max_dip : array
            Fractures maximum dip for each fracture family.
        fractures_min_length : array
            The minimum length of the fractures for each fracture family.
        fractures_max_length : array
            The maximum length of the fractures for each fracture family.

        Examples
        --------
        >>> geol.generate_fractures(fractures_densities = [0.00005,0.0001],
                                    fractures_alpha = [2, 2],
                                    fractures_min_orientation = [340, 70],
                                    fractures_max_orientation = [20, 110],
                                    fractures_min_dip = [80, 40],
                                    fractures_max_dip = [90, 50],
                                    fractures_min_length : [ 100, 100],
                                    fractures_max_length : [4000, 4000])
        >>> print(geol.fractures)
        """

        self.data['fractures'] = {}
        self.data['fractures']['data']  = np.zeros((self.grid.nx, self.grid.ny, self.grid.nz))
        self.fractures = []
        fracture_id = 0

        ## Redefine fracturation domain
        lenmax = max(fractures_max_length)
        xd0, xd1, yd0, yd1, zd0, zd1 = self.grid.xmin, self.grid.xmax, self.grid.ymin, self.grid.ymax, self.grid.zmin, self.grid.zmax
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
        self.fractures_numbers = np.array(fractures_densities) * area

        ## Calculate fractures in each family
        for frac_family in range(len(fractures_densities)):

            # Define all the constants required to randomly draw the length in distribution
            palpha = (1-fractures_alpha[frac_family])
            invpalpha = 1/palpha
            fmina = fractures_min_length[frac_family]**palpha
            frangea = fractures_max_length[frac_family]**palpha - fmina

            # Define all the constants required for the orientation distribution
            orientation_min = fractures_min_orientation[frac_family]
            orientation_max = fractures_max_orientation[frac_family]
            if orientation_min > orientation_max:
                orientation_min = orientation_min - 360

            orientation_mean_angle = (orientation_max + orientation_min)  / 2
            orientation_std_angle  = (orientation_max - orientation_mean_angle) / 3
            orientation_mean_angle = math.radians(orientation_mean_angle)
            orientation_std_angle  = math.radians(orientation_std_angle)

            # Define all the constants required for the dip distribution
            if (self.grid.nz > 1):
                dip_min = fractures_min_dip[frac_family]
                dip_max = fractures_max_dip[frac_family]
                if dip_min > dip_max:
                    dip_min = dip_min - 360

                dip_mean_angle = (dip_max + dip_min)  / 2
                dip_std_angle  = (dip_max - dip_mean_angle) / 3
                dip_mean_angle = math.radians(dip_mean_angle)
                dip_std_angle  = math.radians(dip_std_angle)

            # Computes the kappa value for the Von Mises orientation distribution
            orientation_func  = lambda k : orientation_std_angle**2 - 1 + mpmath.besseli(1,k) / mpmath.besseli(0,k)
            orientation_kappa = mpmath.findroot(orientation_func,1)

            if (self.grid.nz > 1):
                dip_func  = lambda k : dip_std_angle**2 - 1 + mpmath.besseli(1,k) / mpmath.besseli(0,k)
                dip_kappa = mpmath.findroot(dip_func,1)

            # Generate poisson number for each fracture family
            real_frac_number = np.random.poisson(self.fractures_numbers[frac_family])

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

                # FRACTURE ORIENTATION/DIP from von Mises distribution
                frac_orientation = np.random.vonmises(orientation_mean_angle, orientation_kappa)
                if (self.grid.nz > 1):
                    frac_dip = np.random.vonmises(dip_mean_angle, dip_kappa)
                else:
                    frac_dip = math.radians(0)

                # Calculate normal Vector
                x = np.sin(frac_orientation + np.radians(90)) * np.cos(frac_dip)
                y = np.cos(frac_orientation + np.radians(90)) * np.cos(frac_dip)
                z = np.sin(frac_dip)
                a = y
                b = -x
                c = -z

                # Store calculated fracture
                self.fractures.append(Fracture(fracture_id, frac_family, [xm,ym,zm], radius, frac_orientation, frac_dip, [a,b,c]))
                fracture_id += 1
        return None

    def _float_eq(self, a, b, tolerance=1e-5):
        """
        Returns True if the difference between a and b
        is lower than tolerance.
        """
        if np.all( np.isnan(a) ) :
            if np.all( np.isnan(b) ) :
                return True
            else :
                return False
        return np.all( abs(a-b) < tolerance )

    def _unit_intersect(self, n, d):
        """
        Computes the intersection between a unit circle
        and a line on a 2D plane.

        Parameters
        ----------
        The unit circle is centered on the origin (x=0, y=0) and has
        a radius of 1.

        The line is defined by two parameters :

        n : numpy array of size 2
            The normal vector perpendicular to the line
            Warning : The norm of n must be 1.

        d : floating point value
            Defines the position of the line, it is a
            signed distance to the origin x=0, y=0.
            It can be positive or negative.

        Returns
        -------
        [X1, X2, Y1, Y2] : numpy array with 4 floating point values
            the coordinates of the two points of intersection when the
            exists. Returns 4 np.NaN if there is no intersection

        """
        nx, ny = n[0], n[1] # For readibility

        if not self._float_eq( np.linalg.norm(n), 1) :
            print("WARNING - unitcintl function : norm of n must be equal to 1")

        if np.abs(d) > 1 : # Case with no intersection
            return np.array( [ np.NaN,np.NaN,np.NaN,np.NaN ] )

        sd = np.sqrt(1 - d**2)

        if ny !=0 :
            X1 = d * nx + ny * sd
            X2 = d * nx - ny * sd
            Y1 = d * ny - nx * sd
            Y2 = d * ny + nx * sd
        else :
            X1 = X2 = d
            Y1 = sd
            Y2 = -sd

        return np.array( [X1,X2,Y1,Y2] )

    def _disk_zplane_intersect(self, center, n, R, zi):
        """
        Computes the intersection between a disk
        and a horizontal plane (constant z plane) in 3D.

        Parameters
        ----------
        The disk is defined by three parameters :

        center : numpy array of size 3
            the 3D coordinates of the center of the disk.

        n : numpy array of size 3
            the normal vector perpendicular to the disk
            Warning : The norm of n must be 1.

        R : floating point value
            Radius of the disk

        zi : floating point value
            Position of the horizontal plane along the z axis

        Returns
        -------
        [x1, y1, x2, y2] : numpy array with 4 floating point values
            the coordinates of the two extremities of the intersection
            between the plane and disk if there is an intersection.
            Returns 4 np.NaN if there is no intersection.
        """

        xc, yc, zc = center[0], center[1], center[2] # For readibility

        tau = R**2 - (zi - zc)**2

        if tau <= 0 : # Avoid computing square root of negative number
            #print("The plane does not touch the sphere")
            x1, y1, x2, y2 = np.NaN, np.NaN, np.NaN, np.NaN

        else :
            tau = np.sqrt( tau )
            b = n[2] * (zc - zi) / tau / np.sqrt(1 - n[2]**2)

            if b > 1 :
                #print("The plane does not touch the disk")
                x1, y1, x2, y2 = np.NaN, np.NaN, np.NaN, np.NaN

            else :
                n2 = n.copy()[0:2] # Projection on horizontal plane
                n2 /= np.linalg.norm(n2)
                xint = self._unit_intersect(n2, b)

                x1 = tau * xint[0] + xc
                x2 = tau * xint[1] + xc
                y1 = tau * xint[2] + yc
                y2 = tau * xint[3] + yc

        return np.array([x1, y1, x2, y2])

    def _disk_xplane_intersect(self, center, n, R, xi):
        """
        Computes the intersection between a disk
        and a vertical plane (constant x plane) in 3D.

        Parameters
        ----------
        The disk is defined by three parameters :

        center : numpy array of size 3
            the 3D coordinates of the center of the disk.

        n : numpy array of size 3
            the normal vector perpendicular to the disk
            Warning : The norm of n must be 1.

        R : floating point value
            Radius of the disk

        xi : float
            Position of the vertical plane along the x axis

        Returns
        -------
        [y1, z1, y2, z2] : numpy array with 4 floating point values
            the coordinates of the two extremities of the intersection
            between the plane and disk if there is an intersection.
            Returns 4 np.NaN if there is no intersection.
        """

        xc, yc, zc = center[0], center[1], center[2] # For readibility

        tau = R**2 - (xi - xc)**2

        if tau <= 0 : # Avoid computing square root of negative number
            #print("The plane does not touch the sphere")
            y1, z1, y2, z2 = np.NaN, np.NaN, np.NaN, np.NaN

        else :
            tau = np.sqrt( tau )
            b = n[0] * (xc - xi) / tau / np.sqrt(1 - n[0]**2)

            if b > 1 :
                #print("The plane does not touch the disk")
                y1, z1, y2, z2 = np.NaN, np.NaN, np.NaN, np.NaN

            else :
                n2 = np.array( [n[1], n[2]] ) # Projection on vertical x plane
                n2 /= np.linalg.norm(n2)
                xint = self._unit_intersect(n2, b)

                y1 = tau * xint[0] + yc
                y2 = tau * xint[1] + yc
                z1 = tau * xint[2] + zc
                z2 = tau * xint[3] + zc

        return np.array([y1, z1, y2, z2])

    def _disk_yplane_intersect(self, center, n, R, yi):
        """
        Computes the intersection between a disk
        and a vertical plane (constant y plane) in 3D.

        Parameters
        ----------
        The disk is defined by three parameters :

        center : numpy array of size 3
            the 3D coordinates of the center of the disk.

        n : numpy array of size 3
            the normal vector perpendicular to the disk
            Warning : The norm of n must be 1.

        R : floating point value
            Radius of the disk

        yi : float
            Position of the vertical plane along the y axis

        Returns
        -------
        [x1, z1, x2, z2] : numpy array with 4 floating point values
            the coordinates of the two extremities of the intersection
            between the plane and disk if there is an intersection.
            Returns 4 np.NaN if there is no intersection.
        """

        xc, yc, zc = center[0], center[1], center[2] # For readibility

        tau = R**2 - (yi - yc)**2
        if tau <= 0 : # Avoid computing square root of negative number
            #print("The plane does not touch the sphere")
            x1, z1, x2, z2 = np.NaN, np.NaN, np.NaN, np.NaN

        else :
            tau = np.sqrt( tau )
            b = n[0] * (yc - yi) / tau / np.sqrt(1 - n[0]**2)

            if b > 1 :
                #print("The plane does not touch the disk")
                x1, z1, x2, z2 = np.NaN, np.NaN, np.NaN, np.NaN


            else :
                n2 = np.array( [n[0], n[2]] ) # Projection on vertical x plane
                n2 /= np.linalg.norm(n2)
                xint = self._unit_intersect(n2, b)

                x1 = tau * xint[0] + xc
                x2 = tau * xint[1] + xc
                z1 = tau * xint[2] + zc
                z2 = tau * xint[3] + zc

        return np.array([x1, z1, x2, z2])

    def _rst2d(self, m, xs, xe, ys, ye):
        """
        Rasterize a line on a 2D plane.

        Parameters
        ----------
        m : 2D numpy array

        xs, xe, ys, ye : integer values
            indices of positions on the grid m
            starting and ending location of the line

        Returns
        -------
        m : 2D numpy array
            with value 1 where the grid is touched by the line
        """

        nx, ny = m.shape

        if xe < xs : # Ensuires dx always positive
            xe, xs = xs, xe
            ye, ys = ys, ye

        dx = xe - xs
        dy = ye - ys

        if xs > nx : # The line is entirely out of the grid
            return

        if (dx + np.abs(dy)) == 0 : # Case with a single pixel
            if (ys >= 0) and (ys < ny) :
                m[xs, ys] = 1
            return

        if dx >= abs(dy): # line with greater horizontal than vertical extension
            i = np.arange(max(0, xs), min(xe, nx) )
            j = np.int32( np.round((i-xs) * dy / dx + ys ) )

            indx_ok = (j>=0) & (j<ny) # To crop the part that is out of grid
            i = i[indx_ok]
            j = j[indx_ok]

        else:
            if ye < ys : # to ensure that arange works properly
                xe, xs = xs, xe
                ye, ys = ys, ye

            if ys>ny :
                return

            j = np.arange(max(ys, 0) , min(ye,ny) )
            i = np.int32( np.round((j-ys) * dx / dy + xs ) )

            indx_ok = (i>=0) & (i<nx) # To crop the part that is out of grid
            i = i[indx_ok]
            j = j[indx_ok]

        m[i,j] = 1
        return

    def rasterize_fracture_network(self):
        """
        Rasterizes a set of fractures on a 3D grid.

        Returns
        -------
        raster_fractures : 3D numpy array
            With value 1 where the grid is touched by a fracture
        """

        dx, dy, dz = self.grid.dx, self.grid.dy, self.grid.dz
        nx, ny, nz = self.grid.nx, self.grid.ny, self.grid.nz
        x0, y0, z0 = self.grid.x0, self.grid.y0, self.grid.z0

        # Creates empty array to store raster of fracture indicators
        raster_fractures = np.zeros( (nx, ny, nz) )

        # Loop over the fractures
        for f in tqdm(self.fractures, desc="Rasterizing"):

            # Get fracture geometry
            xc, yc, zc    = f.get_position()
            nfx, nfy, nfz = f.get_normal()
            n = [nfx, nfy, nfz]
            R = f.get_radius()

            # Test orientation of the fractures
            if nfz**2 <= (nfx**2 + nfy**2) : # subvertical case
                # Vertical index of the center of the fracture
                kc = ( (zc - z0) / dz  ).astype(int)

                # Projected vertical extension
                vz = R * np.sqrt( 1 - n[2]**2 )

                # Vertical extension in number of cells
                dk = np.floor(vz / dz).astype(int)

                # Loop over the indices of horizontal levels
                for k in range( max(kc-dk,0), min(kc+dk+2,nz) ) :

                    # Corresponding z value
                    zi = z0 + k * dz + dz / 2

                    # Computes the points of intersection of the fracture with horizontal plane
                    intersect = self._disk_zplane_intersect([xc, yc, zc], n, R, zi)

                    # If there is an intersection
                    if np.isfinite( intersect[0] ) :

                        # Get matrix indices from x and y coordinates
                        #i1, j1 = xy2ij(intersect[0], intersect[1])
                        #i2, j2 = xy2ij(intersect[2], intersect[3])
                        i1, j1 = self.grid.get_i(intersect[0]), self.grid.get_j(intersect[1])
                        i2, j2 = self.grid.get_i(intersect[2]), self.grid.get_j(intersect[3])

                        # Rasterize the line
                        self._rst2d(raster_fractures[:,:,k], i1, i2, j1, j2)

            elif nfx**2 < (nfz**2 + nfy**2) : # subhorizontal case 1
                # Horizontal x index of the center of the fracture
                ic = ( (xc - x0) / dx  ).astype(int)

                # Projected x extension
                vx = R * np.sqrt( 1 - n[0]**2 )
                # x extension in number of cells
                di = np.floor(vx / dx).astype(int)

                # Loop over the indices of horizontal levels
                for i in range( max(ic-di,0), min(ic+di+2,nx) ) :

                    # Corresponding x value
                    xi = x0 + i * dx + dx / 2

                    # Computes the points of intersection of the fracture with vertical x plane
                    intersect = self._disk_xplane_intersect([xc, yc, zc], n, R, xi)

                    # If there is an intersection
                    if np.isfinite( intersect[0] ) :

                        # Get matrix indices from y and z coordinates
                        #j1, k1 = self.yz2jk(intersect[0], intersect[1])
                        #j2, k2 = self.yz2jk(intersect[2], intersect[3])
                        j1, k1 = self.grid.get_j(intersect[0]), self.grid.get_k(intersect[1])
                        j2, k2 = self.grid.get_j(intersect[2]), self.grid.get_k(intersect[3])

                        # Rasterize the line
                        self._rst2d(raster_fractures[i,:,:], j1, j2, k1, k2)

            else : # This may not be necessary
                # Horizontal y index of the center of the fracture
                yc = ( (yc - y0) / dx  ).astype(int)
                # Projected y extension
                vy = R * np.sqrt( 1 - n[1]**2 )
                # x extension in number of cells
                dj = np.floor(vy / dy).astype(int)

                # Loop over the indices of horizontal levels
                for j in range( max(jc-dj,0), min(jc+dj+2,nx) ) :

                    # Corresponding x value
                    yi = y0 + j * dy + dy / 2

                    # Computes the points of intersection of the fracture with vertical x plane
                    intersect = self._disk_yplane_intersect([xc, yc, zc], n, R, yi)

                    # If there is an intersection
                    if np.isfinite( intersect[0] ) :

                        # Get matrix indices from y and z coordinates
                        #i1, k1 = self.xz2ik(intersect[0], intersect[1])
                        #i2, k2 = self.xz2ik(intersect[2], intersect[3])
                        i1, k1 = self.grid.get_i(intersect[0]), self.grid.get_k(intersect[1])
                        i2, k2 = self.grid.get_i(intersect[2]), self.grid.get_k(intersect[3])

                        # Rasterize the line
                        self._rst2d(raster_fractures[:,j,:], i1, i2, k1, k2)

            #raster_frac = np.swapaxes(raster_fractures,0,2)
            self.data['fractures']['data'] = raster_fractures
            self.data['fractures']['img']  = np.transpose(raster_fractures, (1,0,2))
            self.data['fractures']['mode'] = 'random'
        return raster_fractures

    def compute_stats_on_data(self, data_key):
        """
        Computes statistics (occurrence, frequency and area) of the geologic data.

        Parameters
        ----------
        data_key : str
            Type of data : 'geology', 'faults', 'fractures'.

        Examples
        --------
        >>> geol.compute_stats_on_data("geology")
        >>> print(geol.data["geology"]["stats"])
        """
        if data_key in ['geology', 'faults', 'fractures']:
            data = self.data[data_key]["data"]
            unique, counts = np.unique(data, return_counts=True)
            occurrence, frequency, area = [], [], []
            for nbr in counts:
                occurrence.append(nbr)
                frequency.append(100*nbr/(self.grid.nx*self.grid.ny*self.grid.nz))
                area.append(nbr*self.grid.dx*self.grid.dy*self.grid.dz)
            self.data[data_key]['stats'] = {'ID':unique, 'occurrence':occurrence, 'frequency':frequency, 'area':area}
        else:
            self.data[data_key]['stats'] = np.nan
        return None


    def _show_fractures_stats(self):
        """

        """
        return None


    def show(self, data="geology", cmap='gray_r'):
        """
        Shows data of the geology manager.

        Parameters
        ----------
        data : str || array, optional
            Data to show : 'geology', 'topography', 'orientationx', 'orientationy' 'faults' or 'fractures'.
            By default, all data are showed.
        cmap : str, optional
            Color map, 'gray_r' by default.

        Examples
        --------
        >>> geol.show()
        """
        fig, ax1 = plt.subplots()
        d = self.data[data]['img']
        """
        origin = None
        if 'img' in self.data[data]:
            d = self.data[data]['img']
        elif data in ['topography', 'surface', 'orientationx', 'orientationy']:
            d = np.transpose(self.data[data]['data'], (1,0))
        else:
            d = np.transpose(self.data[data]['data'], (1,0,2)) # because imshow MxN
        if self.data[data]['mode'] in ['gslib', 'csv']:
            origin="lower"
        """
        im = ax1.imshow(d, extent=self.grid.extent, cmap=cmap)
        fig.colorbar(im, ax=ax1)
        plt.title(data)
        plt.show()
        return None
