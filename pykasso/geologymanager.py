"""
TODO :
- _check_geological_IDs ???
- _compute_surface ???
- Show data ?
- Show fractures stats ?
"""
from .functions import opendatafile, loadpoints
from .fracture  import Fracture
from tqdm import tqdm

import sys
import math
import mpmath
import numpy  as np

class GeologyManager():
    """
    Create a geology manager : a class for managing geology, faults and fractures.

    Parameters
    ----------
    grid : Grid()
        GeologyManager() class needs a Grid() object as argument.
    """

    def __init__(self, grid):
        #data model = {key : {data : var, stats : var, mode : var}}
        self.data  = {}   
        self.grid  = grid

    def set_data_null(self, data_key):
        """
        Set data to 'null' for an indicated data key.

        Parameters
        ----------
        data_key : str
            Type of data : 'geology', 'faults', 'fractures' or 'karst'.
        """
        self.data[data_key] = {}
        self.data[data_key]['data']  = np.zeros((self.grid.nz, self.grid.ny, self.grid.nx))
        occurence = self.grid.nz * self.grid.ny * self.grid.nx
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
        data_key : str
            Type of data : 'geology', 'faults', 'fractures' or 'karst'.
        datafile_location : str
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
        text = opendatafile(datafile_location)

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
        if len(data_lines) != (self.grid.nx*self.grid.ny*self.grid.nz):
            sys.exit("- set_data() - Error : Dimensions declared does not match with gslib file's size.")

        # Get values from text files
        data = np.zeros(len(data_lines))
        try:
            for data_line, k in zip(data_lines, range(len(data_lines))):
                data[k] = data_line.strip().split()[0]
        except:
            print("- set_data() - Unexpected error:", sys.exc_info()[0])
            raise

        return data.reshape((self.grid.nz, self.grid.ny, self.grid.nx))

    def _compute_stats(self, data_key):
        """
        Compute statistics on the imported data.
        """
        stats = {}

        # Count occurencies
        for z in range(self.grid.nz):
            for y in range(self.grid.ny):
                for x in range(self.grid.nx):
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
            frequency.append(100*stats[k]/(self.grid.nx*self.grid.ny*self.grid.nz))
            volume.append(stats[k]*self.grid.dx*self.grid.dy*self.grid.dz)
        return {'ID':ID, 'occurence':occurence, 'frequency':frequency, 'volume':volume}

    # def _check_geological_IDs(self, IDs):
        # for ID in IDs:
            # if ID not in self.data["geology"]["stats"]["ID"]:
                # sys.exit("Detected geological IDs do not match which indicated ones.")
        # self.geological_IDs = IDs
        # return None

    # def _compute_surface(self):
        # """
        # Compute the surface IDs and the digital elevation model.
        # """
        # if self.data["geology"]["mode"] is "null":
            # self.DEM         = np.full((self.grid.ny, self.grid.nx), self.grid.nz * self.grid.dz)
            # self.surface_IDs = np.full((self.grid.ny, self.grid.nx), 0)
        # elif self.data["geology"]["mode"] is "import":
            # self.DEM         = np.zeros((self.grid.ny,self.grid.nx))
            # self.surface_IDs = np.zeros((self.grid.ny,self.grid.nx))
            # for y in range(self.grid.ny):
                # for x in range(self.grid.nx):
                    # for z in range(self.grid.nz):
                        # if self.data["geology"]["data"][z][y][x] == self.geological_IDs[0]:
                            # if z == 0:
                                # self.surface_IDs[y][x] = self.geological_IDs[0]
                                # break
                        # else:
                            # self.DEM[y][x] = z + 1
                            # self.surface_IDs[y][x] = self.data["geology"]["data"][z][y][x]
            # self.DEM = self.DEM * self.grid.dz + self.grid.z0
        # return None

    def generate_fractures(self, fractures_densities, alpha, fractures_min_orientation, fractures_max_orientation, fractures_min_dip, fractures_max_dip, fractures_min_length, fractures_max_length):
        """
        Generate fractures as Fracture() objects.

        Parameters
        ----------
        fractures_densities : list
            Fracture densities for each fracture family.
        alpha : list
            Degree of power law for each fracture family.
        fractures_min_orientation : list
            Fractures minimum orientation for each fracture family.
        fractures_max_orientation : list
            Fractures maximum orientation for each fracture family.
        fractures_min_dip : list
            Fractures minimum dip for each fracture family.
        fractures_max_dip : list
            Fractures maximum dip for each fracture family.
        fractures_min_length : list
            The minimum lenght of the fractures for each fracture family.
        fractures_max_length : list
            The maximum lenght of the fractures for each fracture family.
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

                # FRACTURE ORIENTATION/DIP from von Mises distribution
                frac_orientation = np.random.vonmises(orientation_mean_angle, orientation_kappa)
                if (self.grid.nz > 1):
                    frac_dip = np.random.vonmises(dip_mean_angle, dip_kappa)
                else:
                    frac_dip = math.radians(90)
                
                # Calculate normal Vector
                x = np.sin(frac_orientation) * np.cos(frac_dip)
                y = np.cos(frac_orientation) * np.cos(frac_dip)
                z = np.sin(frac_dip)
                a = y
                b = -x
                c = -z

                # Store calculated fracture
                self.fractures.append(Fracture(fracture_id, frac_family, [xm,ym,zm], radius, frac_orientation, frac_dip, [a,b,c]))
                fracture_id += 1
        return None

    def float_eq(self, a, b, tolerance=1e-5):
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
        
    def unit_intersect(self, n, d):
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
        
        if not self.float_eq( np.linalg.norm(n), 1) :
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
        
    def disk_zplane_intersect(self, center, n, R, zi):
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
                xint = self.unit_intersect(n2, b)
                
                x1 = tau * xint[0] + xc
                x2 = tau * xint[1] + xc
                y1 = tau * xint[2] + yc
                y2 = tau * xint[3] + yc
                
        return np.array([x1, y1, x2, y2])
        
    def disk_xplane_intersect(self, center, n, R, xi):
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
                xint = self.unit_intersect(n2, b)
                
                y1 = tau * xint[0] + yc
                y2 = tau * xint[1] + yc
                z1 = tau * xint[2] + zc
                z2 = tau * xint[3] + zc
                
        return np.array([y1, z1, y2, z2])
        
    def disk_yplane_intersect(self, center, n, R, yi):
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
                xint = self.unit_intersect(n2, b)    
                
                x1 = tau * xint[0] + xc
                x2 = tau * xint[1] + xc
                z1 = tau * xint[2] + zc
                z2 = tau * xint[3] + zc    
                
        return np.array([x1, z1, x2, z2])
        
    def rst2d(self, m, xs, xe, ys, ye):
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
        
        nz, ny, nx = m.shape 
        
        if xe < xs : # Ensuires dx always positive
            xe, xs = xs, xe
            ye, ys = ys, ye
        
        dx = xe - xs
        dy = ye - ys
        
        if xs > nx : # The line is entirely out of the grid
            return
            
        if (dx + np.abs(dy)) == 0 : # Case with a single pixel
            if (ys >= 0) and (ys < ny) :
                m[ys, xs] = 1
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
            
        m[j,i] = 1
        return
        
    def rasterize_fracture_network(self):
        """
        Rasterize a set of fractures on a 3D grid.
        
        Returns
        -------        
        raster_fractures : 3D numpy array
            with value 1 where the grid is touched by a fracture
        """
        
        dx, dy, dz = self.grid.dx, self.grid.dy, self.grid.dz
        nx, ny, nz = self.grid.nx, self.grid.ny, self.grid.nz
        x0, y0, z0 = self.grid.x0, self.grid.y0, self.grid.z0    

        # Creates empty array to store raster of fracture indicators
        raster_fractures = np.zeros( (nz, ny, nx) )

        # Loop over the fractures
        for f in tqdm(self.fractures, desc="Rasterizing"):

            # Get fracture geometry
            xc, yc, zc = f.get_position()
            nfx, nfy, nfz = f.get_normal() # Could provide directly an array ### FM : ??? => n = f.get_normal()
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
                    intersect = self.disk_zplane_intersect([xc, yc, zc], n, R, zi)

                    # If there is an intersection
                    if np.isfinite( intersect[0] ) :

                        # Get matrix indices from x and y coordinates
                        #i1, j1 = xy2ij(intersect[0], intersect[1])
                        #i2, j2 = xy2ij(intersect[2], intersect[3])
                        i1, j1 = self.grid.get_i(intersect[0]), self.grid.get_j(intersect[1])
                        i2, j2 = self.grid.get_i(intersect[2]), self.grid.get_j(intersect[3])

                        # Rasterize the line
                        self.rst2d(raster_fractures[k,:,:], i1, i2, j1, j2)

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
                    intersect = self.disk_xplane_intersect([xc, yc, zc], n, R, xi)

                    # If there is an intersection
                    if np.isfinite( intersect[0] ) :

                        # Get matrix indices from y and z coordinates
                        #j1, k1 = self.yz2jk(intersect[0], intersect[1])
                        #j2, k2 = self.yz2jk(intersect[2], intersect[3])
                        j1, k1 = self.grid.get_j(intersect[0]), self.grid.get_k(intersect[1])
                        j2, k2 = self.grid.get_j(intersect[2]), self.grid.get_k(intersect[3])

                        # Rasterize the line
                        self.rst2d(raster_fractures[:,:,:], j1, j2, k1, k2)

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
                    intersect = self.disk_yplane_intersect([xc, yc, zc], n, R, yi)

                    # If there is an intersection
                    if np.isfinite( intersect[0] ) :

                        # Get matrix indices from y and z coordinates
                        #i1, k1 = self.xz2ik(intersect[0], intersect[1])
                        #i2, k2 = self.xz2ik(intersect[2], intersect[3])
                        i1, k1 = self.grid.get_i(intersect[0]), self.grid.get_k(intersect[1])
                        i2, k2 = self.grid.get_i(intersect[2]), self.grid.get_k(intersect[3])

                        # Rasterize the line
                        self.rst2d(raster_fractures[:,:,:], i1, i2, k1, k2)
        return raster_fractures

    # def _generate_fractures_maps(self):
        # """
        # Draw the fractures map according to generated fracturation.
        # """
        # fractures_maps = np.zeros((len(self.fractures)+1,self.grid.ny,self.grid.nx))

        # for frac_family in self.fractures:
            # for fracture in self.fractures[frac_family]["fractures"]:
                # x0, y0, x1, y1 = fracture.position
                # d = np.sqrt((x1-x0)**2 + (y1-y0)**2)
                # D = np.sqrt((self.grid.dx/2)**2 + (self.grid.dy/2)**2)
                # n = int(d / D)
                # xs = np.linspace(x0, x1, n)
                # ys = np.linspace(y0, y1, n)

                # for (x,y) in zip(xs, ys):
                    # X = int(math.ceil((x - self.grid.x0 - self.grid.dx/2) / self.grid.dx))
                    # Y = int(math.ceil((y - self.grid.y0 - self.grid.dy/2) / self.grid.dy))
                    # fractures_maps[frac_family][Y][X] = 1

            # self.fractures[frac_family]['frac_map'] = fractures_maps[frac_family]
        # self.data['fractures']['data'] = (sum([fractures_maps[i] for i in range(len(self.fractures))]) > 0).astype(int)
        # self.data['fractures']['mode'] = 'generate'
        # return None

    def show_fractures_stats(self):
    # if fractures none
        radius      = {}
        orientation = {}
        dip         = {}
        for family in range(self.fractures[-1].get_family()+1):
            radius[family]      = []
            orientation[family] = []
            dip[family]         = []
            for frac in self.fractures:
                if (family == frac.get_family()):
                    radius[family].append(frac.get_radius())
                    orientation[family].append(frac.get_orientation())
                    dip[family].append(frac.get_dip())
        print(radius)
        print(orientation)
        print(dip)
        return None
        
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

    # def show(self, frac_family=None):
        # """
        # Show the geology, faults, fractures and initial karst network maps.
        # """
        # grid = self._get_pyvista_grid()
        # grid.cell_arrays["values"] = self.data["geology"]["data"].flatten(order="F")

        # p = pv.Plotter(shape=(2, 2))

        # p.subplot(0, 0)
        # p.add_text("Geology", font_size=30)
        # p.add_mesh(grid)

        # p.subplot(0, 1)
        # p.add_text("Faults", font_size=30)
        # p.add_mesh(pv.Cube(), show_edges=True, color="tan")

        # p.subplot(1, 0)
        # p.add_text("Fractures", font_size=30)
        # sphere = pv.Sphere()
        # p.add_mesh(sphere, scalars=sphere.points[:, 2])
        # p.add_scalar_bar("Z")
        #plotter.add_axes()
        # p.add_axes(interactive=True)

        # p.subplot(1, 1)
        # p.add_text("Karst network", font_size=30)
        # p.add_mesh(pv.Cone(), color="g", show_edges=True)
        # p.show_bounds(all_edges=True)

        #Display the window
        # p.show()

def _extents(f):
    """
    Little private function to correctly locate data on plots.
    """
    delta = f[1] - f[0]
    return [f[0] - delta/2, f[-1] + delta/2]
