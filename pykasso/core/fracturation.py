"""
TODO
"""

import math
import mpmath
import numpy as np
from numba import njit # demander à val

# Condition
import pyvista as pv

def generate_fractures(grid, densities, alpha, min_orientation, max_orientation, min_dip, max_dip, min_length, max_length):
    """
    TODO
    Generates fractures as Fracture instances according to the parameters.

    Parameters
    ----------
    densities : array
        Fracture densities for each fracture family.
    alpha : array
        Degree of power law for each fracture family.
    min_orientation : array
        Fractures minimum orientation for each fracture family.
    max_orientation : array
        Fractures maximum orientation for each fracture family.
    min_dip : array
        Fractures minimum dip for each fracture family.
    max_dip : array
        Fractures maximum dip for each fracture family.
    min_length : array
        The minimum length of the fractures for each fracture family.
    max_length : array
        The maximum length of the fractures for each fracture family.

    Examples
    --------
    >>> geol._generate_fractures(densities = [0.00005,0.0001],
                                alpha = [2, 2],
                                min_orientation = [340, 70],
                                max_orientation = [20, 110],
                                min_dip = [80, 40],
                                max_dip = [90, 50],
                                min_length : [ 100, 100],
                                max_length : [4000, 4000])
    >>> print(geol.fractures)
    """

    fractures = []
    fracture_id = 0

    ## Redefine fracturation domain
    lenmax = max(max_length)
    xd0, xd1, yd0, yd1, zd0, zd1 = grid.xmin, grid.xmax, grid.ymin, grid.ymax, grid.zmin, grid.zmax
    Lx, Ly, Lz = xd1-xd0, yd1-yd0, zd1-zd0

    shiftx = min(Lx/2,lenmax/2)
    shifty = min(Ly/2,lenmax/2)
    shiftz = min(Lz/2,lenmax/2)

    Lex = 2 * shiftx + Lx
    Ley = 2 * shifty + Ly
    Lez = 2 * shiftz + Lz

    area = Lex * Ley
    xmin = grid.xlimits[0] - shiftx
    ymin = grid.ylimits[0] - shifty
    zmin = grid.zlimits[0] - shiftz

    ## Total numbers of fractures for each family
    fractures_numbers = np.array(densities) * area

    ## Calculate fractures in each family
    for frac_family in range(len(densities)):

        # Define all the constants required to randomly draw the length in distribution
        palpha = (1-alpha[frac_family])
        invpalpha = 1/palpha
        fmina = min_length[frac_family]**palpha
        frangea = max_length[frac_family]**palpha - fmina

        # Define all the constants required for the orientation distribution
        orientation_min = min_orientation[frac_family]
        orientation_max = max_orientation[frac_family]
        if orientation_min > orientation_max:
            orientation_min = orientation_min - 360

        orientation_mean_angle = (orientation_max + orientation_min)  / 2
        orientation_std_angle  = (orientation_max - orientation_mean_angle) / 3
        orientation_mean_angle = math.radians(orientation_mean_angle)
        orientation_std_angle  = math.radians(orientation_std_angle)

        # Define all the constants required for the dip distribution
        if (grid.nz > 1):
            dip_min = min_dip[frac_family]
            dip_max = max_dip[frac_family]
            if dip_min > dip_max:
                dip_min = dip_min - 360

            dip_mean_angle = (dip_max + dip_min)  / 2
            dip_std_angle  = (dip_max - dip_mean_angle) / 3
            dip_mean_angle = math.radians(dip_mean_angle)
            dip_std_angle  = math.radians(dip_std_angle)

        # Computes the kappa value for the Von Mises orientation distribution
        orientation_func  = lambda k : orientation_std_angle**2 - 1 + mpmath.besseli(1,k) / mpmath.besseli(0,k)
        orientation_kappa = mpmath.findroot(orientation_func,1)

        if (grid.nz > 1):
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
            if (grid.nz > 1):
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
            fractures.append(Fracture(fracture_id, frac_family, [xm,ym,zm], radius, frac_orientation, frac_dip, [a,b,c]))
            fracture_id += 1
    return fractures

####################
### VOXELIZATION ###
####################

def voxelize_fractures(grid, fractures, method):
    """
    TODO
    """
    if method == 'python':
        fracs_array = _voxelize_fractures_python(grid, fractures)
    elif method == 'pyvista':
        fracs_array = _voxelize_fractures_pyvista(grid, fractures)
    elif method == 'c':
        pass
    else:
        print('TODO')

    return fracs_array

#############################################################################################################################
### pyvista ###
###############

def _voxelize_fractures_pyvista(grid, fractures):
    """
    TODO
    """

    ##############################
    ### 1 - Construct the rays ###
    ##############################
    slice_x = slice(0, grid.nx)
    slice_y = slice(0, grid.ny) 
    slice_z = slice(0, grid.nz)
    val_s, val_f = 0, -1
    SLICES = [ # xz, zy, yx
        [[slice_x, val_s  , slice_z], [slice_x, val_f  , slice_z]],
        [[val_s  , slice_y, slice_z], [val_f  , slice_y, slice_z]],
        [[slice_x, slice_y,   val_s], [slice_x, slice_y,   val_f]]
    ]
    DIMS = ['X', 'Y', 'Z']
    POINTS_START = []
    POINTS_FINAL = []
    # RAYS = []

    for (slice_start, slice_final) in SLICES:

        points_start = []
        points_final = []

        for dim in DIMS:
            i_s, j_s, k_s = slice_start
            coordinates = grid._get_property(dim)[i_s, j_s, k_s]
            coordinates = coordinates.flatten(order='F')
            points_start.append(coordinates)

            i_f, j_f, k_f = slice_final
            coordinates = grid._get_property(dim)[i_f, j_f, k_f]
            coordinates = coordinates.flatten(order='F')
            points_final.append(coordinates)

        points_start = [list(tpl) for tpl in zip(*points_start)]
        points_final = [list(tpl) for tpl in zip(*points_final)]

        POINTS_START.append(points_start)
        POINTS_FINAL.append(points_final)

        # for (start, final) in zip(points_start, points_final):
        #     RAYS.append(pv.Line(start, final))

    # flat the lists
    POINTS_START = [item for sublist in POINTS_START for item in sublist]
    POINTS_FINAL = [item for sublist in POINTS_FINAL for item in sublist]
    # convert in numpy arrays
    POINTS_START = np.array(POINTS_START)
    POINTS_FINAL = np.array(POINTS_FINAL)

    #####################################################
    ### 2 - Transform fractures into pyvista polygons ###
    #####################################################
    POLYGONS = []
    polygon_definition = 8 # TODO - Memory settings
    for fracture in fractures:
        x, y, z = fracture.get_position()
        a, b, c = fracture.get_normal()
        rad     = fracture.radius
        POLYGONS.append(pv.Polygon(center=(x, y, z), radius=rad, normal=(a, b, c), n_sides=polygon_definition))

    ###############################
    ### 3 - Intersection points ###
    ###############################
    POINTS = []
    for polygon in POLYGONS:
        for (start, final) in zip(POINTS_START, POINTS_FINAL):
            points, ind = polygon.ray_trace(start, final)
            if points.size != 0:
                POINTS.append(points)

    ############################
    ### 4 - Retrieve indices ###
    ############################
    INDICES = []
    for p in POINTS:
        x, y, z = p[0]
        i, j, k = grid.get_i(x), grid.get_j(y), grid.get_k(z)
        INDICES.append([i, j, k])

    # Remove duplicates
    INDICES = np.array(INDICES)
    INDICES = np.unique(INDICES, axis=0)
    i, j, k = zip(*INDICES)

    ############################
    ### 5 - Fracturation map ###
    ############################
    fracs_array = np.zeros((grid.nx, grid.ny, grid.nz))
    fracs_array[i, j, k] = 1

    return fracs_array

#############################################################################################################################
### python ###
##############

def _float_eq(a, b, tolerance=1e-5):
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

def _unit_intersect(n, d):
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

    if not _float_eq( np.linalg.norm(n), 1) :
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

def _disk_zplane_intersect(center, n, R, zi):
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
            xint = _unit_intersect(n2, b)

            x1 = tau * xint[0] + xc
            x2 = tau * xint[1] + xc
            y1 = tau * xint[2] + yc
            y2 = tau * xint[3] + yc

    return np.array([x1, y1, x2, y2])

def _disk_xplane_intersect(center, n, R, xi):
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
            xint = _unit_intersect(n2, b)

            y1 = tau * xint[0] + yc
            y2 = tau * xint[1] + yc
            z1 = tau * xint[2] + zc
            z2 = tau * xint[3] + zc

    return np.array([y1, z1, y2, z2])

def _disk_yplane_intersect(center, n, R, yi):
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
            xint = _unit_intersect(n2, b)

            x1 = tau * xint[0] + xc
            x2 = tau * xint[1] + xc
            z1 = tau * xint[2] + zc
            z2 = tau * xint[3] + zc

    return np.array([x1, z1, x2, z2])

def _rst2d(m, xs, xe, ys, ye):
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

    if xs >= nx : # The line is entirely out of the grid
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

def _voxelize_fractures_python(grid, fractures):
    """
    Rasterizes a set of fractures on a 3D grid.

    Returns
    -------
    raster_fractures : 3D numpy array
        With value 1 where the grid is touched by a fracture
    """

    dx, dy, dz = grid.dx, grid.dy, grid.dz
    nx, ny, nz = grid.nx, grid.ny, grid.nz
    x0, y0, z0 = grid.x0, grid.y0, grid.z0

    # Creates empty array to store raster of fracture indicators
    raster_fractures = np.zeros( (nx, ny, nz) , dtype=np.int_)

    # Loop over the fractures
    # for f in tqdm(fractures, desc="Rasterizing"):
    for f in fractures:

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
                intersect = _disk_zplane_intersect([xc, yc, zc], n, R, zi)

                # If there is an intersection
                if np.isfinite( intersect[0] ) :

                    # Get matrix indices from x and y coordinates
                    #i1, j1 = xy2ij(intersect[0], intersect[1])
                    #i2, j2 = xy2ij(intersect[2], intersect[3])
                    i1, j1 = grid.get_i(intersect[0]), grid.get_j(intersect[1])
                    i2, j2 = grid.get_i(intersect[2]), grid.get_j(intersect[3])

                    # Rasterize the line
                    _rst2d(raster_fractures[:,:,k], i1, i2, j1, j2)

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
                intersect = _disk_xplane_intersect([xc, yc, zc], n, R, xi)

                # If there is an intersection
                if np.isfinite( intersect[0] ) :

                    # Get matrix indices from y and z coordinates
                    #j1, k1 = yz2jk(intersect[0], intersect[1])
                    #j2, k2 = yz2jk(intersect[2], intersect[3])
                    j1, k1 = grid.get_j(intersect[0]), grid.get_k(intersect[1])
                    j2, k2 = grid.get_j(intersect[2]), grid.get_k(intersect[3])

                    # Rasterize the line
                    _rst2d(raster_fractures[i,:,:], j1, j2, k1, k2)

        else : # This may not be necessary
            # Horizontal y index of the center of the fracture
            yc = ( (yc - y0) / dx  ).astype(int)
            # Projected y extension
            vy = R * np.sqrt( 1 - n[1]**2 )
            # x extension in number of cells
            dj = np.floor(vy / dy).astype(int)

            # Loop over the indices of horizontal levels
            # for j in range( max(jc-dj,0), min(jc+dj+2,nx) ) : # TODO - correction correcte ?
            for j in range( max(yc-dj,0), min(yc+dj+2,nx) ) :

                # Corresponding x value
                yi = y0 + j * dy + dy / 2

                # Computes the points of intersection of the fracture with vertical x plane
                intersect = _disk_yplane_intersect([xc, yc, zc], n, R, yi)

                # If there is an intersection
                if np.isfinite( intersect[0] ) :

                    # Get matrix indices from y and z coordinates
                    #i1, k1 = xz2ik(intersect[0], intersect[1])
                    #i2, k2 = xz2ik(intersect[2], intersect[3])
                    i1, k1 = grid.get_i(intersect[0]), grid.get_k(intersect[1])
                    i2, k2 = grid.get_i(intersect[2]), grid.get_k(intersect[3])

                    # Rasterize the line
                    _rst2d(raster_fractures[:,j,:], i1, i2, k1, k2)

        #raster_frac = np.swapaxes(raster_fractures,0,2)
    return raster_fractures

###############
### CLASSES ###
###############

class Fracture():
    """
    Class modeling fractures as objects.
    """

    def __init__(self, ID, family, position, radius, orientation, dip, normal):
        """
        Creates a fracture according to the parameters.

        Parameters
        ----------
        ID : int
            Fracture id.
        family : int
            Fracture family id.
        position : array
            Position of the center of the fracture [x, y, z].
        radius : float
            Radius of the fracture.
        orientation : float
            Orientation of the fracture.
        dip : float
            Dip of the fracture.
        normal : array
            Normal vector of the fracture [a, b, c].

        Examples
        --------
        >>> frac = pk.Fracture(0, 0, [10,10,10], 5, 180, 90, [1,0,0])
        """
        self.ID          = ID
        self.family      = family
        self.position    = position
        self.radius      = radius
        self.orientation = orientation
        self.dip         = dip
        self.normal      = normal

    def __repr__(self):
        return '[id:{}, fam.:{}, x:{}, y:{}, z:{}, rad:{}, or.:{}, dip.:{}, n:({},{},{})] \n'.format(self.ID,self.family, round(self.position[0],2), round(self.position[1],2), round(self.position[2],2), round(self.radius,2), round(self.orientation,2), round(self.dip,2), round(self.normal[0],2), round(self.normal[1],2), round(self.normal[2],2))

    def get_ID(self):
        """
        Returns the ID of the fracture.

        Examples
        --------
        >>> i = frac.get_ID()
        """
        return self.ID

    def get_family(self):
        """
        Returns the family ID of the fracture.

        Examples
        --------
        >>> family = frac.get_family()
        """
        return self.family

    def get_position(self):
        """
        Returns an array with the (x, y, z) coordinates of the center of the fracture.

        Examples
        --------
        >>> x, y, z = frac.get_position()
        """
        return self.position

    def get_radius(self):
        """
        Returns the radius of the fracture.

        Examples
        --------
        >>> rad = frac.get_radius()
        """
        return self.radius

    def get_orientation(self):
        """
        Returns the orientation of the fracture.

        Examples
        --------
        >>> orien = frac.get_orientation()
        """
        return self.orientation

    def get_dip(self):
        """
        Returns the dip of the fracture.

        Examples
        --------
        >>> dip = frac.get_dip()
        """
        return self.dip

    def get_normal(self):
        """
        Returns an array with the (a, b, c) vector component of the normal of the fracture.

        Examples
        --------
        >>> a, b, c = frac.get_normal()
        """
        return self.normal