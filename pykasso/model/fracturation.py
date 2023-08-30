"""
This module contains fonctions designed for fracture generation and
discretization.
"""

### External dependencies
import mpmath
import numpy as np
import pandas as pd

### Local dependencies
from pykasso.model.geologic_features import Volume

### Typing
from pykasso._typing import Grid, RandomNumberGenerator

## TODO
# equation _solve_kappa, -1 ???
# distribution vonmises arguments : mean, std / loc, scale
# (https://numpy.org/doc/stable/reference/random/generated/numpy.random.normal.html)


class Fractures(Volume):
    """Class modeling the fracturation model."""
    
    def __init__(self, rng: RandomNumberGenerator, *args, **kwargs):
        label = 'fractures'
        super().__init__(label, *args, **kwargs)
        
        self.rng = rng
        self.i = 0
        self.families = pd.DataFrame()
        self.fractures = pd.DataFrame()
        self.fractures_voxelized = {}
        self.model = 'superposition'
        self.models = {}
        
    def generate_fracture_family(self, name, grid, **settings):
        self.i = self.i + 1
        if 'cost' in settings:
            cost = settings['cost']
            del settings['cost']
        
        # Creates family
        frac_family = {
            'id': [self.i],
            'name': [name],
            'cost': [cost],
            # 'geology': pass # TODO - add geology
        }
        frac_family = pd.DataFrame(data=frac_family)
            
        # Generates fractures
        fractures = self.generate_fractures(grid, **settings)
        fractures.insert(0, 'family_id', self.i)
        
        # Updates tables
        self.fractures = pd.concat([self.fractures, fractures])
        self.families = pd.concat([self.families, frac_family])
        
        return None
     
    def generate_model(self, model: str, grid) -> None:
        """Constructs the model for fracturation according to selection."""
        frac_model = np.zeros_like(grid, dtype=np.float64)
        
        ### Voxelizes the fractures
        fractures_voxelized = {}
        for i in self.families['id']:
            fractures = self.fractures[self.fractures['family_id'] == i]
            if fractures.empty:
                fractures_voxelized[i] = np.zeros_like(grid.data_volume)
            else:
                fractures_voxelized[i] = voxelize_fractures(grid, fractures)
            
        ### 1 - Superposition of families
        if model == 'superposition':
            
            # Sorts fracture families by cost
            self.families = self.families.sort_values(['cost', 'id'])
            fractures_family_ids = self.families['id'].values
            
            # Updates the final array
            for family_id in fractures_family_ids:
                frac_model = np.where(
                    fractures_voxelized[family_id] == 1,
                    family_id,
                    frac_model
                )
        
        ### 2 - Sums of families
        elif model == 'sum':
            frac_model = sum([d for d in fractures_voxelized.values()])
        
        ### 3 - Binary grid giving presence/absence of fractures
        elif model == 'binary':
            frac_model = sum([d for d in fractures_voxelized.values()])
            frac_model = np.where(frac_model > 0, 1, 0)
        
        ### Selects final model
        self.data_volume = frac_model.copy()
        
        ######### TODO #########
        # # Constraints model with geology if provided
        # if 'geology' in kwargs:
        #     frac_model_geology = np.zeros_like(frac_model)
        #     for geologic_id in kwargs['geology']:
        #         frac_model_geology = np.where(
        #             Geology.data_volume == geologic_id,
        #             frac_model,
        #             frac_model_geology
        #         )
        #     self.fractures_families['model'] = frac_model_geology
                
        # data = self.fractures_families['model']
        ######### TODO #########
        
        # Sets the costs
        costs = pd.Series(self.families['cost'].values,
                          index=self.families['id']).to_dict()
        self._set_costs(costs)
        
        return None
                
    def generate_fractures(self, grid, density: float, alpha: float,
                           orientation: float, dip: float,
                           length: float,
                           orientation_distribution: str = 'vonmises',
                           dip_distribution: str = 'vonmises',
                           length_distribution: str = 'power'):
        """
        TODO
        """
        
        ######################
        ### INITIALIZATION ###
        ######################
        
        orientation_min = 0
        orientation_max = 360
        dip_min = 0
        dip_max = 360
        
        ### Defines if parameter is a range or an unique value
        if isinstance(orientation, (list)):
            orientation_min, orientation_max = orientation
        else:
            orientation_distribution = 'unique'
        
        if isinstance(dip, (list)):
            dip_min, dip_max = dip
        else:
            dip_distribution = 'unique'
        
        if isinstance(length, (list)):
            length_min, length_max = min(length), max(length)
        else:
            length_distribution = 'unique'
            length_max = length

        ### Redefine fracturation domain
        Lx = grid.xmax - grid.xmin
        Ly = grid.ymax - grid.ymin
        Lz = grid.zmax - grid.zmin

        shift_x = min(Lx / 2, length_max / 2)
        shift_y = min(Ly / 2, length_max / 2)
        shift_z = min(Lz / 2, length_max / 2)

        Lex = 2 * shift_x + Lx
        Ley = 2 * shift_y + Ly
        Lez = 2 * shift_z + Lz

        area = Lex * Ley
        xmin = grid.xmin - shift_x
        ymin = grid.ymin - shift_y
        zmin = grid.zmin - shift_z

        ### Total numbers of fractures
        fractures_numbers = np.array(density) * area

        ### Generate poisson number for each fracture family
        real_frac_number = self.rng.poisson(fractures_numbers)
        
        ###########################
        ### Calculate fractures ###
        ###########################

        ##### FRACTURE CENTER LOCATION

        # from triple uniform distribution
        xm = self._random_random(Lex, xmin, real_frac_number)
        ym = self._random_random(Ley, ymin, real_frac_number)
        zm = self._random_random(Lez, zmin, real_frac_number)

        ##### FRACTURE ORIENTATION
        
        if orientation_min > orientation_max:
            orientation_min = orientation_min - 360

        # No distribution case
        if orientation_distribution == 'unique':
            orientation_fractures = [orientation] * real_frac_number
        
        # Uniform distribution case
        elif orientation_distribution == 'uniform':
            orientation_fractures = self._random_uniform(orientation_min,
                                                         orientation_max,
                                                         real_frac_number)
        
        # Von Mises distribution case
        elif orientation_distribution == 'vonmises':
            # Computes ...
            orient_mean_angle = (orientation_max + orientation_min) / 2
            orient_std_angle = (orientation_max - orient_mean_angle) / 3
            orient_mean_angle_rad = np.radians(orient_mean_angle)
            orient_std_angle_rad = np.radians(orient_std_angle)

            # Computes the kappa value for the Von Mises orient. distribution
            orient_kappa = self._solve_kappa(orient_std_angle_rad)
            orientation_fractures = (
                self._random_vonmises(orient_mean_angle_rad,
                                      orient_kappa,
                                      real_frac_number)
            )
            orientation_fractures = np.degrees(orientation_fractures)
            
        ##### FRACTURE DIP
        
        if dip_min > dip_max:
            dip_min = dip_min - 360

        # No distribution case
        if dip_distribution == 'unique':
            dip_fractures = [dip] * real_frac_number
        
        # Uniform distribution case
        elif dip_distribution == 'uniform':
            dip_fractures = self._random_uniform(dip_min,
                                                 dip_max,
                                                 real_frac_number)
            
        # Von Mises distribution case
        elif dip_distribution == 'vonmises':
            dip_mean_angle = (dip_max + dip_min) / 2
            dip_std_angle = (dip_max - dip_mean_angle) / 3
            dip_mean_angle_rad = np.radians(dip_mean_angle)
            dip_std_angle_rad = np.radians(dip_std_angle)
            
            # Computes the kappa value for the Von Mises orient. distribution
            dip_kappa = self._solve_kappa(dip_std_angle_rad)
            dip_fractures = (
                self._random_vonmises(dip_mean_angle_rad,
                                      dip_kappa,
                                      real_frac_number)
            )
            dip_fractures = np.degrees(dip_fractures)

        ##### FRACTURE LENGHT

        # No distribution case
        if length_distribution == 'unique':
            radius = [length / 2] * real_frac_number

        # Uniform distribution case
        elif length_distribution == 'uniform':
            radius = self._random_uniform(length_min,
                                          length_max,
                                          real_frac_number)

        # Trucated power law distribution case
        elif length_distribution == 'power':
            frac_length = self._random_power(length_min, length_max,
                                             alpha, real_frac_number)
            radius = frac_length / 2
            
        ##### FRACTURE NORMAL VECTOR
        normal = calculate_normal(dip_fractures, orientation_fractures)

        #######################
        ### Returns results ###
        #######################
        
        data = {
            'x': xm,
            'y': ym,
            'z': zm,
            'radius': radius,
            'orientation': orientation_fractures,
            'dip': dip_fractures,
            'normal': normal
        }
        fractures = pd.DataFrame(data=data)

        return fractures

    def _random_random(self, a, b, n):
        """ """
        out = a * self.rng.random(size=n) + b
        return out
    
    def _random_uniform(self, min_value, max_value, n):
        """ """
        out = self.rng.uniform(min_value, max_value, size=n)
        return out
    
    def _solve_kappa(self, std_value):
        """ """
        func = lambda kappa: (std_value**2 - 1
                              + mpmath.besseli(1, kappa)
                              / mpmath.besseli(0, kappa))
        kappa = mpmath.findroot(func, 1)
        return kappa
    
    def _random_vonmises(self, mean_value, kappa, n):
        """ """
        out = self.rng.vonmises(mean_value, kappa, size=n)
        return out
    
    def _random_power(self, min_val, max_val, alpha, n):
        """ """
        palpha = (1 - alpha)
        invpalpha = 1 / palpha
        fmina = min_val**palpha
        frangea = max_val**palpha - fmina
        u = self.rng.random(size=n)
        out = (fmina + u * frangea)**invpalpha
        return out

 
def calculate_normal(dip, orientation):
    """"""
    orientation = np.radians(orientation)
    dip = np.radians(dip)
    x = np.sin(dip) * np.cos(np.radians(90) - orientation)
    y = np.sin(dip) * np.sin(np.radians(90) - orientation)
    z = np.cos(dip)
    a = y
    b = -x
    c = -z
    normal = list(zip(a, b, c))
    return normal


####################
### VOXELIZATION ###
####################


def _float_eq(a, b, tolerance: float = 1e-5):
    """
    Returns True if the difference between a and b is lower than tolerance.

    Parameters
    ----------
    a : _type_
        _description_
    b : _type_
        _description_
    tolerance : float, optional
        _description_, by default 1e-5

    Returns
    -------
    result : bool
    """
    if np.all(np.isnan(a)):
        if np.all(np.isnan(b)):
            return True
        else:
            return False
    return np.all(abs(a - b) < tolerance)


def _unit_intersect(n: np.ndarray, d: float):
    """
    Computes the intersection between a unit circle and a line on a 2D plane.
    The unit circle is centered on the origin (x=0, y=0) and has a radius of 1.
    The line is defined by parameters 'n' and 'd'.

    Parameters
    ----------
    n : numpy.ndarray of size 2
        The normal vector perpendicular to the line.
        Warning : The norm of n must be 1.
    d : float
        Defines the position of the line.
        It is a signed distance to the origin x=0, y=0.
        It can be positive or negative.

    Returns
    -------
    result : numpy.ndarray with 4 floating point values [X1, X2, Y1, Y2]
        The coordinates of the two points of intersection when it exists.
        Returns 4 np.NaN if there is no intersection.
    """
    nx, ny = n[0], n[1]  # For readibility

    if not _float_eq(np.linalg.norm(n), 1):
        print("WARNING - unitcintl function : norm of n must be equal to 1")

    if np.abs(d) > 1:  # Case with no intersection
        return np.array([np.NaN, np.NaN, np.NaN, np.NaN])

    sd = np.sqrt(1 - d**2)

    if ny != 0:
        X1 = d * nx + ny * sd
        X2 = d * nx - ny * sd
        Y1 = d * ny - nx * sd
        Y2 = d * ny + nx * sd
    else:
        X1 = X2 = d
        Y1 = sd
        Y2 = -sd

    return np.array([X1, X2, Y1, Y2])


def _disk_zplane_intersect(center: np.ndarray, n: np.ndarray,
                           R: float, zi: float):
    """
    Computes the intersection between a disk and a horizontal plane (constant
    z plane) in 3D. The disk is defined by three parameters :

    Parameters
    ----------
    center : numpy.ndarray of size 3
        The 3D coordinates of the center of the disk.
    n : numpy.ndarray of size 3
        The normal vector perpendicular to the disk.
        Warning : The norm of n must be 1.
    R : float
        Radius of the disk.
    zi : float
        Position of the horizontal plane along the z-axis.

    Returns
    -------
    result : np.ndarray with 4 floating point values [x1, y1, x2, y2]
        The coordinates of the two extremities of the intersection between the
        plane and disk if there is an intersection. Returns 4 np.NaN if there
        is no intersection.
    """

    xc, yc, zc = center[0], center[1], center[2]  # For readibility

    tau = R**2 - (zi - zc)**2

    if tau <= 0:  # Avoid computing square root of negative number
        # print("The plane does not touch the sphere")
        x1, y1, x2, y2 = np.NaN, np.NaN, np.NaN, np.NaN

    else:
        tau = np.sqrt(tau)
        b = n[2] * (zc - zi) / tau / np.sqrt(1 - n[2]**2)

        if b > 1:
            # print("The plane does not touch the disk")
            x1, y1, x2, y2 = np.NaN, np.NaN, np.NaN, np.NaN

        else:
            n2 = n.copy()[0:2]  # Projection on horizontal plane
            n2 /= np.linalg.norm(n2)
            xint = _unit_intersect(n2, b)

            x1 = tau * xint[0] + xc
            x2 = tau * xint[1] + xc
            y1 = tau * xint[2] + yc
            y2 = tau * xint[3] + yc

    return np.array([x1, y1, x2, y2])


def _disk_xplane_intersect(center: np.ndarray, n: np.ndarray,
                           R: float, xi: float):
    """
    Computes the intersection between a disk and a vertical plane (constant x
    plane) in 3D. The disk is defined by three parameters :

    Parameters
    ----------
    center : np.ndarray of size 3
        The 3D coordinates of the center of the disk.
    n : np.ndarray of size 3
        The normal vector perpendicular to the disk.
        Warning : The norm of n must be 1.
    R : float
        Radius of the disk.
    xi : float
        Position of the vertical plane along the x-axis.

    Returns
    -------
    result : np.ndarray with 4 floating point values [y1, z1, y2, z2]
        The coordinates of the two extremities of the intersection between the
        plane and disk if there is an intersection. Returns 4 np.NaN if there
        is no intersection.
    """

    xc, yc, zc = center[0], center[1], center[2]  # For readibility

    tau = R**2 - (xi - xc)**2

    if tau <= 0:  # Avoid computing square root of negative number
        # print("The plane does not touch the sphere")
        y1, z1, y2, z2 = np.NaN, np.NaN, np.NaN, np.NaN

    else:
        tau = np.sqrt(tau)
        b = n[0] * (xc - xi) / tau / np.sqrt(1 - n[0]**2)

        if b > 1:
            # print("The plane does not touch the disk")
            y1, z1, y2, z2 = np.NaN, np.NaN, np.NaN, np.NaN

        else:
            n2 = np.array([n[1], n[2]])  # Projection on vertical x plane
            n2 /= np.linalg.norm(n2)
            xint = _unit_intersect(n2, b)

            y1 = tau * xint[0] + yc
            y2 = tau * xint[1] + yc
            z1 = tau * xint[2] + zc
            z2 = tau * xint[3] + zc

    return np.array([y1, z1, y2, z2])


def _disk_yplane_intersect(center: np.ndarray, n: np.ndarray,
                           R: float, yi: float):
    """
    Computes the intersection between a disk and a vertical plane (constant y
    plane) in 3D. The disk is defined by three parameters :

    Parameters
    ----------
    center : np.ndarray of size 3
        The 3D coordinates of the center of the disk.
    n : np.ndarray of size 3
        The normal vector perpendicular to the disk.
        Warning : The norm of n must be 1.
    R : float
        Radius of the disk.
    yi : float
        Position of the vertical plane along the y-axis.

    Returns
    -------
    result : np.ndarray with 4 floating point values [x1, z1, x2, z2]
        The coordinates of the two extremities of the intersection between the
        plane and disk if there is an intersection. Returns 4 np.NaN if there
        is no intersection.
    """

    xc, yc, zc = center[0], center[1], center[2]  # For readibility

    tau = R**2 - (yi - yc)**2
    if tau <= 0:  # Avoid computing square root of negative number
        # print("The plane does not touch the sphere")
        x1, z1, x2, z2 = np.NaN, np.NaN, np.NaN, np.NaN

    else:
        tau = np.sqrt(tau)
        b = n[0] * (yc - yi) / tau / np.sqrt(1 - n[0]**2)

        if b > 1:
            # print("The plane does not touch the disk")
            x1, z1, x2, z2 = np.NaN, np.NaN, np.NaN, np.NaN

        else:
            n2 = np.array([n[0], n[2]])  # Projection on vertical x plane
            n2 /= np.linalg.norm(n2)
            xint = _unit_intersect(n2, b)

            x1 = tau * xint[0] + xc
            x2 = tau * xint[1] + xc
            z1 = tau * xint[2] + zc
            z2 = tau * xint[3] + zc

    return np.array([x1, z1, x2, z2])


def _rst2d(m: np.ndarray, xs: int, xe: int, ys: int, ye: int):
    """Rasterizes a line on a 2D plane.

    Parameters
    ----------
    m : numpy.ndarray of dim 2
    xs, xe, ys, ye : int
        Indices of positions on the grid m.
        Starting and ending location of the line.

    Returns
    -------
    m : numpy.ndarray of dim 2
        With value 1 where the grid is touched by the line.
    """

    nx, ny = m.shape

    if xe < xs:  # Ensuires dx always positive
        xe, xs = xs, xe
        ye, ys = ys, ye

    dx = xe - xs
    dy = ye - ys

    if xs >= nx:  # The line is entirely out of the grid
        return

    if (dx + np.abs(dy)) == 0:  # Case with a single pixel
        if (ys >= 0) and (ys < ny):
            m[xs, ys] = 1
        return

    if dx >= abs(dy):  # line with greater horizontal than vertical extension
        i = np.arange(max(0, xs), min(xe, nx))
        j = np.int32(np.round((i - xs) * dy / dx + ys))

        indx_ok = (j >= 0) & (j < ny)  # To crop the part that is out of grid
        i = i[indx_ok]
        j = j[indx_ok]

    else:
        if ye < ys:  # to ensure that arange works properly
            xe, xs = xs, xe
            ye, ys = ys, ye

        if ys > ny:
            return

        j = np.arange(max(ys, 0), min(ye, ny))
        i = np.int32(np.round((j - ys) * dx / dy + xs))

        indx_ok = (i >= 0) & (i < nx)  # To crop the part that is out of grid
        i = i[indx_ok]
        j = j[indx_ok]

    m[i, j] = 1
    return


def voxelize_fractures(grid, fractures):
    """Rasterizes a set of fractures on a 3D grid.
    
    Parameters
    ----------
    grid : _type_
        _description_
    fractures : _type_
        _description_

    Returns
    -------
    raster_fractures : numpy.ndarray of dim 3
        With value 1 where the grid is touched by a fracture.
    """

    dx, dy, dz = grid.dx, grid.dy, grid.dz
    nx, ny, nz = grid.nx, grid.ny, grid.nz
    x0, y0, z0 = grid.x0, grid.y0, grid.z0

    # Creates empty array to store raster of fracture indicators
    raster_fractures = np.zeros((nx, ny, nz), dtype=np.int_)

    # Loop over the fractures
    for (i, f) in fractures.iterrows():

        # Get fracture geometry
        xc, yc, zc = f['x'], f['y'], f['z']
        nfx, nfy, nfz = f['normal']
        nfx, nfy, nfz = float(nfx), float(nfy), float(nfz)
        n = [nfx, nfy, nfz]
        R = f['radius']

        # Test orientation of the fractures
        if nfz**2 <= (nfx**2 + nfy**2):  # subvertical case
            # Vertical index of the center of the fracture
            kc = ((zc - z0) / dz)
            kc = int(kc)

            # Projected vertical extension
            vz = R * np.sqrt(1 - n[2]**2)

            # Vertical extension in number of cells
            dk = np.floor(vz / dz).astype(int)

            # Loop over the indices of horizontal levels
            for k in range(max(kc - dk, 0), min(kc + dk + 2, nz)):

                # Corresponding z value
                zi = z0 + k * dz + dz / 2

                # Computes the points of intersection of the fracture with
                # horizontal plane
                intersect = _disk_zplane_intersect([xc, yc, zc], n, R, zi)

                # If there is an intersection
                if np.isfinite(intersect[0]):

                    # Get matrix indices from x and y coordinates
                    i1, j1 = grid.get_i(intersect[0]), grid.get_j(intersect[1])
                    i2, j2 = grid.get_i(intersect[2]), grid.get_j(intersect[3])

                    # Rasterize the line
                    _rst2d(raster_fractures[:, :, k], i1, i2, j1, j2)

        elif nfx**2 < (nfz**2 + nfy**2):  # subhorizontal case 1
            # Horizontal x index of the center of the fracture
            ic = ((xc - x0) / dx)
            ic = int(ic)

            # Projected x extension
            vx = R * np.sqrt(1 - n[0]**2)
            # x extension in number of cells
            di = np.floor(vx / dx).astype(int)

            # Loop over the indices of horizontal levels
            for i in range(max(ic - di, 0), min(ic + di + 2, nx)):

                # Corresponding x value
                xi = x0 + i * dx + dx / 2

                # Computes the points of intersection of the fracture with
                # vertical x plane
                intersect = _disk_xplane_intersect([xc, yc, zc], n, R, xi)

                # If there is an intersection
                if np.isfinite(intersect[0]):

                    # Get matrix indices from y and z coordinates
                    j1, k1 = grid.get_j(intersect[0]), grid.get_k(intersect[1])
                    j2, k2 = grid.get_j(intersect[2]), grid.get_k(intersect[3])

                    # Rasterize the line
                    _rst2d(raster_fractures[i, :, :], j1, j2, k1, k2)

        else:  # This may not be necessary
            # Horizontal y index of the center of the fracture
            jc = ((yc - y0) / dx)
            jc = int(jc)

            # Projected y extension
            vy = R * np.sqrt(1 - n[1]**2)
            # x extension in number of cells
            dj = np.floor(vy / dy).astype(int)

            # Loop over the indices of horizontal levels
            for j in range(max(jc - dj, 0), min(jc + dj + 2, nx)):

                # Corresponding x value
                yi = y0 + j * dy + dy / 2

                # Computes the points of intersection of the fracture with
                # vertical x plane
                intersect = _disk_yplane_intersect([xc, yc, zc], n, R, yi)

                # If there is an intersection
                if np.isfinite(intersect[0]):

                    # Get matrix indices from y and z coordinates
                    i1, k1 = grid.get_i(intersect[0]), grid.get_k(intersect[1])
                    i2, k2 = grid.get_i(intersect[2]), grid.get_k(intersect[3])

                    # Rasterize the line
                    _rst2d(raster_fractures[:, j, :], i1, i2, k1, k2)

        # raster_frac = np.swapaxes(raster_fractures,0,2)
    return raster_fractures
