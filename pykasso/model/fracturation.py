"""
This module contains fonctions designed for fracture generation and
discretization.
"""

### External dependencies
import mpmath
import numpy as np
import pandas as pd

### Local dependencies
from pykasso.model.geologic_features import GeologicFeature
from pykasso.core._namespaces import DEFAULT_FMM_COSTS
### Typing
from pykasso.core.grid import Grid
from numpy.random import Generator


class Fractures(GeologicFeature):
    """Class modeling the fracturation model."""
    
    def __init__(self,
                 grid: Grid,
                 rng: Generator,
                 *args,
                 **kwargs):
        """
        Class modeling the fracturation model.

        Parameters
        ----------
        grid : Grid
            pyKasso's ``Grid`` of the model.
        rng : Generator
            Random Generator Number of numpy.
        """
        # Set super constructor parameters
        label = 'fractures'
        dim = 3
        super().__init__(grid, label, dim, *args, **kwargs)
        
        # Set class attributes
        self.rng = rng
        self.i = self.stats.index.max()
        self.families = pd.DataFrame()
        self.fractures = pd.DataFrame()
        self.fractures_voxelized = {}
        
    def set_names(self,
                  names: dict[int, str],
                  default_name: str = 'family {}',
                  ) -> None:
        return super().set_names(names, default_name)

    def set_costs(self,
                  costs: dict[int, str],
                  default_cost: float = DEFAULT_FMM_COSTS['fractures'],
                  ) -> None:
        return super().set_costs(costs, default_cost)
    
    def set_model(self,
                  model: dict[int, str],
                  default_model: bool = True,
                  ) -> None:
        model.setdefault(0, False)
        return super().set_model(model, default_model)
        
    def generate_fracture_family(self,
                                 name: str,
                                 settings: dict,
                                 cost: float = DEFAULT_FMM_COSTS['fractures'],
                                 ) -> None:
        """
        Create a new fracture family.
        Populate the `self.families` and `self.fractures` dataframe attributes.

        Parameters
        ----------
        name : str
            Name of the fracture family.
        settings : dict
            Parameters for the fracture generation.
        cost : float, optional
            Travel cost of the fracture family. By default the default travel
            cost value for fractures.
        """
        # Update fracture families counter
        self.i = self.i + 1
        
        # Create a dataframe storing the new family settings
        frac_family = {
            'name': [name],
            'density': settings['density'],
            'cost': cost,
            # 'geology': [0,1,2] # TODO - add geology
        }
        frac_family = pd.DataFrame(data=frac_family, index=[self.i])
        
        # Generate fractures
        fractures = self.generate_fractures(**settings)
        fractures.insert(0, 'family_id', self.i)
        
        # Update tables
        self.fractures = pd.concat([self.fractures, fractures])
        self.families = pd.concat([self.families, frac_family])
        
        # Update dicts
        self.names |= self.families['name'].to_dict()
        self.costs |= self.families['cost'].to_dict()
        index = self.families.index.to_list()
        self.model |= {i: True for i in index}
        return None
     
    def compute_model(self) -> None:
        """
        Construct the array for the fracturation model by voxelizing the
        fractures.
        
        Voxelize each fracture family and populate the `fractures_voxelized`
        dictionary attribute. Sort fracture families by travel cost and
        combine all the voxels according to their ranking.
        """
        # Voxelize the fractures
        fractures_voxelized = {}
        for i in self.families.index:
            fractures = self.fractures[self.fractures['family_id'] == i]
            if fractures.empty:
                fractures_voxelized[i] = np.zeros_like(self.grid.data_volume)
            else:
                fractures_voxelized[i] = voxelize_fractures(self.grid,
                                                            fractures)
                
        # Sort fracture families by cost
        self.families = self.families.sort_values(['cost'])
        fractures_family_ids = self.families.index
        
        # Construct the fracturation array
        frac_model = np.zeros_like(self.grid, dtype=np.float64)
        for family_id in fractures_family_ids:
            frac_model = np.where(
                fractures_voxelized[family_id] == 1,
                family_id,
                frac_model
            )
    
        # Update the attributes
        self.data_volume = frac_model.copy()
        self.fractures_voxelized = fractures_voxelized
                
        return None
                
    def generate_fractures(self,
                           density: float,
                           orientation: float,
                           dip: float,
                           length: float,
                           orientation_distribution: str = 'vonmises',
                           dip_distribution: str = 'vonmises',
                           length_distribution: str = 'power',
                           **kwargs: dict,
                           ) -> pd.DataFrame:
        """TODO"""
        
        ######################
        ### INITIALIZATION ###
        ######################
        
        ### Set optional parameters
        alpha = kwargs.get('alpha', 2)  # TODO : Default value ??
        
        ### Set angle parameters
        orientation_min = 0
        orientation_max = 360
        dip_min = 0
        dip_max = 360
        
        ### Define if parameter is a range or an unique value
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
        Lx = self.grid.xmax - self.grid.xmin
        Ly = self.grid.ymax - self.grid.ymin
        Lz = self.grid.zmax - self.grid.zmin

        shift_x = min(Lx / 2, length_max / 2)
        shift_y = min(Ly / 2, length_max / 2)
        shift_z = min(Lz / 2, length_max / 2)

        Lex = 2 * shift_x + Lx
        Ley = 2 * shift_y + Ly
        Lez = 2 * shift_z + Lz

        area = Lex * Ley
        xmin = self.grid.xmin - shift_x
        ymin = self.grid.ymin - shift_y
        zmin = self.grid.zmin - shift_z

        ### Total numbers of fractures
        fractures_numbers = np.array(density) * area

        ### Generate poisson number for each fracture family
        real_frac_number = self.rng.poisson(fractures_numbers)
        
        ###########################
        ### Calculate fractures ###
        ###########################

        ##### FRACTURE CENTER LOCATION

        # Get fracture position from triple uniform distribution
        xm = Lex * self.rng.random(size=real_frac_number) + xmin
        ym = Ley * self.rng.random(size=real_frac_number) + ymin
        zm = Lez * self.rng.random(size=real_frac_number) + zmin

        ##### FRACTURE ORIENTATION
        
        if orientation_min > orientation_max:
            orientation_min = orientation_min - 360
        
        # No distribution case
        if orientation_distribution == 'unique':
            orientation_fractures = [orientation] * real_frac_number
        
        # Uniform distribution case
        elif orientation_distribution == 'uniform':
            orientation_fractures = self._uniform(orientation_min,
                                                  orientation_max,
                                                  real_frac_number)
        
        # Von Mises distribution case
        elif orientation_distribution == 'vonmises':
            orientation_fractures = self._vonmises(orientation_min,
                                                   orientation_max,
                                                   real_frac_number)
        
        ##### FRACTURE DIP
        if dip_min > dip_max:
            dip_min = dip_min - 360

        # No distribution case
        if dip_distribution == 'unique':
            dip_fractures = [dip] * real_frac_number
        
        # Uniform distribution case
        elif dip_distribution == 'uniform':
            dip_fractures = self._uniform(dip_min,
                                          dip_max,
                                          real_frac_number)
            
        # Von Mises distribution case
        elif dip_distribution == 'vonmises':
            dip_fractures = self._vonmises(dip_min,
                                           dip_max,
                                           real_frac_number)

        ##### FRACTURE LENGHT

        # No distribution case
        if length_distribution == 'unique':
            radius = [length / 2] * real_frac_number

        # Uniform distribution case
        elif length_distribution == 'uniform':
            radius = self._uniform(length_min,
                                   length_max,
                                   real_frac_number)

        # Trucated power law distribution case
        elif length_distribution == 'power':
            frac_length = self._power(length_min, length_max,
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
    
    def _uniform(self,
                 value_min: float,
                 value_max: float,
                 n: int,
                 ) -> np.ndarray:
        """TODO"""
        out = self.rng.uniform(low=value_min, high=value_max, size=n)
        return out
    
    @staticmethod
    def _vonmises_calculate_mu(theta_min: float,
                               theta_max: float,
                               ) -> float:
        """
        Calculate the mu parameter from the von Mises distribution.
        """
        mu = (theta_min + theta_max) / 2
        return mu
    
    @staticmethod
    def _vonmises_calculate_kappa(theta_max: float,
                                  mu: float,
                                  ) -> float:
        """
        Calculate the kappa parameter from the von Mises distribution.
        """
        std = (theta_max - mu) / 3
        kappa = Fractures._vonmises_solve_kappa(std)
        return kappa

    @staticmethod
    def _vonmises_solve_kappa(std: float):
        """
        Solve the variance equation for von Mises distribution by finding the
        corresponding kappa.
        """
        def func(kappa):
            a = std**2 - 1
            b = mpmath.besseli(1, kappa) / mpmath.besseli(0, kappa)
            return a + b
        kappa = mpmath.findroot(func, 0)
        return kappa
    
    def _vonmises(self,
                  theta_min: float,
                  theta_max: float,
                  n: int,
                  ) -> np.ndarray:
        """
        https://en.wikipedia.org/wiki/Von_Mises_distribution
        """
        theta_min_rad = np.radians(theta_min)
        theta_max_rad = np.radians(theta_max)
        mu = Fractures._vonmises_calculate_mu(theta_min_rad, theta_max_rad)
        kappa = Fractures._vonmises_calculate_kappa(theta_max_rad, mu)
        out = self.rng.vonmises(mu, kappa, size=n)
        out = np.degrees(out)
        return out
    
    def _power(self,
               value_min: float,
               value_max: float,
               alpha: float,
               n: int,
               ) -> np.ndarray:
        """
        TODO
        """
        palpha = (1 - alpha)
        invpalpha = 1 / palpha
        fmina = value_min**palpha
        frangea = value_max**palpha - fmina
        u = self.rng.random(size=n)
        out = (fmina + u * frangea)**invpalpha
        return out


def calculate_normal(dip, orientation):
    """
    TODO
    """
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
        Returns 4 np.nan if there is no intersection.
    """
    nx, ny = n[0], n[1]  # For readibility

    if not _float_eq(np.linalg.norm(n), 1):
        print("WARNING - unitcintl function : norm of n must be equal to 1")

    if np.abs(d) > 1:  # Case with no intersection
        return np.array([np.nan, np.nan, np.nan, np.nan])

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
        plane and disk if there is an intersection. Returns 4 np.nan if there
        is no intersection.
    """

    xc, yc, zc = center[0], center[1], center[2]  # For readibility

    tau = R**2 - (zi - zc)**2

    if tau <= 0:  # Avoid computing square root of negative number
        # print("The plane does not touch the sphere")
        x1, y1, x2, y2 = np.nan, np.nan, np.nan, np.nan

    else:
        tau = np.sqrt(tau)
        b = n[2] * (zc - zi) / tau / np.sqrt(1 - n[2]**2)

        if b > 1:
            # print("The plane does not touch the disk")
            x1, y1, x2, y2 = np.nan, np.nan, np.nan, np.nan

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
        plane and disk if there is an intersection. Returns 4 np.nan if there
        is no intersection.
    """

    xc, yc, zc = center[0], center[1], center[2]  # For readibility

    tau = R**2 - (xi - xc)**2

    if tau <= 0:  # Avoid computing square root of negative number
        # print("The plane does not touch the sphere")
        y1, z1, y2, z2 = np.nan, np.nan, np.nan, np.nan

    else:
        tau = np.sqrt(tau)
        b = n[0] * (xc - xi) / tau / np.sqrt(1 - n[0]**2)

        if b > 1:
            # print("The plane does not touch the disk")
            y1, z1, y2, z2 = np.nan, np.nan, np.nan, np.nan

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
        plane and disk if there is an intersection. Returns 4 np.nan if there
        is no intersection.
    """

    xc, yc, zc = center[0], center[1], center[2]  # For readibility

    tau = R**2 - (yi - yc)**2
    if tau <= 0:  # Avoid computing square root of negative number
        # print("The plane does not touch the sphere")
        x1, z1, x2, z2 = np.nan, np.nan, np.nan, np.nan

    else:
        tau = np.sqrt(tau)
        b = n[0] * (yc - yi) / tau / np.sqrt(1 - n[0]**2)

        if b > 1:
            # print("The plane does not touch the disk")
            x1, z1, x2, z2 = np.nan, np.nan, np.nan, np.nan

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


def voxelize_fractures(grid,
                       fractures: pd.DataFrame,
                       ) -> np.ndarray:
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
