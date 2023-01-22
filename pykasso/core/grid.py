"""
Module to model the grid.
"""

import math
import numpy as np
from matplotlib.path import Path

class Grid():
    """
    Class modeling the structured grid of the studied domain.
    Warning: none of the attributes are designed to be modified.
    """

    def __init__(self, x0:float, y0:float, z0:float, nx:int, ny:int, nz:int, dx:float, dy:float, dz:float) -> None:
        """
        Creates a structured grid according to the parameters.
        The calculation points are located in the center of the cells.
        The origin (x0, y0, z0) is located at the bottom left corner of the grid.

        Parameters
        ----------
        x0 : float
            x-coordinate of centerpoint of bottom left cell.
        y0 : float
            y-coordinate of centerpoint of bottom left cell.
        z0 : float
            z-coordinate of centerpoint of bottom left cell.
        nx : int
            Number of cells in the x-direction.
        ny : int
            Number of cells in the y-direction.
        nz : int
            Number of cells in the z-direction.
        dx : float
            Cell width in the x-direction.
        dy : float
            Cell width in the y-direction.
        dz : float
            Cell width in the z-direction.

        Examples
        --------
        >>> grid = pk.Grid(0, 0, 0, 10, 10, 10, 1, 1, 1)
        """
        # TODO
        # Memory usage settings
        precision = np.float64

        # Initialization
        self.__x0 = x0
        self.__y0 = y0
        self.__z0 = z0
        self.__nx = nx
        self.__ny = ny
        self.__nz = nz
        self.__dx = dx
        self.__dy = dy
        self.__dz = dz

        # Calculating the meshgrid
        self.__x = np.arange(self.__x0, self.__x0+(self.__nx)*self.__dx, self.__dx, dtype=precision) # 1D array of centerpoints of each cell along x axis
        self.__y = np.arange(self.__y0, self.__y0+(self.__ny)*self.__dy, self.__dy, dtype=precision) # 1D array of centerpoints of each cell along y axis
        self.__z = np.arange(self.__z0, self.__z0+(self.__nz)*self.__dz, self.__dz, dtype=precision) # 1D array of centerpoints of each cell along z axis
        self.__X, self.__Y, self.__Z = np.meshgrid(self.__x, self.__y, self.__z, indexing="ij")      # 3D array of dim (nx, ny, nz) with xyz coord of each cell's centerpoint

        # Calculating the limits of the grid
        self.__xlimits = [self.__x0-self.__dx/2, self.__x[-1]+self.__dx/2]
        self.__ylimits = [self.__y0-self.__dy/2, self.__y[-1]+self.__dy/2]
        self.__zlimits = [self.__z0-self.__dz/2, self.__z[-1]+self.__dz/2]
        self.__extent  = [self.__xlimits[0], self.__xlimits[1], self.__ylimits[0], self.__ylimits[1]] # coordinates of extent for plt.imshow()

        # Declaring some useful variables
        self.__xmin = self.__xlimits[0] # coordinates of leftmost edge of leftmost cells
        self.__xmax = self.__xlimits[1] # coordinates of rightmost edge of rightmost cells
        self.__ymin = self.__ylimits[0] # coordinates of bottom edge of bottom cells
        self.__ymax = self.__ylimits[1] # coordinates of top edge of top cells
        self.__zmin = self.__zlimits[0] # coordinates of deepest edge of deeptest cells
        self.__zmax = self.__zlimits[1] # coordinates of shallowest edge of shallowest cells

        # Misc
        self.__area   = (nx * dx) * (ny * dy)
        self.__volume = self.__area * (nz * dz)
        self.__nodes = nx * ny * nz
        self.__node_volume = dx * dy * dz
        x_path = [self.__xlimits[0], self.__xlimits[0], self.__xlimits[1], self.__xlimits[1], self.__xlimits[0]]
        y_path = [self.__ylimits[0], self.__ylimits[1], self.__ylimits[1], self.__ylimits[0], self.__ylimits[0]]
        self.__path = Path(list(zip(x_path, y_path)))

    def __repr__(self):
        return "Grid({}, {}, {}, {}, {}, {}, {}, {}, {})".format(self.x0, self.y0, self.z0, self.nx, self.ny, self.nz, self.dx, self.dy, self.dz)

    def __str__(self):
        return "Grid\n[x0, y0, z0] : ({}, {}, {})\n[nx, ny, nz] : ({}, {}, {})\n[dx, dy, dz] : ({}, {}, {})".format(self.x0, self.y0, self.z0, self.nx, self.ny, self.nz, self.dx, self.dy, self.dz)

    ###############
    ### GETTERS ###
    ###############

    @property
    def x0(self):
        return self.__x0

    @property
    def y0(self):
        return self.__y0

    @property
    def z0(self):
        return self.__z0

    @property
    def nx(self):
        return self.__nx

    @property
    def ny(self):
        return self.__ny

    @property
    def nz(self):
        return self.__nz

    @property
    def dx(self):
        return self.__dx

    @property
    def dy(self):
        return self.__dy

    @property
    def dz(self):
        return self.__dz

    @property
    def x(self):
        return self.__x

    @property
    def y(self):
        return self.__y

    @property
    def z(self):
        return self.__z

    @property
    def X(self):
        return self.__X

    @property
    def Y(self):
        return self.__Y

    @property
    def Z(self):
        return self.__Z

    @property
    def xlimits(self):
        return self.__xlimits

    @property
    def ylimits(self):
        return self.__ylimits

    @property
    def zlimits(self):
        return self.__zlimits

    @property
    def extent(self):
        return self.__extent

    @property
    def xmin(self):
        return self.__xmin

    @property
    def xmax(self):
        return self.__xmax

    @property
    def ymin(self):
        return self.__ymin

    @property
    def ymax(self):
        return self.__ymax

    @property
    def zmin(self):
        return self.__zmin

    @property
    def zmax(self):
        return self.__zmax

    @property
    def area(self):
        return self.__area

    @property
    def volume(self):
        return self.__volume

    @property
    def nodes(self):
        return self.__nodes

    @property
    def node_volume(self):
        return self.__node_volume

    @property
    def path(self):
        return self.__path

    ###############
    ### METHODS ###
    ###############

    def get_i(self, x: float) -> int:
        """
        Retrieves i-index with x-coordinate.

        Parameters
        ----------
        x : float
            x-coordinate.

        Returns
        -------
        result : int
            i-index.

        Examples
        --------
        >>> i = grid.get_i(3.141)
        """
        x = np.array(x)
        out = np.ceil((x - self.x0 - self.dx/2) / self.dx)
        out = out.astype('int32')
        return out

    def get_j(self, y: float) -> int:
        """
        Retrieves j-index with y-coordinate.

        Parameters
        ----------
        y : float
            y-coordinate.

        Returns
        -------
        result : int
            j-index.

        Examples
        --------
        >>> j = grid.get_j(3.141)
        """
        y = np.array(y)
        out = np.ceil((y - self.y0 - self.dy/2) / self.dy)
        out = out.astype('int32')
        return out

    def get_k(self, z: float) -> int:
        """
        Retrieves k-index with z-coordinate.

        Parameters
        ----------
        z : float
            z-coordinate.

        Returns
        -------
        result : int
            k-index.

        Examples
        --------
        >>> k = grid.get_k(3.141)
        """
        z = np.array(z)
        out = np.ceil((z - self.z0 - self.dz/2) / self.dz)
        out = out.astype('int32')
        return out

    def get_x(self, i: int) -> float:
        """
        Retrieves x-coordinate with i-index.
        Returns none type if i-index is outside the grid.

        Parameters
        ----------
        i : int
            i-index.

        Returns
        -------
        result : float
            x-coordinate.

        Examples
        --------
        >>> x = grid.get_x(42)
        """
        i = np.array(i)
        i = i.astype('int32')
        out = np.where((i < 0) | (i > (self.nx-1)), None, self.x[i])
        return out

    def get_y(self, j: int) -> float:
        """
        Retrieves y-coordinate with j-index.
        Returns none type if j-index is outside the grid.

        Parameters
        ----------
        j : int
            j-index.

        Returns
        -------
        result : float
            y-coordinate.

        Examples
        --------
        >>> y = grid.get_y(42)
        """
        j = np.array(j)
        j = j.astype('int32')
        out = np.where((j < 0) | (j > (self.ny-1)), None, self.y[j])
        return out

    def get_z(self, k: int) -> float:
        """
        Retrieves z-coordinate with k-index.
        Returns none type if k-index is outside the grid.

        Parameters
        ----------
        k : int
            k-index.

        Returns
        -------
        result : float
            z-coordinate.

        Examples
        --------
        >>> z = grid.get_z(42)
        """
        k = np.array(k)
        k = k.astype('int32')
        out = np.where((k < 0) | (k > (self.nz-1)), None, self.z[k])
        return out

    def is_inbox(self, x: float, y: float, z: float) -> bool:
        """
        Tests if a (x, y, z) point is inside the grid.

        Parameters
        ----------
        x : float
            x-coordinate.
        y : float
            y-coordinate.
        z : float
            z-coordinate.

        Returns
        -------
        result : bool
            Boolean.

        Examples
        --------
        >>> test = grid.is_inbox(0.53, 11.4, 6.0)
        """
        i, j, k = self.get_i(x), self.get_j(y), self.get_k(z)
        return bool( (k-self.nz+1)*k <= 0 ) and ( (j-self.ny+1)*j <= 0 ) and ( (i-self.nx+1)*i <= 0 )

    def _get_property(self, property: str):
        """
        TODO
        """
        if property == 'X':
            return self.X
        elif property == 'Y':
            return self.Y
        elif property == 'Z':
            return self.Z
        else:
            print('property has not been recognized')
            return None
