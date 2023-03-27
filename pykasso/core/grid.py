"""
This module contains a class modeling a structured grid.
"""

### External dependencies
import numpy as np
from matplotlib.path import Path

### Typing
from pykasso._typing import Path as pyplotPath


class Grid():
    """
    This class models the three dimensional structured grid of the studied
    domain.
    """

    def __init__(self, x0: float, y0: float, z0: float, nx: int, ny: int,
                 nz: int, dx: float, dy: float, dz: float) -> None:
        """
        Constructs a three dimensional structured grid of nx-nodes-length,
        ny-nodes-width and nz-nodes-depth. Node centers are located in the
        center of the dx-length, dy-width and dz-depth cells. The (x0,y0,z0)
        origin is based at the bottom left corner of the grid.
        
        .. note::
            Its attributes are protected and thus not designed to be modified.

        Parameters
        ----------
        x0 : float
            x-coordinate of centerpoint of bottom left cell.
        y0 : float
            y-coordinate of centerpoint of bottom left cell.
        z0 : float
            z-coordinate of centerpoint of bottom left cell.
        nx : int
            Number of cells in the x-axis.
        ny : int
            Number of cells in the y-axis.
        nz : int
            Number of cells in the z-axis.
        dx : float
            Cell length in the x-axis.
        dy : float
            Cell width in the y-axis.
        dz : float
            Cell depth in the z-axis.
            
        Notes
        -----
            .

        Examples
        --------
        >>> import pykasso as pk
        >>> grid = pk.Grid(0, 0, 0, 10, 10, 10, 10, 20, 30)
        """
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

        # Calculates the 1D arrays of centerpoints of each cell along each axis
        # Then calculates 3D arrays of meshgrids
        self.__x = np.arange(self.__x0, self.__x0 + (self.__nx) * self.__dx,
                             self.__dx, dtype=np.float64)
        self.__y = np.arange(self.__y0, self.__y0 + (self.__ny) * self.__dy,
                             self.__dy, dtype=np.float64)
        self.__z = np.arange(self.__z0, self.__z0 + (self.__nz) * self.__dz,
                             self.__dz, dtype=np.float64)
        self.__X, self.__Y, self.__Z = np.meshgrid(self.__x, self.__y,
                                                   self.__z, indexing="ij")

        # Calculates the dimensions limits
        self.__xlimits = [self.__x0 - self.__dx / 2,
                          self.__x[-1] + self.__dx / 2]
        self.__ylimits = [self.__y0 - self.__dy / 2,
                          self.__y[-1] + self.__dy / 2]
        self.__zlimits = [self.__z0 - self.__dz / 2,
                          self.__z[-1] + self.__dz / 2]
        # Coordinates of extent for plt.imshow()
        self.__extent = [self.__xlimits[0], self.__xlimits[1],
                         self.__ylimits[0], self.__ylimits[1]]

        # Declares minimum and maximum attributes for each dimension
        self.__xmin = self.__xlimits[0]
        self.__xmax = self.__xlimits[1]
        self.__ymin = self.__ylimits[0]
        self.__ymax = self.__ylimits[1]
        self.__zmin = self.__zlimits[0]
        self.__zmax = self.__zlimits[1]

        # Declares some useful attributes
        self.__area = (nx * dx) * (ny * dy)  # total area of the xy-surface
        self.__volume = self.__area * (nz * dz)  # total volume of the grid
        self.__nodes = nx * ny * nz  # total number of nodes in the grid
        self.__node_volume = dx * dy * dz  # volume of one single node
        self.__data_volume = np.ones((nx, ny, nz))  # array modeling the grid
        self.__shape = (nx, ny, nz)  # shape of the modeling array
        
        # Other attributes used in further geometrical operations
        x_path = [self.__xlimits[0], self.__xlimits[0], self.__xlimits[1],
                  self.__xlimits[1], self.__xlimits[0]]
        y_path = [self.__ylimits[0], self.__ylimits[1], self.__ylimits[1],
                  self.__ylimits[0], self.__ylimits[0]]
        self.__path = Path(list(zip(x_path, y_path)))
        self.__surface_coordinates = list(zip(x_path[:-1], y_path[:-1]))
        
    def __str__(self) -> str:
        txt = ("pyKasso's grid"
               "\n[x0, y0, z0] : ({}, {}, {})"
               "\n[nx, ny, nz] : ({}, {}, {})"
               "\n[dx, dy, dz] : ({}, {}, {})"
               .format(self.x0, self.y0, self.z0,
                       self.nx, self.ny, self.nz,
                       self.dx, self.dy, self.dz))
        return txt

    ###############
    ### GETTERS ###
    ###############

    @property
    def x0(self) -> float:
        """Gets the x-coordinate origin of the grid."""
        return self.__x0

    @property
    def y0(self) -> float:
        """Gets the y-coordinate origin of the grid."""
        return self.__y0

    @property
    def z0(self) -> float:
        """Gets the z-coordinate origin of the grid."""
        return self.__z0

    @property
    def nx(self) -> int:
        """Gets the number of cells in the x-axis."""
        return self.__nx

    @property
    def ny(self) -> int:
        """Gets the number of cells in the y-axis."""
        return self.__ny

    @property
    def nz(self) -> int:
        """Gets the number of cells in the z-axis."""
        return self.__nz

    @property
    def dx(self) -> float:
        """Gets the cell length in the x-axis."""
        return self.__dx

    @property
    def dy(self) -> float:
        """Gets the cell width in the y-axis."""
        return self.__dy

    @property
    def dz(self) -> float:
        """Gets the cell depth in the z-axis."""
        return self.__dz

    @property
    def x(self) -> np.ndarray:
        """
        Gets the one dimensional array of centerpoints of each cell along
        x-axis.
        """
        return self.__x

    @property
    def y(self) -> np.ndarray:
        """
        Gets the one dimensional array of centerpoints of each cell along
        y-axis.
        """
        return self.__y

    @property
    def z(self) -> np.ndarray:
        """
        Gets the one dimensional array of centerpoints of each cell along
        z-axis.
        """
        return self.__z

    @property
    def X(self) -> np.ndarray:
        """Gets the x-axis meshgrid."""
        return self.__X

    @property
    def Y(self) -> np.ndarray:
        """Gets the y-axis meshgrid."""
        return self.__Y

    @property
    def Z(self) -> np.ndarray:
        """Gets the z-axis meshgrid."""
        return self.__Z

    @property
    def xlimits(self) -> list:
        """Gets the x-axis border limits."""
        return self.__xlimits

    @property
    def ylimits(self) -> list:
        """Gets the y-axis border limits."""
        return self.__ylimits

    @property
    def zlimits(self) -> list:
        """Gets the z-axis border limits."""
        return self.__zlimits

    @property
    def extent(self) -> list:
        """Gets the [xmin, xmax, ymin, ymax] limits of the model."""
        return self.__extent

    @property
    def xmin(self) -> float:
        """Gets the lower limit of x-axis."""
        return self.__xmin

    @property
    def xmax(self) -> float:
        """Gets the upper limit of x-axis."""
        return self.__xmax

    @property
    def ymin(self) -> float:
        """Gets the lower limit of y-axis."""
        return self.__ymin

    @property
    def ymax(self) -> float:
        """Gets the upper limit of y-axis."""
        return self.__ymax

    @property
    def zmin(self) -> float:
        """Gets the lower limit of z-axis."""
        return self.__zmin

    @property
    def zmax(self) -> float:
        """Gets the upper limit of z-axis."""
        return self.__zmax

    @property
    def area(self) -> float:
        """Gets the total area of the xy-surface."""
        return self.__area
        
    @property
    def volume(self) -> float:
        """Gets the total volume of the grid."""
        return self.__volume

    @property
    def nodes(self) -> int:
        """Gets the total number of nodes in the grid."""
        return self.__nodes

    @property
    def node_volume(self) -> float:
        """Gets the volume of one single node."""
        return self.__node_volume
    
    @property
    def data_volume(self) -> np.ndarray:
        """Gets the array modeling the grid."""
        return self.__data_volume
    
    @property
    def shape(self) -> tuple:
        """Gets the dimensions of the grid."""
        return self.__shape

    @property
    def path(self) -> Path:
        """Gets the path."""
        return self.__path
    
    @property
    def surface_coordinates(self) -> list:
        """Gets the surface coordinates."""
        return self.__surface_coordinates

    ###############
    ### METHODS ###
    ###############

    def get_i(self, x: float) -> int:
        """
        Retrieves i-index from x-coordinate.

        Parameters
        ----------
        x : float
            x-coordinate.

        Returns
        -------
        out : int
            i-index.

        Examples
        --------
        >>> grid = pk.Grid(0, 0, 0, 10, 10, 10, 10, 20, 30)
        >>> grid.get_i(70)
        7
        """
        x = np.array(x)
        out = np.ceil((x - self.x0 - self.dx / 2) / self.dx)
        out = out.astype('int32')
        return out

    def get_j(self, y: float) -> int:
        """
        Retrieves j-index from y-coordinate.

        Parameters
        ----------
        y : float
            y-coordinate.

        Returns
        -------
        out : int
            j-index.

        Examples
        --------
        >>> grid = pk.Grid(0, 0, 0, 10, 10, 10, 10, 20, 30)
        >>> grid.get_j(70)
        3
        """
        y = np.array(y)
        out = np.ceil((y - self.y0 - self.dy / 2) / self.dy)
        out = out.astype('int32')
        return out

    def get_k(self, z: float) -> int:
        """
        Retrieves k-index from z-coordinate.

        Parameters
        ----------
        z : float
            z-coordinate.

        Returns
        -------
        out : int
            k-index.

        Examples
        --------
        >>> grid = pk.Grid(0, 0, 0, 10, 10, 10, 10, 20, 30)
        >>> grid.get_k(70)
        2
        """
        z = np.array(z)
        out = np.ceil((z - self.z0 - self.dz / 2) / self.dz)
        out = out.astype('int32')
        return out
    
    def get_indices(self, x: float, y: float, z: float = None) -> tuple:
        """
        Retrieves both i and j-index from x and y-coordinates. Returns also
        k-index when z-coordinate is provided.

        Parameters
        ----------
        x : float
            x-coordinate.
        y : float
            y-coordinate.
        z : float, optional
            z-coordinate, by default None.

        Returns
        -------
        out : tuple
            i and j-index, or i, j and k-index.
            
        Examples
        --------
        >>> grid = pk.Grid(0, 0, 0, 10, 10, 10, 10, 20, 30)
        >>> grid.get_indices(70, 70)
        (7, 3)
        >>> grid.get_indices(70, 70, 70)
        (7, 3, 2)
        """
        if z is None:
            out = (self.get_i(x), self.get_j(y))
            return out
        else:
            out = (self.get_i(x), self.get_j(y), self.get_k(z))
            return out

    def get_x(self, i: int) -> float:
        """
        Retrieves x-coordinate from i-index.

        Parameters
        ----------
        i : int
            i-index.

        Returns
        -------
        out : float
            x-coordinate.

        Examples
        --------
        >>> grid = pk.Grid(0, 0, 0, 10, 10, 10, 10, 20, 30)
        >>> grid.get_x(5)
        50
        """
        i = np.array(i)
        i = i.astype('int32')
        out = np.where((i < 0) | (i > (self.nx - 1)), None, self.x[i])
        out = float(out)
        return out

    def get_y(self, j: int) -> float:
        """
        Retrieves y-coordinate from j-index.

        Parameters
        ----------
        j : int
            j-index.

        Returns
        -------
        out : float
            y-coordinate.

        Examples
        --------
        >>> grid = pk.Grid(0, 0, 0, 10, 10, 10, 10, 20, 30)
        >>> grid.get_y(5)
        100
        """
        j = np.array(j)
        j = j.astype('int32')
        out = np.where((j < 0) | (j > (self.ny - 1)), None, self.y[j])
        out = float(out)
        return out

    def get_z(self, k: int) -> float:
        """
        Retrieves z-coordinate from k-index.

        Parameters
        ----------
        k : int
            k-index.

        Returns
        -------
        out : float
            z-coordinate.

        Examples
        --------
        >>> grid = pk.Grid(0, 0, 0, 10, 10, 10, 10, 20, 30)
        >>> grid.get_z(5)
        150
        """
        k = np.array(k)
        k = k.astype('int32')
        out = np.where((k < 0) | (k > (self.nz - 1)), None, self.z[k])
        out = float(out)
        return out
    
    def get_coordinates(self, i: int, j: int, k: int = None) -> tuple:
        """
        Retrieves both x and y-coordinates from i and j-index. Returns also
        z-coordinate when k-index is provided.

        Parameters
        ----------
        i : int
            i-index.
        j : int
            j-index.
        k : int
            k-index.

        Returns
        -------
        out : tuple
            x and y-coordinates, or x, y and z-coordinates.
            
        Examples
        --------
        >>> grid = pk.Grid(0, 0, 0, 10, 10, 10, 10, 20, 30)
        >>> grid.get_coordinates(5, 5)
        (50.0, 100.0)
        >>> grid.get_coordinates(5, 5, 5)
        (50.0, 100.0, 150.0)
        """
        if k is None:
            out = (self.get_x(i), self.get_y(j))
            return out
        else:
            out = (self.get_x(i), self.get_y(j), self.get_z(k))
            return out

    def is_inbox(self, x: float, y: float, z: float) -> bool:
        """
        Returns true if a (x, y, z)-coordinate point is inside the grid,
        otherwise false.

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
        out : bool

        Examples
        --------
        >>> grid = pk.Grid(0, 0, 0, 10, 10, 10, 10, 20, 30)
        >>> grid.is_inbox(70, 70, 70)
        True
        >>> grid.is_inbox(70, 70, 1000)
        False
        """
        i, j, k = self.get_i(x), self.get_j(y), self.get_k(z)
        out = (bool((k - self.nz + 1) * k <= 0)
               and ((j - self.ny + 1) * j <= 0)
               and ((i - self.nx + 1) * i <= 0))
        return out
