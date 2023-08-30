"""
Module defining the structured grid.
"""

### External dependencies
import numpy as np
from shapely.geometry import Polygon

### TODO
# Correct documentation


class Grid():
    """Defining of the pyKasso grid.
    
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

        # Calculate the 1D arrays of centerpoints of each cell along each axis
        # 3D meshgrids are available through 'X', 'Y' and 'Z' properties.
        self.__x = np.arange(self.__x0, self.__x0 + (self.__nx) * self.__dx,
                             self.__dx, dtype=np.float64)
        self.__y = np.arange(self.__y0, self.__y0 + (self.__ny) * self.__dy,
                             self.__dy, dtype=np.float64)
        self.__z = np.arange(self.__z0, self.__z0 + (self.__nz) * self.__dz,
                             self.__dz, dtype=np.float64)

        # Calculate the dimensions limits
        self.__xlimits = (self.__x0 - self.__dx / 2,
                          self.__x[-1] + self.__dx / 2)
        self.__ylimits = (self.__y0 - self.__dy / 2,
                          self.__y[-1] + self.__dy / 2)
        self.__zlimits = (self.__z0 - self.__dz / 2,
                          self.__z[-1] + self.__dz / 2)
        # Coordinates of extent for plt.imshow()
        self.__extent = (self.__xlimits[0], self.__xlimits[1],
                         self.__ylimits[0], self.__ylimits[1])

        # Declare minimum and maximum attributes for each dimension
        self.__xmin = self.__xlimits[0]
        self.__xmax = self.__xlimits[1]
        self.__ymin = self.__ylimits[0]
        self.__ymax = self.__ylimits[1]
        self.__zmin = self.__zlimits[0]
        self.__zmax = self.__zlimits[1]

        # Declare some useful attributes
        self.__area = (nx * dx) * (ny * dy)  # total area of the xy-surface
        self.__volume = self.__area * (nz * dz)  # total volume of the grid
        self.__nodes = nx * ny * nz  # total number of nodes in the grid
        self.__node_area = dx * dy  # volume of one single node
        self.__node_volume = dx * dy * dz  # volume of one single node
        self.__shape = (nx, ny, nz)  # shape of the modeling array
        
        # Other attributes used in further geometrical operations
        x_path = [self.__xlimits[0], self.__xlimits[0], self.__xlimits[1],
                  self.__xlimits[1], self.__xlimits[0]]
        y_path = [self.__ylimits[0], self.__ylimits[1], self.__ylimits[1],
                  self.__ylimits[0], self.__ylimits[0]]
        self.__surface_coordinates = tuple(zip(x_path[:-1], y_path[:-1]))
        self.__polygon = Polygon(self.__surface_coordinates)
        
    def __repr__(self) -> str:
        return str(self.get_grid_parameters())
    
    def __str__(self) -> str:
        txt = ("Grid"
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
        """Get the x-coordinate origin of the grid."""
        return self.__x0

    @property
    def y0(self) -> float:
        """Get the y-coordinate origin of the grid."""
        return self.__y0

    @property
    def z0(self) -> float:
        """Get the z-coordinate origin of the grid."""
        return self.__z0

    @property
    def nx(self) -> int:
        """Get the number of cells in the x-axis."""
        return self.__nx

    @property
    def ny(self) -> int:
        """Get the number of cells in the y-axis."""
        return self.__ny

    @property
    def nz(self) -> int:
        """Get the number of cells in the z-axis."""
        return self.__nz

    @property
    def dx(self) -> float:
        """Get the cell length in the x-axis."""
        return self.__dx

    @property
    def dy(self) -> float:
        """Get the cell width in the y-axis."""
        return self.__dy

    @property
    def dz(self) -> float:
        """Get the cell depth in the z-axis."""
        return self.__dz

    @property
    def x(self) -> np.ndarray:
        """
        Get the one dimensional array of centerpoints of each cell along
        x-axis.
        """
        return self.__x

    @property
    def y(self) -> np.ndarray:
        """
        Get the one dimensional array of centerpoints of each cell along
        y-axis.
        """
        return self.__y

    @property
    def z(self) -> np.ndarray:
        """
        Get the one dimensional array of centerpoints of each cell along
        z-axis.
        """
        return self.__z

    @property
    def xlimits(self) -> list:
        """Get the x-axis border limits."""
        return self.__xlimits

    @property
    def ylimits(self) -> list:
        """Get the y-axis border limits."""
        return self.__ylimits

    @property
    def zlimits(self) -> list:
        """Get the z-axis border limits."""
        return self.__zlimits

    @property
    def extent(self) -> list:
        """Get the [xmin, xmax, ymin, ymax] limits of the model."""
        return self.__extent

    @property
    def xmin(self) -> float:
        """Get the lower limit of x-axis."""
        return self.__xmin

    @property
    def xmax(self) -> float:
        """Get the upper limit of x-axis."""
        return self.__xmax

    @property
    def ymin(self) -> float:
        """Get the lower limit of y-axis."""
        return self.__ymin

    @property
    def ymax(self) -> float:
        """Get the upper limit of y-axis."""
        return self.__ymax

    @property
    def zmin(self) -> float:
        """Get the lower limit of z-axis."""
        return self.__zmin

    @property
    def zmax(self) -> float:
        """Get the upper limit of z-axis."""
        return self.__zmax

    @property
    def area(self) -> float:
        """Get the total area of the xy-surface."""
        return self.__area
        
    @property
    def volume(self) -> float:
        """Get the total volume of the grid."""
        return self.__volume

    @property
    def nodes(self) -> int:
        """Get the total number of nodes in the grid."""
        return self.__nodes
    
    @property
    def node_area(self) -> float:
        """Get the area of one single node."""
        return self.__node_area

    @property
    def node_volume(self) -> float:
        """Get the volume of one single node."""
        return self.__node_volume
    
    @property
    def data_volume(self) -> np.ndarray:
        """Get the array modeling the grid."""
        data_volume = np.ones((self.nx, self.ny, self.nz))
        return data_volume
    
    @property
    def shape(self) -> tuple:
        """Get the dimensions of the grid."""
        return self.__shape
    
    @property
    def surface_coordinates(self) -> list:
        """Get the surface coordinates."""
        return self.__surface_coordinates

    @property
    def polygon(self) -> Polygon:
        """Get the surface polygon."""
        return self.__polygon
    
    ###############
    ### METHODS ###
    ###############
    
    def get_grid_parameters(self) -> dict:
        """
        Get the grid parameters as a dict.
        
        Returns
        -------
        grid_parameters : dict
        """
        grid_parameters = {
            'x0': self.x0,
            'y0': self.y0,
            'z0': self.z0,
            'nx': self.nx,
            'ny': self.ny,
            'nz': self.nz,
            'dx': self.dx,
            'dy': self.dy,
            'dz': self.dz,
        }
        return grid_parameters
    
    def get_meshgrids(self) -> tuple:
        """
        Get the x, y and z-axis meshgrids.
        
        Returns
        -------
        out : tuple
            X-, Y- and Z-axis meshgrids.
        """
        X, Y, Z = np.meshgrid(self.x, self.y, self.z, indexing="ij")
        out = (X, Y, Z)
        return out

    def get_i(self, x: float) -> int:
        """
        Retrieve i-index from x-coordinate.

        Parameters
        ----------
        x : float
            x-coordinate.

        Returns
        -------
        i : int
            i-index.

        Examples
        --------
        >>> grid = Grid(0, 0, 0, 10, 10, 10, 10, 20, 30)
        >>> grid.get_i(70)
        7
        """
        x = np.array(x)
        i = np.floor((x - (self.x0 - (self.dx / 2))) / self.dx)
        i = i.astype('int32')
        return i

    def get_j(self, y: float) -> int:
        """
        Retrieve j-index from y-coordinate.

        Parameters
        ----------
        y : float
            y-coordinate.

        Returns
        -------
        j : int
            j-index.

        Examples
        --------
        >>> grid = Grid(0, 0, 0, 10, 10, 10, 10, 20, 30)
        >>> grid.get_j(70)
        3
        """
        y = np.array(y)
        j = np.floor((y - (self.y0 - (self.dy / 2))) / self.dy)
        j = j.astype('int32')
        return j

    def get_k(self, z: float) -> int:
        """
        Retrieve k-index from z-coordinate.

        Parameters
        ----------
        z : float
            z-coordinate.

        Returns
        -------
        k : int
            k-index.

        Examples
        --------
        >>> grid = Grid(0, 0, 0, 10, 10, 10, 10, 20, 30)
        >>> grid.get_k(70)
        2
        """
        z = np.array(z)
        k = np.floor((z - (self.z0 - (self.dz / 2))) / self.dz)
        k = k.astype('int32')
        return k
    
    def get_indices(self, point: tuple) -> tuple:
        """
        Retrieve both i and j-index from x and y-coordinates. Returns also
        k-index when z-coordinate is provided.

        Parameters
        ----------
        point : tuple
            (x, y) or (x, y, z) coordinates tuple.

        Returns
        -------
        out : tuple
            i and j-index, or i, j and k-index.
            
        Examples
        --------
        >>> grid = Grid(0, 0, 0, 10, 10, 10, 10, 20, 30)
        >>> grid.get_indices((70, 70))
        (7, 3)
        >>> grid.get_indices((70, 70, 70))
        (7, 3, 2)
        """
        # Perform some preliminary tests
        if not isinstance(point, np.ndarray):
            point = np.array(point)
        
        if len(point.shape) == 1:
            point = np.reshape(point, (-1, point.shape[0]))
        
        # TODO - to remove ???
        # print(point)
        # print(point.shape)
        # print(len(point.shape))

        if point.shape[1] == 2:
            x, y = list(zip(*point))
            out = (self.get_i(x), self.get_j(y))
        elif point.shape[1] == 3:
            x, y, z = list(zip(*point))
            out = (self.get_i(x), self.get_j(y), self.get_k(z))
        else:
            pass
        return out

    def get_x(self, i: int) -> float:
        """
        Retrieve x-coordinate from i-index.

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
        >>> grid = Grid(0, 0, 0, 10, 10, 10, 10, 20, 30)
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
        Retrieve y-coordinate from j-index.

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
        >>> grid = Grid(0, 0, 0, 10, 10, 10, 10, 20, 30)
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
        Retrieve z-coordinate from k-index.

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
        >>> grid = Grid(0, 0, 0, 10, 10, 10, 10, 20, 30)
        >>> grid.get_z(5)
        150
        """
        k = np.array(k)
        k = k.astype('int32')
        out = np.where((k < 0) | (k > (self.nz - 1)), None, self.z[k])
        out = float(out)
        return out
    
    def get_coordinates(self, indices: tuple) -> tuple:
        """
        Retrieve both x and y-coordinates from i and j-index. Returns also
        z-coordinate when k-index is provided.

        Parameters
        ----------
        point : tuple
            (i, j) or (i, j, k) indices tuple.

        Returns
        -------
        out : tuple
            x and y-coordinates, or x, y and z-coordinates.
            
        Examples
        --------
        >>> grid = Grid(0, 0, 0, 10, 10, 10, 10, 20, 30)
        >>> grid.get_coordinates((5, 5))
        (50.0, 100.0)
        >>> grid.get_coordinates((5, 5, 5))
        (50.0, 100.0, 150.0)
        """
        if len(indices) == 2:
            i, j = indices
            out = (self.get_x(i), self.get_y(j))
        elif len(indices) == 3:
            i, j, k = indices
            out = (self.get_x(i), self.get_y(j), self.get_z(k))
        else:
            pass
        return out
        
    def is_x_valid(self, x: float) -> bool:
        """
        Return true if the x-coordinate is inside the grid domain.

        Parameters
        ----------
        x : float
            x-coordinate.

        Returns
        -------
        out : bool
        """
        i = self.get_i(x)
        out = bool((i - self.nx + 1) * i <= 0)
        return out
    
    def is_y_valid(self, y: float) -> bool:
        """
        Return true if the y-coordinate is inside the grid domain.

        Parameters
        ----------
        y : float
            y-coordinate.

        Returns
        -------
        out : bool
        """
        j = self.get_j(y)
        out = bool((j - self.ny + 1) * j <= 0)
        return out
    
    def is_z_valid(self, z: float) -> bool:
        """
        Return true if the z-coordinate is inside the grid domain.

        Parameters
        ----------
        z : float
            z-coordinate.

        Returns
        -------
        out : bool
        """
        k = self.get_k(z)
        out = bool((k - self.nz + 1) * k <= 0)
        return out

    def is_inbox(self, point: tuple) -> bool:
        """
        Return true if a (x, y)-point is within the surface of the grid or if
        a (x, y, z)-point is inside the grid.

        Parameters
        ----------
        point : tuple
            (x, y) or (x, y, z) coordinates tuple.

        Returns
        -------
        out : bool

        Examples
        --------
        >>> grid = Grid(0, 0, 0, 10, 10, 10, 10, 20, 30)
        >>> grid.is_inbox((70, 70, 70))
        True
        >>> grid.is_inbox((70, 70, 1000))
        False
        """
        if len(point) == 2:
            x, y = point
            out = self.is_x_valid(x) and self.is_y_valid(y)
        elif len(point) == 3:
            x, y, z = point
            out = (self.is_x_valid(x) and self.is_y_valid(y)
                   and self.is_z_valid(z))
        else:
            pass
        return out
