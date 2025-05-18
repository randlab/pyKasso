"""
Module defining the structured grid.
"""

### External dependencies
import numpy as np
from shapely.geometry import Polygon

### Typing
from typing import Union


class Grid():
    """
    Class modeling a structured grid.
    
    This class models the three dimensional structured grid of the studied
    domain.
    
    Attributes
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
    """

    def __init__(
        self,
        x0: Union[int, float],
        y0: Union[int, float],
        z0: Union[int, float],
        nx: int,
        ny: int,
        nz: int,
        dx: Union[int, float],
        dy: Union[int, float],
        dz: Union[int, float],
    ) -> None:
        """
        Construct a three dimensional structured grid of nx-nodes-length,
        ny-nodes-width and nz-nodes-depth. Node centers are located in the
        center of the dx-length, dy-width and dz-depth cells. The (x0,y0,z0)
        origin is based at the bottom left corner of the grid.
        
        .. note::
            Its attributes are protected and thus not designed to be modified.

        Parameters
        ----------
        x0 : Union[int, float]
            x-coordinate of centerpoint of bottom left cell.
        y0 : Union[int, float]
            y-coordinate of centerpoint of bottom left cell.
        z0 : Union[int, float]
            z-coordinate of centerpoint of bottom left cell.
        nx : int
            Number of cells in the x-axis.
        ny : int
            Number of cells in the y-axis.
        nz : int
            Number of cells in the z-axis.
        dx : Union[int, float]
            Cell length in the x-axis.
        dy : Union[int, float]
            Cell width in the y-axis.
        dz : Union[int, float]
            Cell depth in the z-axis.
        """
        
        ### Initialization
        self.__x0 = float(x0)
        self.__y0 = float(y0)
        self.__z0 = float(z0)
        self.__nx = nx
        self.__ny = ny
        self.__nz = nz
        self.__dx = float(dx)
        self.__dy = float(dy)
        self.__dz = float(dz)

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
        msg = ("Grid(x0={}, y0={}, z0={}, nx={}, ny={}, nz={}, dx={}, dy={},"
               " dz={})".format(self.x0, self.y0, self.z0, self.nx, self.ny,
                                self.nz, self.dx, self.dy, self.dz))
        return msg
    
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
    def xlimits(self) -> tuple[float, float]:
        """Get the x-axis border limits."""
        return self.__xlimits

    @property
    def ylimits(self) -> tuple[float, float]:
        """Get the y-axis border limits."""
        return self.__ylimits

    @property
    def zlimits(self) -> tuple[float, float]:
        """Get the z-axis border limits."""
        return self.__zlimits

    @property
    def extent(self) -> tuple[float, float, float, float]:
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
        return np.ones((self.nx, self.ny, self.nz))
    
    @property
    def shape(self) -> tuple[int, int, int]:
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
        Return a dictionary containing the grid parameters.
        
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
    
    def get_meshgrids(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Return a tuple with the X, Y and Z-axis meshgrids.
        
        Returns
        -------
        out : (np.ndarray, np.ndarray, np.ndarray)
            (X, Y, Z) meshgrids
        """
        X, Y, Z = np.meshgrid(self.x, self.y, self.z, indexing="ij")
        out = (X, Y, Z)
        return out

    def get_i(self, x: Union[float, list, tuple, np.ndarray]) -> np.ndarray:
        """
        Retrieve i-index from x-coordinate.

        Parameters
        ----------
        x : array_like
            x-coordinate.

        Returns
        -------
        i : np.ndarray
            i-index.
        """
        x = np.array(x)
        i = np.floor((x - (self.x0 - (self.dx / 2))) / self.dx)
        i = i.astype('int32')
        return i

    def get_j(self, y: Union[float, list, tuple, np.ndarray]) -> np.ndarray:
        """
        Retrieve j-index from y-coordinate.

        Parameters
        ----------
        y : array_like
            y-coordinate.

        Returns
        -------
        j : np.ndarray
            j-index.
        """
        y = np.array(y)
        j = np.floor((y - (self.y0 - (self.dy / 2))) / self.dy)
        j = j.astype('int32')
        return j

    def get_k(self, z: Union[float, list, tuple, np.ndarray]) -> np.ndarray:
        """
        Retrieve k-index from z-coordinate.

        Parameters
        ----------
        z : array_like
            z-coordinate.

        Returns
        -------
        k : np.ndarray
            k-index.
        """
        z = np.array(z)
        k = np.floor((z - (self.z0 - (self.dz / 2))) / self.dz)
        k = k.astype('int32')
        return k
    
    def get_indices(
        self,
        coordinates: Union[
            list,
            tuple[float, float],
            tuple[float, float, float],
            np.ndarray
        ]
    ) -> np.ndarray:
        """
        Retrieve both i and j-index from x and y-coordinates. Returns also
        k-index when z-coordinate is provided.

        Parameters
        ----------
        coordinates : array_like
            (x, y) or (x, y, z) coordinates tuple.

        Returns
        -------
        out : (2,) or (3,) np.ndarray
            i and j-index, or i, j and k-index.
            
        Raises
        ------
        TypeError
            The shape of the ``coordinates`` parameter is invalid.
        """
        # Convert type in np.ndarray
        if not isinstance(coordinates, np.ndarray):
            coordinates = np.array(coordinates)
        
        # Shape the array in a (2,) array
        if len(coordinates.shape) == 1:
            coordinates = np.reshape(coordinates, (-1, coordinates.shape[0]))

        # Unzip coordinates
        if coordinates.shape[1] == 2:
            x = coordinates[:, 0]
            y = coordinates[:, 1]
            out = self.get_i(x), self.get_j(y)
        elif coordinates.shape[1] == 3:
            x = coordinates[:, 0]
            y = coordinates[:, 1]
            z = coordinates[:, 2]
            out = self.get_i(x), self.get_j(y), self.get_k(z)
        else:
            print(coordinates.shape)
            msg = ("Shape of the `coordinates` parameter is invalid. Only (2, "
                   "n) and (3, n) array like are valid.")
            raise TypeError(msg)
        
        # # Output a more convinient format when only one entry
        # if len(out) == 1:
        #     out = out[0]
        return out

    def get_x(self, i: Union[int, list, tuple, np.ndarray]) -> np.ndarray:
        """
        Retrieve x-coordinate from i-index.

        Parameters
        ----------
        i : array_like
            i-index.

        Returns
        -------
        out : np.ndarray
            x-coordinate.
        """
        i = np.array(i)
        i = i.astype('int32')
        out = np.where((i < 0) | (i > (self.nx - 1)), None, self.x[i])
        return out

    def get_y(self, j: Union[int, list, tuple, np.ndarray]) -> np.ndarray:
        """
        Retrieve y-coordinate from j-index.

        Parameters
        ----------
        j : array_like
            j-index.

        Returns
        -------
        out : np.ndarray
            y-coordinate.
        """
        j = np.array(j)
        j = j.astype('int32')
        out = np.where((j < 0) | (j > (self.ny - 1)), None, self.y[j])
        return out

    def get_z(self, k: Union[int, list, tuple, np.ndarray]) -> np.ndarray:
        """
        Retrieve z-coordinate from k-index.

        Parameters
        ----------
        k : array_like
            k-index.

        Returns
        -------
        out : np.ndarray
            z-coordinate.
        """
        k = np.array(k)
        k = k.astype('int32')
        out = np.where((k < 0) | (k > (self.nz - 1)), None, self.z[k])
        return out
    
    def get_coordinates(
        self,
        indices: Union[
            list,
            tuple[int, int],
            tuple[int, int, int],
            np.ndarray
        ]
    ) -> np.ndarray:
        """
        Retrieve both x and y-coordinates from i and j-index. Returns also
        z-coordinate when k-index is provided.

        Parameters
        ----------
        point : array_like
            (i, j) or (i, j, k) indices tuple.

        Returns
        -------
        out : (2,) or (3,) np.ndarray
            x and y-coordinates, or x, y and z-coordinates.
            
        Raises
        ------
        TypeError
            The shape of the ``indices`` parameter is invalid.
        """
        # Convert type in np.ndarray
        if not isinstance(indices, np.ndarray):
            indices = np.array(indices)
        
        # Shape the array in a (2,) array
        if len(indices.shape) == 1:
            indices = np.reshape(indices, (-1, indices.shape[0]))

        # Unzip coordinates
        if indices.shape[1] == 2:
            i = indices[:, 0]
            j = indices[:, 1]
            out = self.get_x(i), self.get_y(j)
        elif indices.shape[1] == 3:
            i = indices[:, 0]
            j = indices[:, 1]
            k = indices[:, 2]
            out = self.get_x(i), self.get_y(j), self.get_z(k)
        else:
            msg = ("Shape of the `indices` parameter is invalid. Only (2, "
                   "n) and (3, n) array like are valid.")
            raise TypeError(msg)
        
        return out
        
    def is_x_valid(
        self,
        x: Union[float, list, tuple, np.ndarray],
    ) -> np.ndarray:
        """
        Return true if the ``x`` coordinate is inside the grid domain.

        Parameters
        ----------
        x : array_like
            x-coordinate.

        Returns
        -------
        out : np.ndarray as boolean
        """
        i = self.get_i(x)
        out = (i - self.nx + 1) * i <= 0
        return out
    
    def is_y_valid(
        self,
        y: Union[float, list, tuple, np.ndarray],
    ) -> np.ndarray:
        """
        Return true if the ``y`` coordinate is inside the grid domain.

        Parameters
        ----------
        y : array_like
            y-coordinate.

        Returns
        -------
        out : np.ndarray as boolean
        """
        j = self.get_j(y)
        out = (j - self.ny + 1) * j <= 0
        return out
    
    def is_z_valid(
        self,
        z: Union[float, list, tuple, np.ndarray],
    ) -> np.ndarray:
        """
        Return true if the ``z`` coordinate is inside the grid domain.

        Parameters
        ----------
        z : array_like
            z-coordinate.

        Returns
        -------
        out : np.ndarray as boolean
        """
        k = self.get_k(z)
        out = (k - self.nz + 1) * k <= 0
        return out

    def is_inbox(
        self,
        coordinates: Union[
            list,
            tuple[float, float],
            tuple[float, float, float],
            np.ndarray,
        ]
    ) -> np.ndarray:
        """
        Return true if a (x, y)-point is within the surface of the grid or if
        a (x, y, z)-point is inside the grid.

        Parameters
        ----------
        point : array_like
            (x, y) or (x, y, z) coordinates tuple.

        Returns
        -------
        out : (2,) or (3,) np.ndarray as boolean
        """
        # Convert type in np.ndarray
        if not isinstance(coordinates, np.ndarray):
            coordinates = np.array(coordinates)
        
        # Shape the array in a (2,) array
        if len(coordinates.shape) == 1:
            coordinates = np.reshape(coordinates, (-1, coordinates.shape[0]))
            
        # Unzip coordinates
        if coordinates.shape[1] == 2:
            x = coordinates[:, 0]
            y = coordinates[:, 1]
            i, j = self.is_x_valid(x), self.is_y_valid(y)
            out = np.logical_and(i, j)
        elif coordinates.shape[1] == 3:
            x = coordinates[:, 0]
            y = coordinates[:, 1]
            z = coordinates[:, 2]
            i = self.is_x_valid(x)
            j = self.is_y_valid(y)
            k = self.is_z_valid(z)
            out = np.logical_and.reduce([i, j, k])
        else:
            msg = ("Shape of the `coordinates` parameter is invalid. Only (2, "
                   "n) and (3, n) array like are valid.")
            raise TypeError(msg)
        
        # Output a more convinient format when only one entry
        if len(out) == 1:
            out = out[0]
        return out
