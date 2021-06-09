import math
import numpy as np

class Grid():
    """
    Class modeling the structured grid of the studied domain.
    """

    def __init__(self, x0, y0, z0, nx, ny, nz, dx, dy, dz):
        """
        Create a structured grid according to the parameters.
        The calculation points are located in the center of the cells.
        The origin (x0, y0, z0) is located at the left-bottom of the grid.

        Parameters
        ----------
        x0 : float
            x-coordinate of centerpoint of bottom left cell.
        y0 : float
            y-coordinate of centerpoint of bottom left cell.
        z0 : float
            z-coordinate of centerpoint of bottom left cell.
        nx : int
            Number of cells in x direction.
        ny : int
            Number of cells in y direction.
        nz : int
            Number of cells in z direction.
        dx : float
            Cell width in x direction.
        dy : float
            Cell width in y direction.
        dz : float
            Cell width in z direction.

        Examples
        --------
        >>> grid = pk.Grid(0, 0, 0, 10, 10, 10, 1, 1, 1)
        """
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.dx = dx
        self.dy = dy
        self.dz = dz

        # Calculating the meshgrid
        self.x = np.arange(self.x0, self.x0+(self.nx)*self.dx, self.dx, dtype=np.float_) # 1D array of centerpoints of each cell along x axis
        self.y = np.arange(self.y0, self.y0+(self.ny)*self.dy, self.dy, dtype=np.float_) # 1D array of centerpoints of each cell along y axis
        self.z = np.arange(self.z0, self.z0+(self.nz)*self.dz, self.dz, dtype=np.float_) # 1D array of centerpoints of each cell along z axis
        self.X, self.Y, self.Z = np.meshgrid(self.x, self.y, self.z, indexing="ij")      # 3D array of dim (nx, ny, nz) with xyz coord of each cell's centerpoint
        # https://numpy.org/doc/stable/reference/generated/numpy.meshgrid.html

        # Calculating the limits of the grid
        self.xlimits = [self.x0-self.dx/2, self.x[-1]+self.dx/2]
        self.ylimits = [self.y0-self.dy/2, self.y[-1]+self.dy/2]
        self.zlimits = [self.z0-self.dz/2, self.z[-1]+self.dz/2]
        self.extent  = [self.xlimits[0], self.xlimits[1], self.ylimits[0], self.ylimits[1]] # coordinates of extent for plt.imshow()

        # Declaring some useful variables
        self.xmin = self.xlimits[0] # coordinates of leftmost edge of leftmost cells
        self.xmax = self.xlimits[1] # coordinates of rightmost edge of rightmost cells
        self.ymin = self.ylimits[0] # coordinates of bottom edge of bottom cells
        self.ymax = self.ylimits[1] # coordinates of top edge of top cells
        self.zmin = self.zlimits[0] # coordinates of deepest edge of deeptest cells
        self.zmax = self.zlimits[1] # coordinates of shallowest edge of shallowest cells

    def __str__(self):
        return "[x0, y0, z0] : ({}, {}, {})\n[nx, ny, nz] : ({}, {}, {})\n[dx, dy, dz] : ({}, {}, {})".format(self.x0, self.y0, self.z0, self.nx, self.ny, self.nz, self.dx, self.dy, self.dz)

    def get_i(self, x):
        """
        Retrieve i-index with x-coordinate.

        Parameters
        ----------
        x : float
            x-coordinate to test.

        Returns
        -------
        result : int
            i-index.

        Examples
        --------
        >>> i = grid.get_i(3.141)
        """
        return int(math.ceil((x - self.x0 - self.dx/2) / self.dx))

    def get_j(self, y):
        """
        Retrieve j-index with y-coordinate.

        Parameters
        ----------
        y : float
            y-coordinate to test.

        Returns
        -------
        result : int
            j-index.

        Examples
        --------
        >>> j = grid.get_j(3.141)
        """
        return int(math.ceil((y - self.y0 - self.dy/2) / self.dy))

    def get_k(self, z):
        """
        Retrieve k-index with z-coordinate.

        Parameters
        ----------
        z : float
            z-coordinate to test.

        Returns
        -------
        result : int
            k-index.

        Examples
        --------
        >>> k = grid.get_k(3.141)
        """
        return int(math.ceil((z - self.z0 - self.dz/2) / self.dz))

    def get_x(self, i):
        """
        Retrieve x-coordinate with i-index.
        Return none type if i-index is outside the grid.

        Parameters
        ----------
        i : int
            i-index to test.

        Returns
        -------
        result : float
            x-coordinate.

        Examples
        --------
        >>> x = grid.get_x(42)
        """
        if (i < 0) or (i > (self.nx-1)):
            return None
        else:
            return self.x[i]

    def get_y(self, j):
        """
        Retrieve y-coordinate with j-index.
        Return none type if j-index is outside the grid.

        Parameters
        ----------
        j : int
            j-index to test.

        Returns
        -------
        result : float
            y-coordinate.

        Examples
        --------
        >>> y = grid.get_y(42)
        """
        if (j < 0) or (j > (self.ny-1)):
            return None
        else:
            return self.y[j]

    def get_z(self, k):
        """
        Retrieve z-coordinate with k-index.
        Return none type if k-index is outside the grid.

        Parameters
        ----------
        k : int
            k-index to test.

        Returns
        -------
        result : float
            z-coordinate.

        Examples
        --------
        >>> z = grid.get_z(42)
        """
        if (k < 0) or (k > (self.nz-1)):
            return None
        else:
            return self.z[k]

    def is_inbox(self, x, y, z):
        """
        Test if a (x, y, z) point is inside the grid.

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
            1 if True.

        Examples
        --------
        >>> bool = grid.inbox(0.53, 11.4, 6.0)
        """
        i, j, k = self.get_i(x), self.get_j(y), self.get_k(z)
        return ( (k-self.nz)*k <= 0 ) and ( (j-self.ny)*j <= 0 ) and ( (i-self.nx)*i <= 0 )
