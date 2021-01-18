import math
import numpy as np

""" ??? is it true ???
self.z0 = z0 # z coord of centerpoint of bottom left cell
"""

class Grid():
    """
    Create a grid for the simulation.
    The calculation points are located in the center of the cells, and the origin coordinates are located at the left-bottom of the grid.

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

    Notes
    -----
    - x0, y0 and z0 are considered to be in the center of the first node.
	"""

    def __init__(self, x0, y0, z0, nx, ny, nz, dx, dy, dz):
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
        self.X, self.Y, self.Z = np.meshgrid(self.x, self.y, self.z) # 3D array of dim (nx, ny, nz) with xyz coord of each cell's centerpoint

        # Calculating the limits of the grid
        self.xlimits = [self.x0-self.dx/2, self.x[-1]+self.dx/2]
        self.ylimits = [self.y0-self.dy/2, self.y[-1]+self.dy/2]
        self.zlimits = [self.z0-self.dz/2, self.z[-1]+self.dz/2]
        self.extent  = [self.xlimits[0], self.xlimits[1], self.ylimits[0], self.ylimits[1]] # coordinates of extent for use in plt.imshow()
    
    def __str__(self):
        return "[x0, y0, z0] : ({}, {}, {})\n[nx, ny, nz] : ({}, {}, {})\n[dx, dy, dz] : ({}, {}, {})".format(self.x0, self.y0, self.z0, self.nx, self.ny, self.nz, self.dx, self.dy, self.dz)
        
    def get_i(self, x):
        """
        Retrieve i index with x-coordinate.
        """
        return int(math.ceil((x - self.x0 - self.dx/2) / self.dx)) 
    
    def get_j(self, y):
        """
        Retrieve j index with x-coordinate.
        """
        return int(math.ceil((y - self.y0 - self.dy/2) / self.dy))

    def get_k(self, z):
        """
        Retrieve k index with x-coordinate.
        """
        return int(math.ceil((z - self.z0 - self.dz/2) / self.dz))
    
    def get_x(self, i):
        """
        Retrieve x-coordinate with i index.
        Return none if i index is outside the grid.
        """
        if (i < 0) or (i > (self.nx-1)):
            return None
        else:
            return self.x[i]
    
    def get_y(self, j):
        """
        Retrieve y-coordinate with i index.
        Return none if j index is outside the grid.
        """
        if (j < 0) or (j > (self.ny-1)):
            return None
        else:
            return self.y[j]
            
    def get_z(self, k):
        """
        Retrieve z-coordinate with i index.
        Return none if k index is outside the grid.
        """
        if (k < 0) or (k > (self.nz-1)):
            return None
        else:
            return self.z[k]
            
    def inbox(self, x, y, z):
        """
        Return 1 if x, y and z-coordinate are inside the grid.
        """
        i,j,k = self.get_i(x), self.get_j(y), self.get_k(z)
        return ( (k-self.nz)*k <= 0 ) and ( (j-self.ny)*j <= 0 ) and ( (i-self.nx)*i <= 0 )