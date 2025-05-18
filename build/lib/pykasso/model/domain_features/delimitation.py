import numpy as np
from shapely.geometry import Point, Polygon
from ...core.grid import Grid


class Delimitation():
    """
    Class modeling the vertical limits of the study site.
    """
    
    def __init__(
        self,
        vertices: list,
        grid: Grid,
    ) -> None:
        """
        Construct the delimitation, the vertical limits of the study site.

        Parameters
        ----------
        vertices : list
            List of coordinates representing the vertices of the boundary
            polygon : [[x0,y0], ..., [xn, yn]]. The list must contain at least
            3 vertices.
        grid : Grid
            Grid of the model.
        """
        self.label = 'delimitation'
        self.vertices = vertices
        
        ### Set the polygon with shapely
        path_vertices = self.vertices.copy()
        self.polygon = Polygon(path_vertices)
        
        ### Sets the mask array with a numpy-array
        row, col = np.indices((grid.nx, grid.ny))
        X, Y, Z = grid.get_meshgrids()
        pts = np.column_stack((X[row, col, 0].flatten(),
                               Y[row, col, 0].flatten()))
        msk = [self.polygon.contains(Point(x, y)) for (x, y) in pts]
        msk = np.array(msk).reshape((grid.nx, grid.ny)).astype(int)
        self.data_volume = np.repeat(msk[:, :, np.newaxis], grid.nz, axis=2)
