"""
TODO
"""

### Internal dependencies
import os
import pickle

### External dependencies
import yaml
import numpy as np
import pandas as pd

### Local dependencies
from pykasso.core.grid import Grid


class DataReader():
    """
    Multiple format data reader class.
    
    Supported formats:
        - gslib :
        - csv :
        - txt :
        - npy :
        - jpg :
        - png :
        - asc : https://desktop.arcgis.com/en/arcmap/latest/manage-data/raster-
        and-images/esri-ascii-raster-format.htm
        - grd : http://peterbird.name/guide/grd_format.htm
    """
    
    def __init__(self, *args, **kwargs):
        """Creates a data reader."""
        
    def _set_data_full_2D(self, grid, value: float) -> np.ndarray:
        """Sets data to a 2D-matrice full of the provided value."""
        dtype = np.int8
        out = np.full((grid.nx, grid.ny), value, dtype=dtype)
        return out
    
    def _set_data_full_3D(self, grid, value: float) -> np.ndarray:
        """Sets data to a 3D-matrice full of the provided value."""
        dtype = np.int8
        out = np.full((grid.nx, grid.ny, grid.nz), value, dtype=dtype)
        return out
    
    def _set_data_from_file(self, grid, location: str, extend: bool = False,
                            axis: str = 'z') -> np.ndarray:
        # Gets the extension
        extension = location.split('.')[-1]
        
        # Selects the appropriate filling function
        if extension == 'gslib':
            data = self._set_data_from_gslib(grid, location)
        elif extension == 'csv':
            data = np.genfromtxt(location, delimiter=',').T
        elif extension == 'txt':
            data = np.genfromtxt(location).T
        elif extension == 'npy':
            data = np.load(data)
        elif extension in ['jpg', 'png']:
            data = self._set_data_from_image(location)
        elif extension == 'asc':
            data = np.genfromtxt(location, skip_header=6)
        elif extension == 'grd':
            data = np.genfromtxt(location, skip_header=2)
        else:
            # TODO
            print('TODO : error, extension file not recognized')
        
        # According to axis, repeats data along if necessary
        if extend:
            if (axis.lower() == 'x'):
                data = np.repeat(data[np.newaxis, :, :], grid.nx, axis=0)
            elif (axis.lower() == 'y'):
                data = np.repeat(data[:, np.newaxis, :], grid.ny, axis=1)
            elif (axis.lower() == 'z'):
                data = np.repeat(data[:, :, np.newaxis], grid.nz, axis=2)
            
        return data
    
    def _set_data_from_gslib(self, grid, location: str) -> np.ndarray:
        """
        Sets data from a gslib file.

        Filling method :
        1) x-axis from West to East
        2) y-axis from South to North
        3) z-axis from Bottom to Top
        """
        data = np.genfromtxt(location, skip_header=3)
        if len(data) == grid.nodes:
            data = np.reshape(data, (grid.nx, grid.ny, grid.nz), order='F')
        else:
            data = np.reshape(data, (grid.nx, grid.ny), order='F')
        return data
    
    def _set_data_from_image(self, location: str) -> np.ndarray:
        """
        Sets data from a .jpg or .png file. The size of the image must be the
        same as the size of the grid. If nz > 1, the layer is horizontally
        repeated.
        
        .. info::
            This method usage is not recommended, it should be used only for
            quick testing.
        """
        import PIL
        from PIL import Image
        
        # Reads image
        pil_image = PIL.Image.open(location)
        pil_image = pil_image.convert('L')
        pil_image = pil_image.transpose(Image.FLIP_TOP_BOTTOM)
        data = np.asarray(pil_image).T
        return data
    

class ProjectReader():
    """
    TODO
    """
    
    def __init__(self, project_directory: str, *args, **kwargs) -> None:
        # Reads the project_state.yaml file
        self.project_directory = project_directory.strip("/")
        self.project_state = None
        self._update_project_state()
        x0, y0, z0 = self._get_grid_origin()
        nx, ny, nz = self._get_grid_dimensions()
        dx, dy, dz = self._get_grid_spacing()
        self.grid = Grid(x0, y0, z0, nx, ny, nz, dx, dy, dz)
        
    def _update_project_state(self) -> None:
        path = self.project_directory + "/outputs/project_state.yaml"
        with open(path, "r") as f:
            self.project_state = yaml.safe_load(f)
    
    def _get_grid_origin(self) -> tuple:
        x0 = self.project_state['grid']['x0']
        y0 = self.project_state['grid']['y0']
        z0 = self.project_state['grid']['z0']
        return (x0, y0, z0)
    
    def _get_grid_dimensions(self) -> tuple:
        nx = self.project_state['grid']['nx']
        ny = self.project_state['grid']['ny']
        nz = self.project_state['grid']['nz']
        return (nx, ny, nz)
    
    def _get_grid_spacing(self) -> tuple:
        dx = self.project_state['grid']['dx']
        dy = self.project_state['grid']['dy']
        dz = self.project_state['grid']['dz']
        return (dx, dy, dz)

    def _read_pickle(self, path: str) -> dict:
        with open(path, 'rb') as handle:
            out = pickle.load(handle)
            return out
        
    def _get_simulation_data(self, n: int) -> dict:
        simulation_directory = self.project_state['simulation_locations'][n]
        simulation_data_path = simulation_directory + 'results.pickle'
        simulation_data = self._read_pickle(simulation_data_path)
        return simulation_data


class GSLIB_Reader():
    """
    TODO
    """
    
    def __init__(self, path):
        if not os.path.exists(path):
            msg = "The path is not valid. '{}' does not exist.".format(path)
            raise FileNotFoundError(msg)
        else:
            self.location = path

        with open(path) as f:
            lines = f.readlines()
            
        self.comment = lines[0].strip()
        
        try:
            self.vars_nbr = int(lines[1])
        except Exception as err:
            raise err
            
        self.vars_names = []
        for i in range(self.vars_nbr):
            self.vars_names.append(lines[i + 2].strip())
        
        try:
            skip = 2 + self.vars_nbr
            data = np.genfromtxt(path, skip_header=skip)
            self.data = pd.DataFrame()
            if self.vars_nbr == 1:
                self.data[self.vars_names[0]] = data
            else:
                for i, name in zip(range(self.vars_nbr), self.vars_names):
                    self.data[name] = data[:, i]
        except Exception as err:
            raise err
        
    def __str__(self):
        txt = "GSLIB File Reader \n"
        txt += "File location         : {}\n".format(self.location)
        txt += "Comment line          : {}\n".format(self.comment)
        txt += "Number of variable(s) : {}\n".format(self.vars_nbr)
        txt += "Names of variable(s)  : {}\n".format(self.vars_names)
        txt += "Length of data        : {}".format(len(self.data))
        return txt
     
    def export_column(self, column, dim):
        
        data = self.data.loc[:, column].values
        
        if len(dim) == 2:
            nx, ny = dim
            data = np.reshape(data, (nx, ny), order='F')
        
        elif len(dim) == 3:
            nx, ny, nz = dim
            data = np.reshape(data, (nx, ny, nz), order='F')
            
        return data
