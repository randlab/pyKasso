"""
Module defining a class able to read external data.
"""

### Internal dependencies
import linecache as lc

### External dependencies
import numpy as np
import pandas as pd
import rasterio

### Local dependencies
from ..core._namespaces import (
    VALID_EXTENSIONS_DATA,
    VALID_EXTENSIONS_DATAFRAME,
    VALID_EXTENSIONS_IMAGE,
)
from ..core.grid import Grid

### Typing
from typing import Union


class DataReader():
    """
    Multiple format data reader class.
    
    Supported formats:
        - gslib, vox,
        - csv,
        - txt,
        - npy,
        - jpg, png,
        - grd : http://peterbird.name/guide/grd_format.htm
        - tif / asc
            - asc : https://desktop.arcgis.com/en/arcmap/latest/manage-data/
            raster-and-images/esri-ascii-raster-format.htm
    """
    
    def __init__(self, grid=None, *args, **kwargs):
        """
        Create a data reader.
        """
        if grid is None:
            self.has_grid = False
            self.grid = None
        else:
            self.has_grid = True
            self.grid = grid
            
    def _requires_grid(self):
        """
        Check that a grid is declared in order to use the method correctly.
        """
        if self.has_grid is not True:
            msg = 'To be used, this method requires a declared grid.'
            raise ValueError(msg)
        
    @staticmethod
    def _get_extension_file(filename: str,
                            valid_extensions: list[str]
                            ) -> str:
        """
        Get the extension of a filename and check its validity.
        """
        # Retrieve file extension
        extension = filename.split('.')[-1]
        
        # Remove potential capital letters
        extension = extension.lower()
        
        # Check extension validity
        if extension not in valid_extensions:
            msg = ("File with a '.{}' extension are not supported. Supported "
                   "extensions: {}".format(extension, valid_extensions))
            raise TypeError(msg)
        
        return extension
    
    #####################
    ### Get DataFrame ###
    #####################
    
    @staticmethod
    def get_dataframe_from_file(filename: str,
                                **kwargs: dict,
                                ) -> pd.DataFrame:
        """
        Convert data from a file into a pandas dataframe.

        Parameters
        ----------
        filename : str
            The path to the file that needs to be processed. This should
            include the file name and its extension.
        """
        # Get extension file
        valid_extensions = VALID_EXTENSIONS_DATAFRAME
        extension = DataReader._get_extension_file(filename, valid_extensions)
        
        ### Get adequate reader parameters
        
        # GSLIB
        if extension == 'gslib':
            df = DataReader._get_dataframe_from_gslib(filename, **kwargs)
        
        # VOX
        elif extension == 'vox':
            df = DataReader._get_dataframe_from_vox(filename, **kwargs)
        
        return df
    
    @staticmethod
    def _get_dataframe_from_gslib(filename: str,
                                  **kwargs: dict,
                                  ) -> pd.DataFrame:
        """
        Convert data from a gslib file into a pandas dataframe.
        """
        # Retrieve the number of variables in the gslib file
        n_var = int(lc.getline(filename, 2).strip())
        
        # Retrieve the names of the columns
        names = [lc.getline(filename, i).strip() for i in range(3, n_var + 3)]
        
        # Read the file
        kwargs.setdefault('sep', ' ')
        kwargs.setdefault('skiprows', n_var + 2)
        kwargs.setdefault('names', names)
        df = pd.read_csv(filename, **kwargs)
        
        return df
    
    @staticmethod
    def _get_dataframe_from_vox(filename: str,
                                **kwargs: dict,
                                ) -> pd.DataFrame:
        """
        Convert data from a vox file into a pandas dataframe.
        """
        # Read the file
        kwargs.setdefault('sep', ' ')
        kwargs.setdefault('header', 1)
        df = pd.read_csv(filename, **kwargs)
        
        return df
    
    ######################
    ### Get Data Array ###
    ######################
            
    def get_data_from_file(
        self,
        filename: str,
        extend: bool = False,
        axis: str = 'z',
        usecol: Union[int, str] = None,
        **kwargs
    ) -> np.ndarray:
        """
        Get data from a file.

        Parameters
        ----------
        filename : str
            The path to the file that needs to be processed. This should
            include the file name and its extension.
        extend : bool, default: False
            If ``True``, a 2D dataset will be expanded in 3D in the selected ``axis``.
        axis : str, default: 'z'
            The axis in which data should be expanded if necessary.
        usecol : Union[int, str], default: None
            The rank of the column to consider.

        Returns
        -------
        np.ndarray
            Numpy array containing the file data.
        """
        # Get extension file
        valid_extensions = VALID_EXTENSIONS_DATA
        extension = DataReader._get_extension_file(filename, valid_extensions)
        
        ### Select the appropriate filling function
        
        # GSLIB
        if extension == 'gslib':
            # Get data in dataframe
            df = DataReader.get_dataframe_from_file(filename)
            # Select the right column
            if usecol is None:
                usecol = 0
            if isinstance(usecol, int):
                df = df.iloc[:, usecol]
            elif isinstance(usecol, str):
                df = df[usecol]
            # Transform dataframe into array
            data = self._get_data_from_gslib_df(df)
        
        # VOX
        elif extension == 'vox':
            # Get data in dataframe
            df = DataReader.get_dataframe_from_file(filename)
            # Select the right column
            if usecol is None:
                usecol = 3
            if isinstance(usecol, int):
                df = df.iloc[:, [0, 1, 2, usecol]]
            elif isinstance(usecol, str):
                df = df[['X', 'Y', 'Z', usecol]]
                
            # Transform dataframe into array
            data = self._get_data_from_vox_df(df)
                        
        # CSV
        elif extension == 'csv':
            kwargs.setdefault('delimiter', ',')
            data = np.genfromtxt(filename, **kwargs).T
        
        # TXT
        elif extension == 'txt':
            data = np.genfromtxt(filename, **kwargs)
            
        # NPY
        elif extension == 'npy':
            data = np.load(filename)
        
        # RASTER
        elif extension in ['tif', 'tiff', 'asc']:
            data = rasterio.open(filename).read(1)
            data = np.rot90(data, k=3)
            
        # IMAGES
        elif extension in VALID_EXTENSIONS_IMAGE:
            data = DataReader._get_data_from_image(filename)
            
        ### Control the data dimension
        if self.has_grid:
            # TODO
            pass
        
        ### According to axis, repeat data along if necessary
        if self.has_grid and extend:
            if (axis.lower() == 'x'):
                data = np.repeat(data[np.newaxis, :, :], self.grid.nx, axis=0)
            elif (axis.lower() == 'y'):
                data = np.repeat(data[:, np.newaxis, :], self.grid.ny, axis=1)
            elif (axis.lower() == 'z'):
                data = np.repeat(data[:, :, np.newaxis], self.grid.nz, axis=2)
            
        return data
    
    #######################
    ### REQUIRES A GRID ###
    #######################
    
    def _get_data_from_gslib_df(self,
                                df: pd.DataFrame,
                                ) -> np.ndarray:
        """
        Transform a pandas dataframe to a numpy array.
        """
        self._requires_grid()
        
        # Transform dataframe into array
        data = df.to_numpy()
        
        # Reshape the array with the grid shape
        if len(data) == self.grid.nodes:
            new_shape = (self.grid.nx, self.grid.ny, self.grid.nz)
            data = np.reshape(data, new_shape, order='F')
        else:
            new_shape = (self.grid.nx, self.grid.ny)
            data = np.reshape(data, new_shape, order='F')
        return data
    
    def _get_data_from_vox_df(self,
                              df: pd.DataFrame,
                              ) -> np.ndarray:
        """
        Transform a pandas dataframe to a numpy array.
        """
        self._requires_grid()
        
        # Filter values out of grid
        xyz = df.iloc[:, [0, 1, 2]].to_numpy()
        df = df.assign(is_inbox=self.grid.is_inbox(xyz))
        df = df[df['is_inbox']]

        # Retrieve valid coordinates and data
        xyz = df.iloc[:, [0, 1, 2]]
        d = df.iloc[:, 3]

        # Retrieve corresponding grid indices
        i, j, k = self.grid.get_indices(xyz)

        # Create the data array
        data = np.zeros_like(self.grid.data_volume)
        data[i, j, k] = d
        
        return data
    
    def _get_data_full_2D(self, value: float) -> np.ndarray:
        """
        Set data to a 2D-matrice full of the provided value.
        """
        self._requires_grid()
        
        dtype = np.int8
        out = np.full((self.grid.nx, self.grid.ny), value, dtype=dtype)
        return out
    
    def _get_data_full_3D(self, value: float) -> np.ndarray:
        """
        Set data to a 3D-matrice full of the provided value.
        """
        self._requires_grid()
        
        dtype = np.int8
        shape = (self.grid.nx, self.grid.ny, self.grid.nz)
        out = np.full(shape, value, dtype=dtype)
        return out
    
    ###############################
    ### DOES NOT REQUIRE A GRID ###
    ###############################
    
    @staticmethod
    def _get_data_from_image(filename: str) -> np.ndarray:
        """
        Set data from an image file. The size of the image must be the
        same as the size of the grid. If nz > 1, the layer is horizontally
        repeated.
        
        .. info::
            This method usage is not recommended, it should be used only for
            quick testing.
        """
        import PIL
        from PIL import Image
        
        # Reads image
        pil_image = PIL.Image.open(filename)
        pil_image = pil_image.convert('L')
        pil_image = pil_image.transpose(Image.FLIP_TOP_BOTTOM)
        data = np.asarray(pil_image).T
        return data
    
    #########################
    ### Only Read Methods ###
    #########################
    
    @staticmethod
    def read_vox(
        filename: str,
        usecol: Union[int, str] = None,
    ) -> np.ndarray:
        """
        Read a vox file.

        Parameters
        ----------
        filename : str
            The path to the file that needs to be processed. This should
            include the file name and its extension.
        usecol : Union[int, str], default: None
            The rank of the column to consider.

        Returns
        -------
        np.ndarray
            Numpy array of the collected data.
        """
        
        # Read file
        df = DataReader._get_dataframe_from_vox(filename)
        
        # Retrieve axis
        x = df['X'].unique()
        y = df['Y'].unique()
        z = df['Z'].unique()
        
        # Retrieve grid parameters
        nx, ny, nz = len(x), len(y), len(z)
        x0, y0, z0 = x.min(), y.min(), z.min()
        dx, dy, dz = (x[1] - x[0]), (y[1] - y[0]), (z[1] - z[0])
        
        # Create a grid
        grid = Grid(x0, y0, z0, nx, ny, nz, dx, dy, dz)
        
        # Retrieve corresponding grid indices
        xyz = df.iloc[:, [0, 1, 2]]
        if usecol is None:
            usecol = 3
        if isinstance(usecol, int):
            d = df.iloc[:, usecol]
        elif isinstance(usecol, str):
            d = df[usecol]
        i, j, k = grid.get_indices(xyz)
        
        # Construct the array
        data = np.zeros_like(grid.data_volume)
        data[i, j, k] = d
        
        return data
