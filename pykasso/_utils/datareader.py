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
from pykasso.core._namespaces import (VALID_EXTENSIONS_DATA,
                                      VALID_EXTENSIONS_DATAFRAME,
                                      VALID_EXTENSIONS_IMAGE)

### Typing
from typing import Union


class DataReader():
    """
    Multiple format data reader class.
    
    Supported formats:
        - gslib, vox, 
        - csv :
        - txt :
        - npy :
        - jpg / png
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
        TODO
        """
        if self.has_grid is not True:
            msg = 'To be used, this method requires a declared grid.'
            raise ValueError(msg)
        
    @staticmethod
    def _get_extension_file(filename: str,
                            valid_extensions: list[str]
                            ) -> str:
        """TODO"""
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
        """TODO"""
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
        """TODO"""
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
        """TODO"""
        # Read the file
        kwargs.setdefault('sep', ' ')
        kwargs.setdefault('header', 1)
        df = pd.read_csv(filename, **kwargs)
        
        return df
    
    ######################
    ### Get Data Array ###
    ######################
            
    def get_data_from_file(self,
                           filename: str,
                           extend: bool = False,
                           axis: str = 'z',
                           usecol: Union[int, str] = None,
                           **kwargs
                           ) -> np.ndarray:
        """
        TODO
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
            data = np.genfromtxt(filename, **kwargs).T
            
        # NPY
        elif extension == 'npy':
            data = np.load(data)
            
        # # RASTER # TODO ???
        # # elif extension == 'asc':
        #     # data = np.genfromtxt(fname, skip_header=6)
        # # elif extension == 'grd':
        #     # data = np.genfromtxt(fname, skip_header=2)
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
        """TODO"""
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
        """TODO"""
        self._requires_grid()

        # Filter values out of grid
        df['is_inbox'] = self.grid.is_inbox(df[['X', 'Y', 'Z']])  # TODO
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
        """Set data to a 2D-matrice full of the provided value."""
        self._requires_grid()
        
        dtype = np.int8
        out = np.full((self.grid.nx, self.grid.ny), value, dtype=dtype)
        return out
    
    def _get_data_full_3D(self, value: float) -> np.ndarray:
        """Set data to a 3D-matrice full of the provided value."""
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
        Sets data from an image file. The size of the image must be the
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
