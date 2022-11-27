"""
TODO
"""

import PIL
import numpy as np
from .fracturation import generate_fractures, voxelize_fractures

class GeologicFeature():
    """
    Class modeling a single geologic data of the studied domain.
    """
    def __init__(self, label, name, data, grid):
        """
        Creates a geologic feature.
        This class is designed to describe a particular geologic feature.

        Parameters
        ----------
        "TODO"

        Examples
        --------
        >>>
        """
        self.label = label
        self.name = name

        if isinstance(data, np.ndarray):
            self.data = data
        else:
            # Selects the right methods to load the data
            # TODO
            extension = data.split('.')[-1]
            if extension == 'gslib':
                self.data = self._set_data_from_gslib(grid, data)
            elif extension == 'npy':
                self.data = self._set_data_from_pickle(data)
            elif extension in ['png', 'jpg']:
                self.data = self._set_data_from_image(grid, data)
            else:
                self.data = self._set_data_ones(grid)


    def _set_data_ones(self, grid):
        """
        Sets data to a matrice full of ones.
        i.e. : np.ones((nx,ny,nz))
        """
        return np.ones((grid.nx, grid.ny, grid.nz), dtype=np.int8)


    def _set_data_from_gslib(self, grid, data):
        """
        Sets data from a gslib file.

        Filling method :
        1) x-dimension - West to East
        2) y-dimension - South to North
        3) z-dimension - Bottom to Top
        """
        data = np.genfromtxt(data, skip_header=3, dtype=np.int8)
        if len(data) == grid.nx * grid.ny * grid.nz:
            data = np.reshape(data, (grid.nx, grid.ny, grid.nz), order='F') #reshape to xy grid using Fortran ordering
        else:
            data = np.reshape(data, (grid.nx, grid.ny), order='F')
            data = np.repeat(data[:, :, np.newaxis], grid.nz, axis=2)
        return data


    def _set_data_from_pickle(self, data):
        """
        TODO
        """
        return np.load(data)


    def _set_data_from_image(self, grid, data):
        """
        Sets data from an image.
        The size of the image should be the same that the size of the grid, but this is optional.
        If nz > 1, the layer is horizontally repeated.
        This method usage is not recommended, it should be used only for quick testing.
        """
        # Read image
        pil_image = PIL.Image.open(data)
        # Resize the image according to grid parameters
        pil_image = pil_image.resize((grid.ny, grid.nx))
        npy_image = np.asarray(pil_image)
        # Read colors
        npy_image = (npy_image[:,:,0] == 0)*1
        npy_image = np.repeat(npy_image[:, :, np.newaxis], grid.nz, axis=2)

        # TODO - sens des images, à checker
        # image = np.transpose(image, (1,0)) #imread return image with mxn format so we also need to transpose it
        # #image = np.flipud(image)           #we need to flip it since the data reading start from bottom
        # image = np.fliplr(image)
        
        return npy_image


    def _compute_surface(self, topography):
        """
        TODO
        """
        index    = np.sum(topography, axis=2) - 1
        row, col = np.indices((self.data.shape[0], self.data.shape[1]))
        self.surface = self.data[row, col, index]
        return None


    # def _compute_statistics(self, data):
    #     """
    #     TODO
    #     """
    #     unique, counts = np.unique(data, return_counts=True)
    #     occurrence, frequency, volume = [], [], []
    #     for nbr in counts:
    #         occurrence.append(nbr)
    #         frequency.append(100*nbr/(self.grid.nx*self.grid.ny*self.grid.nz))
    #         volume.append(nbr*self.grid.dx*self.grid.dy*self.grid.dz)
    #     return {'ID':unique, 'occurrence':occurrence, 'frequency':frequency, 'volume':volume}


###########################################################################################################
### SUB CLASSES ###
###################

class Geology(GeologicFeature):
    """
    TODO
    """

    def __init__(self, name, data, grid):
        """
        TODO
        """
        label = 'Geology'
        super().__init__(label, name, data, grid)
        # self.surface = self._compute_surface(self.topography, new_geologic_feature.data)


class Karst(GeologicFeature):
    """
    TODO
    """
    
    def __init__(self, name, data, grid):
        """
        TODO
        """
        label = 'Karst'
        super().__init__(label, name, data, grid)

        # TODO
        # Validation de l'array


class Field(GeologicFeature):
    """
    TODO
    """
    
    def __init__(self, name, data, grid):
        """
        TODO
        """
        label = 'Field'
        super().__init__(label, name, data, grid)


class Faults(GeologicFeature):
    """
    TODO
    """
    
    def __init__(self, name, data, grid):
        """
        TODO
        """
        label = 'Faults'
        super().__init__(label, name, data, grid)

        # TODO
        # Validation de l'array


class Fractures(GeologicFeature):
    """
    TODO
    """
    
    def __init__(self, name, data, grid, settings={}):
        """
        TODO
        """
        label = 'Fractures'

        if data == 'random':
            hypothesis, densities, orientation_min, orientation_max, dip_min, dip_max, alpha, length_min, length_max = settings.values()
            fractures = generate_fractures(grid, densities, alpha, orientation_min, orientation_max, dip_min, dip_max, length_min, length_max)
            # fracs_array = voxelize_fractures(grid, fractures, 'pyvista')
            fracs_array = voxelize_fractures(grid, fractures, 'python')
            data = fracs_array
            self.fractures = fractures

        super().__init__(label, name, data, grid)