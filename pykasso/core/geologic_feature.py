"""
TODO
"""

import PIL
import numpy as np
import pandas as pd

from .fracturation import generate_fractures, voxelize_fractures

# TODO
# - Documentation
# - Gestion du format .grd (http://peterbird.name/guide/grd_format.htm)
# - Stats with mask
# - Retrieve fmm-costs dict

class GeologicFeature():
    """
    Class modeling a single geologic data of the studied domain.
    """
    def __init__(self, label, data, grid, **kwargs):
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

        if 'axis' in kwargs:
            axis = kwargs['axis']

        # Selects the right methods to load the data
            # TODO
        if isinstance(data, np.ndarray):
            self.data = data
        else:
            extension = data.split('.')[-1]
            if extension == 'gslib':
                self.data = self._set_data_from_gslib(grid, data)
            elif extension == 'npy':
                self.data = self._set_data_from_pickle(data)
            elif extension in ['png', 'jpg']:
                self.data = self._set_data_from_image(grid, data, axis)
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
        data = np.genfromtxt(data, skip_header=3)
        if len(data) == grid.nx * grid.ny * grid.nz:
            data = np.reshape(data, (grid.nx, grid.ny, grid.nz), order='F') #reshape to xy grid using Fortran ordering
        else:
            data = np.reshape(data, (grid.nx, grid.ny), order='F')
            
            if self.label == 'Topography':
                pass
            else:
                data = np.repeat(data[:, :, np.newaxis], grid.nz, axis=2)
        return data


    def _set_data_from_grd(self):
        """
        TODO
        """
        pass


    def _set_data_from_pickle(self, data):
        """
        TODO
        """
        return np.load(data)


    def _set_data_from_image(self, grid, data, axis):
        """
        Sets data from an image.
        The size of the image should be the same that the size of the grid, but this is optional.
        If nz > 1, the layer is horizontally repeated.
        This method usage is not recommended, it should be used only for quick testing.
        """
        
        # Read image
        pil_image = PIL.Image.open(data).convert('L')

        if (axis.lower() == 'x'):
            pil_image = pil_image.resize((grid.ny, grid.nz))
        elif (axis.lower() == 'y'):
            pil_image = pil_image.resize((grid.nx, grid.nz))
        elif (axis.lower() == 'z'):
            pil_image = pil_image.resize((grid.nx, grid.ny))
        
        # npy_image = (npy_image[:,:] == 0)*1
        npy_image = np.asarray(pil_image).T + 1000
        n_colors = np.unique(npy_image)
        for i, color in enumerate(np.flip(n_colors)):
            npy_image = np.where(npy_image == color, i, npy_image)
            
        if (axis.lower() == 'x'):
            npy_image = np.repeat(npy_image[np.newaxis, :, :], grid.nx, axis=0)
        elif (axis.lower() == 'y'):
            npy_image = np.repeat(npy_image[:, np.newaxis, :], grid.ny, axis=1)
        elif (axis.lower() == 'z'):
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
        index = np.sum(topography, axis=2) - 1
        index = index.astype(np.int)
        i, j = np.indices((self.data.shape[0], self.data.shape[1]))
        self.surface = self.data[i, j, index]
        return None


    def _compute_statistics(self, grid):
        """
        TODO
        """
        values, counts = np.unique(self.data, return_counts=True)
        stats = {
            # 'id'     : values,
            'counts' : counts,
            'freq'   : counts / grid.nodes,
            'volume' : counts * grid.node_volume,
        }
        self.stats = pd.DataFrame(data=stats, index=values)
        # TODO - stats with mask
        return None


###########################################################################################################
### SUB CLASSES ###
###################

class Topography(GeologicFeature):
    """
    TODO
    """
    def __init__(self, data, grid, **kwargs):
        """
        TODO
        """
        label = 'Topography'

        super().__init__(label, data, grid, **kwargs)
        
        if not np.all((self.data == 0) | (self.data == 1)):
            self.surface_to_volume(grid)

    
    def surface_to_volume(self, grid):
        """
        TODO
        """
        k = grid.get_k(self.data)

        topography = np.zeros((grid.nx, grid.ny, grid.nz))
        for z in range(grid.nz):
            topography[:, :, z] = z
            topography[:, :, z] = np.where(topography[:, :, z] >= k, 1, 0)

        self.data = topography

        return None


    def _compute_topographic_surface(self):
        """
        TODO
        """
        self.surface_indices = np.sum(self.data, axis=2) - 1
        return None
        



class Geology(GeologicFeature):
    """
    TODO
    """

    def __init__(self, data, grid, **kwargs):
        """
        TODO
        """
        label = 'Geology'
        super().__init__(label, data, grid, **kwargs)



class Karst(GeologicFeature):
    """
    TODO
    """
    
    def __init__(self, data, grid, **kwargs):
        """
        TODO
        """
        label = 'Karst'
        super().__init__(label, data, grid, **kwargs)

        # TODO
        # Validation de l'array


class Field(GeologicFeature):
    """
    TODO
    """
    
    def __init__(self, data, grid, **kwargs):
        """
        TODO
        """
        label = 'Field'
        super().__init__(label, data, grid, **kwargs)


class Faults(GeologicFeature):
    """
    TODO
    """
    
    def __init__(self, data, grid, **kwargs):
        """
        TODO
        """
        label = 'Faults'
        super().__init__(label, data, grid, **kwargs)

        # TODO
        # Validation de l'array


class Fractures(GeologicFeature):
    """
    TODO
    """
    
    def __init__(self, data, grid, **kwargs):
        """
        TODO
        """
        label = 'Fractures'

        if 'settings' in kwargs:
            fractures_families_settings = kwargs['settings'].values()

            fractures_df = pd.DataFrame()
            fractures_npy = {}

            # Generates one array for each fractures family
            for (i, fractures_family_settings) in enumerate(fractures_families_settings):
                fractures_df_ = generate_fractures(grid=grid, rng=kwargs['rng'], **fractures_family_settings)
                fractures_df_.insert(0, 'family_id', i+1)
                fractures_df = pd.concat([fractures_df, fractures_df_])

                fractures_npy[i+1] = voxelize_fractures(grid, fractures_df_)
                
            # Sums fractures families
            total = sum([d for d in fractures_npy.values()])
            fractures_npy['t'] = total

            self.fractures = fractures_df
            self.FRACTURES = fractures_npy

            data = total

        super().__init__(label, data, grid, **kwargs)