"""
TODO
"""

import PIL
import numpy as np
import pandas as pd

from .fracturation import generate_fractures, voxelize_fractures

# TODO
# - Gestion du format .grd (http://peterbird.name/guide/grd_format.htm)


class GeologicFeature():
    """
    Class modeling a three dimensional geologic feature of the studied domain.
    """
    def __init__(self, label:str, dim:int, data, grid, **kwargs):
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
        self.dim   = dim
    
        # Selects the right attribute to fill
        if self.dim == 2:
            attribute = 'data_surface'
        elif self.dim == 3:
            attribute = 'data_volume'
        
        data = self._set_data(data, grid, **kwargs)
        setattr(self, attribute, data)


    def _set_data(self, data, grid, **kwargs):
        """
        Selects the right methods to load the data.
        """
        if 'axis' in kwargs:
            axis = kwargs['axis']
            
        if isinstance(data, np.ndarray):
            return data
        else:
            extension = data.split('.')[-1]
            if extension == 'gslib':
                return self._set_data_from_gslib(grid, data)
            elif extension == 'npy':
                return self._set_data_from_pickle(data)
            elif extension in ['png', 'jpg']:
                return self._set_data_from_image(grid, data, axis)
            elif extension == 'csv':
                return self._set_data_from_csv(grid, data)
            elif extension == 'txt':
                return self._set_data_from_txt(grid, data)
            else:
                return self._set_data_ones(grid)
                
                
    def _set_data_ones(self, grid):
        """
        Sets data to a matrice full of ones.
        i.e. : np.ones((nx,ny,nz))
        """
        if self.dim == 2:
            return np.ones((grid.nx, grid.ny), dtype=np.int8)
        elif self.dim == 3:
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
            
            if self.dim == 3:
                data = np.repeat(data[:, :, np.newaxis], grid.nz, axis=2)
        return data
    
    
    def _set_data_from_csv(self, grid, data):
        """ 
        TODO
        """
        data = np.genfromtxt(data, delimiter=',').T
        if self.dim == 3:
            data = np.repeat(data[:, :, np.newaxis], grid.nz, axis=2)
        return data
    
    
    def _set_data_from_txt(self, grid, data):
        """ 
        TODO
        """
        data = np.genfromtxt(data)
        if self.dim == 3:
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
        from PIL import Image
        
        # Reads image
        pil_image = PIL.Image.open(data).convert('L').transpose(Image.FLIP_TOP_BOTTOM)
        npy_image = np.asarray(pil_image).T
        
        if self.label in ['faults', 'fractures']:
            npy_image = (npy_image[:,:] == 0)*1
        else:        
            npy_image = npy_image + 1000
            n_colors = np.unique(npy_image)
            for i, color in enumerate(np.flip(n_colors)):
                npy_image = np.where(npy_image == color, i+1, npy_image)
            
        if self.dim == 3:

            if (axis.lower() == 'x'):
                npy_image = np.repeat(npy_image[np.newaxis, :, :], grid.nx, axis=0)
            elif (axis.lower() == 'y'):
                npy_image = np.repeat(npy_image[:, np.newaxis, :], grid.ny, axis=1)
            elif (axis.lower() == 'z'):
                npy_image = np.repeat(npy_image[:, :, np.newaxis], grid.nz, axis=2)

        return npy_image

#############################################################################

class Surface(GeologicFeature):
    """ 
    TODO
    """
    def __init__(self, label, data, grid, **kwargs):
        """
        TODO
        """
        dim = 2
        super().__init__(label, dim, data, grid, **kwargs)
    
    
    def _surface_to_volume(self, condition, grid):
        """
        TODO
        """
        k = grid.get_k(self.data_surface)
        data_volume = np.zeros((grid.nx, grid.ny, grid.nz))
        for z in range(grid.nz):
            data_volume[:, :, z] = z
            if condition == '>=':
                data_volume[:, :, z] = np.where(data_volume[:, :, z] >= k, 1, 0)
            elif condition == '=':
                data_volume[:, :, z] = np.where(data_volume[:, :, z] == k, 1, 0) 
            elif condition == '<=':
                data_volume[:, :, z] = np.where(data_volume[:, :, z] <= k, 1, 0)
        return data_volume

   
class Volume(GeologicFeature):
    """ 
    TODO
    """
    def __init__(self, label, data, grid, **kwargs):
        """
        TODO
        """
        dim = 3
        super().__init__(label, dim, data, grid, **kwargs)
        self._compute_statistics(grid)
        self.costs = kwargs['costs']
        
    def _compute_statistics(self, grid):
        values, counts = np.unique(self.data_volume, return_counts=True)
        stats = {
            'counts' : counts,
            'freq'   : counts / grid.nodes,
            'volume' : counts * grid.node_volume,
        }
        self.stats = pd.DataFrame(data=stats, index=values)
        return None
    
        
#############################################################################
### 2D Objects ###
##################

class Bedding(Surface):
    """
    TODO
    """
    def __init__(self, data, grid, **kwargs):
        label = 'bedding'
        
        super().__init__(label, data, grid, **kwargs)
        
        self._surface_to_volume('=', grid)

#############################################################################
### 3D Objects ###
##################

class Geology(Volume):
    """
    TODO
    """
    def __init__(self, data, grid, **kwargs):
        """
        TODO
        """
        label = 'geology'
        super().__init__(label, data, grid, **kwargs)


class Faults(Volume):
    """
    TODO
    """
    def __init__(self, data, grid, **kwargs):
        """
        TODO
        """
        label = 'faults'
        super().__init__(label, data, grid, **kwargs)


class Fractures(Volume):
    """
    TODO
    """
    def __init__(self, data, grid, Geology, **kwargs):
        """
        TODO
        """
        label = 'fractures'

        if 'settings' in kwargs:
            fractures_families_settings = kwargs['settings']
            
            self.table     = pd.DataFrame()
            self.fractures_families = {}

            # Generates one array for each fractures family
            for (i, fractures_family_settings) in enumerate(fractures_families_settings):
                
                settings = fractures_families_settings[fractures_family_settings]
                settings['density'] = float(settings['density']) # TODO - à mettre dans validations
                if 'cost' in settings:
                    del settings['cost']
                generated_fractures = generate_fractures(grid=grid, rng=kwargs['rng'], **settings)
                generated_fractures.insert(0, 'family_id', i+1)
                self.table = pd.concat([self.table, generated_fractures])

                self.fractures_families[i+1] = voxelize_fractures(grid, generated_fractures)
            
            # Constructs the model for fracturation
            frac_model = np.zeros_like(grid.data_volume)
            fractures_family_ids = list(self.fractures_families.keys())
            fractures_family_ids = sorted(fractures_family_ids, key=lambda family_id: kwargs['costs'][family_id], reverse=True)
            for family_id in fractures_family_ids:
                frac_model = np.where(self.fractures_families[family_id] == 1, family_id, frac_model)
            self.fractures_families['model'] = frac_model
                
            # Sums fractures families
            frac_sum = sum([d for d in self.fractures_families.values()])
            self.fractures_families['sum'] = frac_sum
            
            # Constraints model with geology if provided
            if 'geology' in kwargs:
                frac_model_geology = np.zeros_like(frac_model)
                for geologic_id in kwargs['geology']:
                    frac_model_geology = np.where(Geology.data_volume == geologic_id, frac_model, frac_model_geology)
                self.fractures_families['model'] =frac_model_geology
                    
            data = self.fractures_families['model']
                
        super().__init__(label, data, grid, **kwargs)
        
        
#############################################################

class ConceptualModel():
    """ 
    TODO
    """
    def __init__(self, conceptual_model, conceptual_model_table):
        """ 
        TODO - réfléchir aux noms.
        """
        self.data_volume = conceptual_model
        self.table       = conceptual_model_table
        # self.conceptual_model_overview       = simple_conceptual_model
        # self.conceptual_model_overview_table = simple_conceptual_model_table