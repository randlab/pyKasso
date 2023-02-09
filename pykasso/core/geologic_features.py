"""
TODO
"""

import sys
import PIL
import numpy as np
import pandas as pd

from .fracturation import generate_fractures, voxelize_fractures

# TODO
# - Gestion du format .grd (http://peterbird.name/guide/grd_format.htm)

this = sys.modules['pykasso.core.sks']

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
            elif extension in ['csv', 'txt']:
                return self._set_data_from_csv(data)
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
    
    
    def _set_data_from_csv(self, data):
        """ 
        TODO
        """
        return np.genfromtxt(data, delimiter=',')


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
        
        # Reads image
        pil_image = PIL.Image.open(data).convert('L')

        if self.dim == 3:

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
            
        if self.dim == 3:

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
        # self._compute_data_surface(domain)
        
    def _compute_statistics(self, grid):
        """
        TODO
        """
        values, counts = np.unique(self.data_volume, return_counts=True)
        stats = {
            'counts' : counts,
            'freq'   : counts / grid.nodes,
            'volume' : counts * grid.node_volume,
        }
        self.stats = pd.DataFrame(data=stats, index=values)
        return None
    
    # def _compute_data_surface(self, domain):
    #     """
    #     TODO
    #     """
    #     # TODO à terminer
    #     i,j,k = np.indices(domain.data_volume.shape)
    #     test = domain.data_volume.astype('bool')
    #     surface_indices = k.max(axis=2, initial=-1, where=test)
    #     i,j = np.indices(surface_indices.shape)
    #     self.data_surface = self.data_volume[i, j, surface_indices]
    #     return None
    
    def _set_fmm_costs(self, kwargs):
        """ 
        TODO
        """
        # TODO - Sets the fmm costs
        if 'costs' in kwargs:
            self.costs = kwargs['costs']
            for i in self.stats.index:
                if i not in self.costs:
                    # TODO
                    msg = "_set_fmm_costs - TODO"
                    raise KeyError(msg)
        else:
            # TODO msg log
            self.costs = {}
            if self.label == 'geology':
                for i in self.stats.index :
                    if i == 0:
                        continue
                    self.costs[i] = this.default_fmm_costs['aquifer']
            elif self.label == 'fractures':
                for i in self.stats.index :
                    if i == 0:
                        continue
                    self.costs[i] = this.default_fmm_costs[self.label] + ((i - 1) * (1/100))
            else:
                self.costs[1] = this.default_fmm_costs[self.label]
                

        
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
        self._set_fmm_costs(kwargs)


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
        self._set_fmm_costs(kwargs)


class Fractures(Volume):
    """
    TODO
    """
    def __init__(self, data, grid, **kwargs):
        """
        TODO
        """
        label = 'fractures'

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
        self._set_fmm_costs(kwargs)
        
        
### TODO ??

# class Karst(Volume):
#     """
#     TODO
#     """
    
#     def __init__(self, data, grid, **kwargs):
#         """
#         TODO
#         """
#         label = 'Karst'
#         dim = 3
#         super().__init__(label, dim, data, grid, **kwargs)

#         # TODO
#         # Validation de l'array



#############################################################

class ConceptualModel():
    """ 
    TODO
    """
    def __init__(self, simple_conceptual_model, simple_conceptual_model_table, conceptual_model, conceptual_model_table):
        """ 
        TODO - réfléchir aux noms.
        """
        self.data_volume = conceptual_model
        self.table       = conceptual_model_table
        # self.conceptual_model_overview       = simple_conceptual_model
        # self.conceptual_model_overview_table = simple_conceptual_model_table