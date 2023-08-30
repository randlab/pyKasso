"""
This module contains classes modeling the domain of the study site.
"""

### Local dependencies
from .geologic_features import Surface

### External dependencies
import numpy as np
from shapely.geometry import Point, Polygon
from scipy.ndimage import binary_dilation

### Typing
from pykasso._typing import Grid, Delimitation, Topography, Bedrock, WaterLevel

### TODO
# Describe each subdomain


class Domain():
    """
    Class modeling the spatial extension of the domain within the study site
    grid locating where karstic network generation is allowed.
    """
    
    ### List of available subdomains
    subdomains = [
        # domain
        "domain",
        "domain_surface",
        "domain_bottom",
        "domain_out",
        # borders
        "domain_borders",
        "domain_borders_sides",
        "domain_borders_surface",
        "domain_borders_bottom",
        # vadose
        "vadose_zone",
        "vadose_borders",
        # phreatic
        "phreatic_zone",
        "phreatic_surface",
        "phreatic_borders_surface",
        # bedrock
        "bedrock",
        "bedrock_",
        "bedrock_vadose",
        "bedrock_phreatic",
    ]
    
    def __init__(self, grid: Grid, delimitation: Delimitation = None,
                 topography: Topography = None, bedrock: Bedrock = None,
                 water_level: WaterLevel = None) -> None:
        """
        Constructs an array modeling the valid spatial extension for karst
        networks calculation.

        Parameters
        ----------
        grid : Grid
            Grid of the model.
        delimitation : Delimitation, optional
            Vertical limits of the model, by default None.
        topography : Topography, optional
            Horizontal upper limits of the model, by default None.
        bedrock : Bedrock, optional
            Horizontal lower limits of the model, by default None.
        water_level : WaterLevel, optional
            Phreatic/vadose limit of the model, by default None.
            
        Notes
        -----
        - TODO : Work in progress
        """
        ### Initialization
        self.grid = grid
        self.data_volume = np.zeros_like(grid.data_volume)
        self.delimitation = delimitation
        self.topography = topography
        self.bedrock = bedrock
        self.water_level = water_level
        self._is_defined = {}
        
        ### Tests if values are defined
        features = ['delimitation', 'topography', 'bedrock', 'water_level']
        for feature in features:
            if eval(feature) is None:
                self._is_defined[feature] = False
            else:
                self._is_defined[feature] = True
  
        ### Computes domain extension
        # Computes volumetric domain extension
        if delimitation is not None:
            self.data_volume += delimitation.data_volume
        if topography is not None:
            self.data_volume += topography.data_volume
        if bedrock is not None:
            self.data_volume += np.logical_not(bedrock.data_volume).astype(int)
        # only keeps cells where all the classes join
        test = self.data_volume == self.data_volume.max()
        self.data_volume = np.where(test, 1, 0)
            
        # Computes apparent surface from domain according to axis
        self.data_surfaces = {
            'x': self.data_volume.max(axis=0),
            'y': self.data_volume.max(axis=1),
            'z': self.data_volume.max(axis=2),
        }
        self.surface = self.data_surfaces['x'].sum() * self.grid.node_area
        
    def get_subdomain(self, subdomain: str) -> np.ndarray:
        """Returns the array modeling the requested subdomain.

        Parameters
        ----------
        subdomain : str
            Name of the requested subdomain. Subdomains available: {}
            
        Returns
        -------
        out : np.ndarray
        """
        ### DOMAIN ###
        if subdomain == 'domain':
            out = self.data_volume
        elif subdomain == 'domain_surface':
            out = self._get_surface_subdomain('up')
        elif subdomain == 'domain_bottom':
            out = self._get_surface_subdomain('down')
        ### BORDERS ###
        elif subdomain == 'domain_borders':
            out = self._get_bordered_subdomain(subdomain)
        elif subdomain == 'domain_borders_sides':
            out = self._get_bordered_subdomain(subdomain)
        elif subdomain == 'domain_borders_surface':
            # 'Domain borders' x 'Face up'
            borders = self._get_bordered_subdomain('domain_borders_sides')
            surf_up = self._get_surface_subdomain('up')
            out = np.logical_and(borders, surf_up)
        elif subdomain == 'domain_borders_bottom':
            # 'Domain borders' x 'Face down'
            borders = self._get_bordered_subdomain('domain_borders_sides')
            surf_dw = self._get_surface_subdomain('down')
            out = np.logical_and(borders, surf_dw)
        ### VADOSE ###
        elif subdomain == 'vadose_zone':
            out = self._get_phreatic_subdomain('vadose_zone')
        elif subdomain == 'vadose_borders':
            # 'Domain borders' x 'Vadose zone'
            borders = self._get_bordered_subdomain('domain_borders_sides')
            vadose_zone = self._get_phreatic_subdomain('vadose_zone')
            out = np.logical_and(borders, vadose_zone)
        ### PHREATIC ###
        elif subdomain == 'phreatic_zone':
            out = self._get_phreatic_subdomain('phreatic_zone')
        elif subdomain == 'phreatic_borders':
            # 'Domain borders' x 'Phreatic zone'
            borders = self._get_bordered_subdomain('domain_borders_sides')
            phreatic_zone = self._get_phreatic_subdomain('phreatic_zone')
            out = np.logical_and(borders, phreatic_zone)
        elif subdomain == 'phreatic_surface':
            out = self._get_phreatic_subdomain('phreatic_surface')
        elif subdomain == 'phreatic_surface_borders':
            # 'Domain borders' x 'Phreatic surface'
            borders = self._get_bordered_subdomain('domain_borders_sides')
            phreatic_surface = self._get_phreatic_subdomain('phreatic_surface')
            out = np.logical_and(borders, phreatic_surface)
        ### BEDROCK ###
        elif subdomain == 'bedrock':
            out = self.bedrock.data_volume
        elif subdomain == 'bedrock_':
            out = self._get_bedrock_subdomain()
        elif subdomain == 'bedrock_vadose':
            # 'Bedrock' x 'Vadose zone'
            bedrock = self._get_bedrock_subdomain()
            vadose_zone = self._get_phreatic_subdomain('vadose_zone')
            out = np.logical_and(bedrock, vadose_zone)
        elif subdomain == 'bedrock_phreatic':
            # 'Bedrock' x 'Phreatic zone'
            bedrock = self._get_bedrock_subdomain()
            phreatic_zone = self._get_phreatic_subdomain('phreatic_zone')
            out = np.logical_and(bedrock, phreatic_zone)
        else:
            # TODO - error
            print('TODO - ERROR domain name')
        return out
    
    def _get_surface_subdomain(self, face_name: str) -> np.ndarray:
        """Computes the selected surfaces from the volumetric domain."""
        
        # retrieves grid indices
        domain = self.data_volume.astype('bool')
        I, J, K = np.indices(domain.shape)
        
        # retrieves z-surface grid indices
        if face_name in ['up', 'down']:
            i_, j_ = np.indices(self.data_surfaces['z'].shape)
            range_ = list(range(self.grid.nz))
            if face_name == 'up':
                face = K.max(axis=2, initial=-1, where=domain)
            elif face_name == 'down':
                face = K.min(axis=2, initial=self.grid.nz, where=domain)
        
        # flattens the indices
        i_ = i_.flatten()
        j_ = j_.flatten()
        k = face.flatten()
        
        # retrieves the valid i,j,k indices
        nodes = list(zip(i_, j_, k))
        valid_nodes = [(i, j, k) for (i, j, k) in nodes if k in range_]
        i, j, k = zip(*valid_nodes)
        
        # colors the array
        volume = np.zeros_like(self.data_volume)
        volume[i, j, k] = 1
        
        return volume
    
    def _get_bordered_subdomain(self, subdomain: str) -> np.ndarray:
        """Computes the required border subdomain."""
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.binary_dilation.html
        
        # Computes the full borders extent
        if subdomain == 'domain_borders':
            k = np.zeros((3, 3, 3), dtype=int)
            k[:, 1, 1], k[1, :, 1], k[1, 1, :] = 1, 1, 1
            dilation = binary_dilation(self.data_volume == 0, k,
                                       border_value=1)
            volume = dilation & self.data_volume
            
        # Computes the x and y-axis borders
        elif subdomain == 'domain_borders_sides':
            k = np.zeros((3, 3), dtype=int)
            k[:, 1], k[1, :] = 1, 1
            volume = np.zeros_like(self.data_volume)
            for z in range(self.grid.nz):
                dilation = binary_dilation(self.data_volume[:, :, z] == 0,
                                           k,
                                           border_value=1)
                test = dilation & self.data_volume[:, :, z]
                volume[:, :, z] = test
        
        return volume
        
    def _get_phreatic_subdomain(self, subdomain: str) -> np.ndarray:
        
        # returns empty array if water level is not defined
        if self.water_level is None:
            volume = np.zeros_like(self.grid.data_volume)
            return volume
        
        # defines bedrock
        if self.bedrock is not None:
            bedrock = self.bedrock.data_volume
        else:
            bedrock = np.zeros_like(self.grid.data_volume)
            
        # defines water level
        water_level = self.water_level.data_volume
        
        # computes subdomain
        if subdomain == 'vadose_zone':
            elems = (self.data_volume, np.logical_not(bedrock),
                     np.logical_not(water_level))
            volume = np.logical_and.reduce(elems)
        elif subdomain == 'phreatic_zone':
            elems = (self.data_volume, np.logical_not(bedrock), water_level)
            volume = np.logical_and.reduce(elems)
        elif subdomain == 'phreatic_surface':
            water_surface = self.water_level._surface_to_volume('=', self.grid)
            volume = np.logical_and(self.data_volume, water_surface)
        return volume
    
    def _get_bedrock_subdomain(self) -> np.ndarray:
        
        # returns empty array if water level is not defined
        if self.bedrock is None:
            volume = np.zeros_like(self.grid.data_volume)
            return volume
            
        # defines bedrock
        roll_value = 2
        volume = np.roll(self.bedrock.data_volume, roll_value, axis=2)
        volume[:, :, 0:2] = 1
        
        return volume
  
    ###############
    ### METHODS ###
    ###############

    def is_3D_point_valid(self, point: tuple) -> bool:
        """
        Returns true if a (x, y, z)-point is inside the domain, otherwise
        false.

        Parameters
        ----------
        point : tuple
            (x, y, z)-point

        Returns
        -------
        out : bool
        """
        if self.grid.is_inbox(point):
            i, j, k = self.grid.get_indices(point)
            out = bool(self.data_volume[i, j, k])
            return out
        else:
            out = False
            return out
     
    # TODO ******************///////////////////******************************
    # def is_coordinate_in_subdomain(self, subdomain: str, x: float, y: float,
    #                                z: float) -> bool:
    #     """
    #     Returns true if a (x, y, z)-coordinate point is inside the subdomain,
    #     otherwise false.

    #     Parameters
    #     ----------
    #     subdomain : str
    #         Name of the subdomain to test: {}
    #     x : float
    #         x-coordinate.
    #     y : float
    #         y-coordinate.
    #     z : float
    #         z-coordinate.

    #     Returns
    #     -------
    #     out : bool
    #     """
    #     subdomain = self.get_subdomain(subdomain)
    #     if self.grid.is_inbox(x, y, z):
    #         i, j, k = self.grid.get_indices(x, y, z)
    #         out = bool(subdomain[i, j, k])
    #         return out
    #     else:
    #         out = False
    #         return out
    # TODO ******************///////////////////******************************
        
    def is_2D_point_valid(self, point: tuple) -> bool:
        """
        Returns true if a z-coordinate exists for a (x, y)-point projected
        inside the domain, otherwise false.

        Parameters
        ----------
        point : tuple
            (x, y)-point

        Returns
        -------
        out : bool
        """
        point = Point(point)
        if self.grid.polygon.contains(point):
            i, j = self.grid.get_indices(point)
            out = bool(self.data_surfaces['z'][i, j])
        elif self.grid.polygon.touches(point):
            x, y = point
            i, j = self.grid.get_indices(x, y)
            test_x = self.grid.is_x_valid(x)
            test_y = self.grid.is_y_valid(y)
            if test_x and test_y:
                out = bool(self.data_surfaces['z'][i, j])
            else:
                out = False
        else:
            out = False
        return out
        
    
#################################
### Domain's embedded classes ###
#################################
    
class Delimitation():
    """
    Class modeling the delimitation, the vertical limits of the study site.
    """
    
    def __init__(self, vertices: list, grid: Grid) -> None:
        """
        Constructs the delimitation, the vertical limits of the study site.

        Parameters
        ----------
        vertices : list
            List of vertices.
        grid : Grid
            Grid of the model.
        """
        self.label = 'delimitation'
        self.vertices = vertices
        
        ### Sets the polygon with shapely
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
        

class Topography(Surface):
    """
    Class modeling the topography, the horizontal upper limit of the study
    site.
    """
    
    def __init__(self, *args, **kwargs) -> None:
        label = 'topography'
        super().__init__(label, *args, **kwargs)
        
        
class Bedrock(Surface):
    """
    Class modeling the bedrock elevation, the horizontal lower limit of the
    study site.
    """
    
    def __init__(self, *args, **kwargs) -> None:
        label = 'bedrock_elevation'
        super().__init__(label, *args, **kwargs)
        
        
class WaterLevel(Surface):
    """
    Class modeling the water level elevation, the phreatic/vadose limit of the
    study site.
    """
    
    def __init__(self, *args, **kwargs) -> None:
        label = 'water_level'
        super().__init__(label, *args, **kwargs)


#####################
### Documentation ###
#####################

subdomains_list = ""
for subdomain in Domain.subdomains:
    txt_elem = """\n\t\t\t- "{}" """.format(subdomain)
    subdomains_list += txt_elem

# Updates documentation
Domain.get_subdomain.__doc__ = (Domain.get_subdomain.__doc__
                                .format(subdomains_list))
# Domain.is_coordinate_in_subdomain.__doc__ = (
#     Domain.is_coordinate_in_subdomain.__doc__.format(subdomains_list)
# )
