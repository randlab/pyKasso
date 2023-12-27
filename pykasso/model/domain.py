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
from pykasso.core.grid import Grid

### TODO
# - Short description of the subdomains
# -

    
class Delimitation():
    """
    Class modeling the vertical limits of the study site.
    """
    
    def __init__(self,
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
        

class Topography(Surface):
    """
    Class modeling the upper horizontal limit of the study site.
    """
    
    def __init__(self,
                 grid: Grid,
                 *args,
                 **kwargs,
                 ) -> None:
        feature = 'topography'
        super().__init__(grid, feature, *args, **kwargs)
        
        
class Bedrock(Surface):
    """
    Class modeling the lower horizontal limit of the study site.
    """
    
    def __init__(self,
                 grid: Grid,
                 *args,
                 **kwargs,
                 ) -> None:
        feature = 'bedrock'
        super().__init__(grid, feature, *args, **kwargs)
        
        
class WaterTable(Surface):
    """
    Class modeling the water level elevation, the phreatic/vadose limit of the
    study site.
    """
    
    def __init__(self,
                 grid: Grid,
                 *args,
                 **kwargs,
                 ) -> None:
        feature = 'water_table'
        super().__init__(grid, feature, *args, **kwargs)


class Domain():
    """
    Class modeling the spatial extension of the domain within the study site
    grid and locating where karstic network generation is allowed.
    """
    
    ### List of available subdomains
    subdomains = [
        "domain",                       #
        "domain_surface",               #
        "domain_bottom",                #
        "domain_out",                   #

        "domain_borders",               #
        "domain_borders_sides",         #
        "domain_borders_surface",       #
        "domain_borders_bottom",        #

        "vadose_zone",                  #
        "vadose_borders",               #

        "phreatic_zone",                #
        "phreatic_surface",             #
        "phreatic_borders_surface",     #

        "bedrock",                      #
        "bedrock_",                     #
        "bedrock_vadose",               #
        "bedrock_phreatic",             #
    ]
    
    def __init__(self,
                 grid: Grid,
                 delimitation: Delimitation = None,
                 topography: Topography = None,
                 bedrock: Bedrock = None,
                 water_table: WaterTable = None,
                 geology: np.ndarray = None,
                 ) -> None:
        """
        Construct an array modeling the valid spatial extension for karst
        conduit network generation.

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
        water_table : WaterTable, optional
            Phreatic/vadose limit of the model, by default None.
        geology : np.ndarray, optional
            Geologic model, by default None.
        """
        ### Initialization
        self.grid = grid
        self.data_volume = np.zeros_like(grid.data_volume)
        self.delimitation = delimitation
        self.topography = topography
        self.bedrock = bedrock
        self.water_table = water_table
        self._is_defined = {}
        
        ### Test if values are defined
        features = ['delimitation', 'topography', 'bedrock', 'water_table']
        for feature in features:
            if eval(feature) is None:
                self._is_defined[feature] = False
            else:
                self._is_defined[feature] = True
  
        ### Compute domain extension
        # Compute volumetric domain extension
        if delimitation is not None:
            self.data_volume += delimitation.data_volume
        
        # When topography and bedrock are both declared, control if data are
        # well inside the grid
        if (topography is not None) and (bedrock is not None):
            topography_is_in_grid = (topography.data_surface >= self.grid.zmin)
            bedrock_is_in_grid = (bedrock.data_surface >= self.grid.zmin)
            valid_domain = np.logical_and(topography_is_in_grid,
                                          bedrock_is_in_grid)
            valid_domain = np.repeat(valid_domain[:, :, np.newaxis],
                                     self.grid.nz,
                                     axis=2)
            self.data_volume += valid_domain
            
        # Add topography
        if topography is not None:
            self.data_volume += topography.data_volume
            
        # Add reversed bedrock
        if bedrock is not None:
            bedrock_vol = np.logical_not(bedrock.data_volume).astype(int)
            self.data_volume += bedrock_vol
            
        # Only keep cells where all the classes join
        test = self.data_volume == self.data_volume.max()
        self.data_volume = np.where(test, 1, 0)
            
        # Compute apparent surface from domain according to axis
        self.data_surfaces = {
            'x': self.data_volume.max(axis=0),
            'y': self.data_volume.max(axis=1),
            'z': self.data_volume.max(axis=2),
        }
        self.surface = self.data_surfaces['x'].sum() * self.grid.node_area
        
        # Union with declared geological model
        if geology is not None:
            test = np.logical_and(geology > 0, self.data_volume)
            test = test.astype(int)
            self.data_volume = np.where(test, 1, 0)
        
    def get_subdomain(self, subdomain: str) -> np.ndarray:
        """
        Return the array modeling the requested subdomain.

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
        elif subdomain == 'phreatic_borders_surface':
            # 'Domain borders' x 'Phreatic surface'
            borders = self._get_bordered_subdomain('domain_borders_sides')
            phreatic_surface = self._get_phreatic_subdomain('phreatic_surface')
            out = np.logical_and(borders.astype(bool),
                                 phreatic_surface.astype(bool))
            out = out.astype(int)
        ### BEDROCK ###
        elif subdomain == 'bedrock':
            if self.bedrock is None:
                out = np.zeros_like(self.grid.data_volume)
            else:
                out = self.bedrock.data_volume
        elif subdomain == 'bedrock_':
            if self.bedrock is None:
                out = np.zeros_like(self.grid.data_volume)
            else:
                out = self._get_bedrock_subdomain()
        elif subdomain == 'bedrock_vadose':
            if self.bedrock is None:
                out = np.zeros_like(self.grid.data_volume)
            else:
                # 'Bedrock' x 'Vadose zone'
                bedrock_r = np.invert(self.bedrock.data_volume.astype(bool))
                bedrock_vadose = self._get_bedrock_subdomain()
                bedrock = np.logical_and(bedrock_r, bedrock_vadose)
                if self._is_defined['water_table']:
                    vadose_zone = self._get_phreatic_subdomain('vadose_zone')
                else:
                    vadose_zone = np.ones_like(bedrock)
                out = np.logical_and(bedrock, vadose_zone)
        elif subdomain == 'bedrock_phreatic':
            # 'Bedrock' x 'Phreatic zone'
            bedrock = self._get_bedrock_subdomain()
            phreatic_zone = self._get_phreatic_subdomain('phreatic_zone')
            out = np.logical_and(bedrock, phreatic_zone)
        else:
            print("ERROR: subdomain '{}' does not exist.".format(subdomain))
            out = None
        return out
    
    def _get_surface_subdomain(self, face_name: str) -> np.ndarray:
        """
        Compute the selected surfaces from the volumetric domain.
        """
        
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
        test = np.isin(k, range_)
        i = i_[test]
        j = j_[test]
        k = k[test]
        
        # nodes = list(zip(i_, j_, k))
        # valid_nodes = [(i, j, k) for (i, j, k) in nodes if k in range_]
        # i, j, k = zip(*valid_nodes)
        
        # colors the array
        volume = np.zeros_like(self.data_volume)
        volume[i, j, k] = 1
        
        return volume
    
    def _get_bordered_subdomain(self, subdomain: str) -> np.ndarray:
        """
        Compute the required border subdomain.
        """
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
        """
        TODO
        """
        # return empty array if water level is not defined
        if self.water_table is None:
            volume = np.zeros_like(self.grid.data_volume)
            return volume
        
        # get bedrock
        if self.bedrock is not None:
            bedrock = self.bedrock.data_volume
        else:
            bedrock = np.zeros_like(self.grid.data_volume)
            
        # get water table
        water_table = self.water_table.data_volume
        
        # compute subdomain
        if subdomain == 'vadose_zone':
            elems = (self.data_volume.astype(bool),
                     np.logical_not(bedrock.astype(bool)),
                     np.logical_not(water_table.astype(bool)))
        elif subdomain == 'phreatic_zone':
            elems = (self.data_volume.astype(bool),
                     np.logical_not(bedrock.astype(bool)),
                     water_table.astype(bool))
        elif subdomain == 'phreatic_surface':
            water_surface = self.water_table._surface_to_volume('=', self.grid)
            elems = (self.data_volume.astype(bool),
                     water_surface.astype(bool))
            
        volume = np.logical_and.reduce(elems)
        volume = volume.astype(int)
        return volume
    
    def _get_bedrock_subdomain(self) -> np.ndarray:
        """
        TODO
        """
        # define bedrock
        roll_value = 2
        volume = np.roll(self.bedrock.data_volume, roll_value, axis=2)
        volume[:, :, 0:2] = 1
        
        return volume
  
    ###############
    ### METHODS ###
    ###############

    def is_3D_point_valid(self, point: tuple) -> bool:
        """
        Return true if a (x, y, z)-point is inside the domain, otherwise
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
        Return true if a z-coordinate exists for a (x, y)-point projected
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
    

#####################
### Documentation ###
#####################

subdomains_list = ""
for subdomain in Domain.subdomains:
    txt_elem = """\n\t\t\t- "{}" """.format(subdomain)
    subdomains_list += txt_elem

# Update documentation
Domain.get_subdomain.__doc__ = (Domain.get_subdomain.__doc__
                                .format(subdomains_list))
# Domain.is_coordinate_in_subdomain.__doc__ = (
#     Domain.is_coordinate_in_subdomain.__doc__.format(subdomains_list)
# )
