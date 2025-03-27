"""
This module contains classes modeling the domain of the study site.
"""

### Local dependencies
from .bedrock import Bedrock
from .delimitation import Delimitation
from .topography import Topography
from .watertable import WaterTable
from ...core.grid import Grid

### External dependencies
import numpy as np
from shapely.geometry import Point, Polygon
from scipy.ndimage import binary_dilation


class Domain():
    """
    Class modeling the spatial extension of the domain within the study site
    grid and locating where karstic network generation is allowed.
    """
    
    ### List of available subdomains
    subdomains = {
        "domain":
            "Cells located in the domain.",
        "domain_surface":
            ("Cells located at the upper interface between the domain and the"
             " outside zone."),
        "domain_bottom":
            ("Cells located at the lower interface between the domain and the"
             " outside zone."),
        # "domain_out": "",

        "domain_borders":
            ("Cells located at the interface between the domain and the"
             " outside zone."),
        "domain_borders_sides":
            ("Cells located at the interface between the domain and the"
             " outside zone, but only in the x- and y-direction."),
        "domain_borders_surface":
            ("Cells located at the upper interface between the domain and the"
             " outside zone, but only in the x- and y-direction."),
        "domain_borders_bottom":
            ("Cells located at the lower interface between the domain and the"
             " outside zone, but only in the x- and y-direction."),

        "vadose_zone":
            ("Cells located in the vadose zone."),
        "vadose_borders":
            ("Cells located in the vadose zone, at the interface between the"
             " vadose zone and the rest of the model, but only in the x- and"
             " y-direction."),

        "phreatic_zone":
            ("Cells located in the phreatic zone."),
        "phreatic_surface":
            ("Cells located in the phreatic zone, at the interface between the"
             " phreatic zone and the vadose zone."),
        "phreatic_borders_surface":
            ("Cells located in the phreatic zone, at the interface between the"
             " phreatic zone and the vadose zone, as well as the outside zone."
             ),

        "bedrock":
            ("Cells located in the zone below the bedrock surface"),
        "bedrock_":
            ("Cells located in the zone defined by the bedrock surface, with"
             " an additional two-cells layer in the upper z-direction."),
        "bedrock_vadose":
            ("Two-cells layer located above the zone defined by the bedrock"
             " surface and intersected by the vadose zone."),
        "bedrock_phreatic":
            ("Two-cells layer located above the zone defined by the bedrock"
             " surface and intersected by the phreatic zone."),
    }
    
    def __init__(
        self,
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
        delimitation : Delimitation, default: None
            Vertical limits of the model.
        topography : Topography, default: None
            Horizontal upper limits of the model.
        bedrock : Bedrock, default: None
            Horizontal lower limits of the model.
        water_table : WaterTable, default: None
            Phreatic/vadose limit of the model.
        geology : np.ndarray, default: None
            Geologic model.
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
        self.surface = self.data_surfaces['x'].sum() * self.grid.node_area # TODO - to control ! rename ?!
        
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
            Name of the requested subdomain. To see the list of available
            subdomains inspect the ``subdomains`` attribute dictionary.
            
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
            out = np.logical_and(borders.astype(bool), surf_up.astype(bool))
            out = out.astype(int)
        elif subdomain == 'domain_borders_bottom':
            # 'Domain borders' x 'Face down'
            borders = self._get_bordered_subdomain('domain_borders_sides')
            surf_dw = self._get_surface_subdomain('down')
            out = np.logical_and(borders.astype(bool), surf_dw.astype(bool))
            out = out.astype(int)
        ### VADOSE ###
        elif subdomain == 'vadose_zone':
            out = self._get_phreatic_subdomain('vadose_zone')
        elif subdomain == 'vadose_borders':
            # 'Domain borders' x 'Vadose zone'
            borders = self._get_bordered_subdomain('domain_borders_sides')
            vadose_zone = self._get_phreatic_subdomain('vadose_zone')
            out = np.logical_and(borders.astype(bool),
                                 vadose_zone.astype(bool))
            out = out.astype(int)
        ### PHREATIC ###
        elif subdomain == 'phreatic_zone':
            out = self._get_phreatic_subdomain('phreatic_zone')
        elif subdomain == 'phreatic_borders':
            # 'Domain borders' x 'Phreatic zone'
            borders = self._get_bordered_subdomain('domain_borders_sides')
            phreatic_zone = self._get_phreatic_subdomain('phreatic_zone')
            out = np.logical_and(borders.astype(bool),
                                 phreatic_zone.astype(bool))
            out = out.astype(int)
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
                bedrock = np.logical_and(bedrock_r.astype(bool),
                                         bedrock_vadose.astype(bool))
                if self._is_defined['water_table']:
                    vadose_zone = self._get_phreatic_subdomain('vadose_zone')
                else:
                    vadose_zone = np.ones_like(bedrock)
                out = np.logical_and(bedrock.astype(bool),
                                     vadose_zone.astype(bool))
                out = out.astype(int)
        elif subdomain == 'bedrock_phreatic':
            # 'Bedrock' x 'Phreatic zone'
            bedrock = self._get_bedrock_subdomain()
            phreatic_zone = self._get_phreatic_subdomain('phreatic_zone')
            out = np.logical_and(bedrock.astype(bool),
                                 phreatic_zone.astype(bool))
            out = out.astype(int)
        else:
            print("ERROR: subdomain '{}' does not exist.".format(subdomain))
            out = None
        return out
    
    def _get_surface_subdomain(self, face_name: str) -> np.ndarray:
        """
        Compute the selected surfaces from the volumetric domain.
        """
        
        # retrieve grid indices
        domain = self.data_volume.astype('bool')
        I, J, K = np.indices(domain.shape)
        
        # retrieve z-surface grid indices
        if face_name in ['up', 'down']:
            i_, j_ = np.indices(self.data_surfaces['z'].shape)
            range_ = list(range(self.grid.nz))
            if face_name == 'up':
                face = K.max(axis=2, initial=-1, where=domain)
            elif face_name == 'down':
                face = K.min(axis=2, initial=self.grid.nz, where=domain)
        
        # flatten the indices
        i_ = i_.flatten()
        j_ = j_.flatten()
        k = face.flatten()
        
        # retrieve the valid i,j,k indices
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
        Compute the phreatic subdomain.
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
        Compute the bedrock subdomain.
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
        Return ``True`` if a (x, y, z)-point is inside the domain, otherwise
        ``False``.

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
     
    def is_2D_point_valid(self, point: tuple) -> bool:
        """
        Return ``True`` if a z-coordinate exists for a (x, y)-point projected
        inside the domain, otherwise ``False``.

        Parameters
        ----------
        point : tuple
            (x, y)-point

        Returns
        -------
        out : bool
        """
        point_object = Point(point)
        if self.grid.polygon.contains(point_object):
            i, j = self.grid.get_indices(point)
            out = bool(self.data_surfaces['z'][i, j])
        elif self.grid.polygon.touches(point_object):
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
