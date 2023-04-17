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


class Domain():
    """
    Class modeling the spatial extension of the domain within the study site
    grid locating where karstic network generation is allowed.
    """
    
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
        
        ### Computes subdomains from the lower and upper surfaces of the domain
        self._set_faces()
        
        ### Computes subdomains from the borders of the domain
        self._set_border_domains()
        
        ### Computes subdomains from water level elevation model
        self._set_phreatic_domains()
        
        ### Computes subdomains from bedrock elevation model
        self._set_bedrock_domains()
            
        ### Names the subdomains
        self.subdomains_names = {
            # domain
            "domain": "self.data_volume",
            "domain_surface": "self.faces['up']",
            "domain_bottom": "self.faces['down']",
            # borders
            "domain_borders": "self.borders['domain_full']",
            "domain_borders_sides": "self.borders['domain_sides']",
            "domain_borders_surface": "self.borders['domain_up']",
            "domain_borders_bottom": "self.borders['domain_down']",
            # phreatic
            "vadose": "self.phreatic['vadose_zone']",
            "vadose_borders": "self.borders['vadose_zone']",
            "phreatic": "self.phreatic['phreatic_zone']",
            "phreatic_surface": "self.phreatic['water_level_surface']",
            "phreatic_borders_surface": "self.borders['water_level_surface']",
        }
    
    def _set_faces(self) -> None:
        """Computes the upper and lower surfaces from the volumetric domain."""
        self.faces = {}  # faces will be stored in this dict attribute
        domain = self.data_volume.astype('bool')
        I, J, K = np.indices(self.data_volume.shape)  # retrieves grid indices
        
        # for face in ['up', 'down', 'east', 'west', 'north', 'soouth']:
        for face_name in ['up', 'down']:
            if face_name in ['up', 'down']:
                # retrieves z-surface grid indices
                i_, j_ = np.indices(self.data_surfaces['z'].shape)
                range_ = list(range(self.grid.nz))
                if face_name == 'up':
                    face = K.max(axis=2, initial=-1, where=domain)
                else:
                    face = K.min(axis=2, initial=self.grid.nz, where=domain)
        
            i_, j_ = i_.flatten(), j_.flatten()
            k = face.flatten()
            nodes = list(zip(i_, j_, k))
            valid_nodes = [(i, j, k) for (i, j, k) in nodes if k in range_]
            i, j, k = zip(*valid_nodes)
            volume = np.zeros_like(self.data_volume)
            volume[i, j, k] = 1
            self.faces[face_name] = volume
        return None

    def _set_border_domains(self) -> None:
        """Computes the borders subdomains."""
        self.borders = {}
        
        # Computes the full borders extent using scipy.binary_dilation
        # algorithm
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.binary_dilation.html
        k = np.zeros((3, 3, 3), dtype=int)
        k[:, 1, 1], k[1, :, 1], k[1, 1, :] = 1, 1, 1
        dilation = binary_dilation(self.data_volume == 0, k, border_value=1)
        self.borders['domain_full'] = dilation & self.data_volume
        
        # x and y-axis borders
        k = np.zeros((3, 3), dtype=int)
        k[:, 1], k[1, :] = 1, 1
        self.borders['domain_sides'] = np.zeros_like(self.data_volume)
        for z in range(self.grid.nz):
            dilation = binary_dilation(self.data_volume[:, :, z] == 0, k,
                                       border_value=1)
            test = dilation & self.data_volume[:, :, z]
            self.borders['domain_sides'][:, :, z] = test
            
        # 'Domain' x 'Face up'
        self.borders['domain_up'] = (
            np.logical_and(self.borders['domain_sides'], self.faces['up'])
        )
        
        # 'Domain' x 'Face down'
        self.borders['domain_down'] = (
            np.logical_and(self.borders['domain_sides'], self.faces['down'])
        )
        return None
        
    def _set_phreatic_domains(self) -> None:
        """Computes the subdomains derivated from the water level elevation."""
        if self.water_level is None:
            vadose_zone = self.data_volume.copy()
            phreatic_zone = np.zeros_like(self.grid.data_volume)
            water_level_surface = np.zeros_like(self.grid.data_volume)
            vadose_zone_borders = self.borders['domain_sides']
            phreatic_zone_borders = np.zeros_like(self.grid.data_volume)
            wls_borders = np.zeros_like(self.grid.data_volume)
        else:
            water_volume = self.water_level.data_volume
            water_surface = self.water_level._surface_to_volume('=', self.grid)
            
            vadose_zone = np.logical_and(self.data_volume,
                                         np.logical_not(water_volume))
            phreatic_zone = np.logical_and(self.data_volume, water_volume)
            water_level_surface = (
                np.logical_and(self.data_volume, water_surface)
            )
            
            # Water level surface border
            vadose_zone_borders = np.logical_and(
                self.borders['domain_sides'], vadose_zone)
            phreatic_zone_borders = np.logical_and(
                self.borders['domain_sides'], phreatic_zone)
            wls_borders = np.logical_and(
                self.borders['domain_sides'], water_level_surface)
        
        # Sets the new domains
        self.phreatic = {
            'vadose_zone': vadose_zone,
            'phreatic_zone': phreatic_zone,
            'water_level_surface': water_level_surface,
            
        }
        self.borders['vadose_zone'] = vadose_zone_borders
        self.borders['phreatic_zone'] = phreatic_zone_borders
        self.borders['water_level_surface'] = wls_borders
        
        return None
    
    def _set_bedrock_domains(self) -> None:
        """ """
        if self.bedrock is None:
            bedrock = np.zeros_like(self.grid.data_volume)
            bedrock_vadose = np.zeros_like(self.grid.data_volume)
            bedrock_phreatic = np.zeros_like(self.grid.data_volume)
        else:
            roll_value = 2
            bedrock = (
                np.roll(self.bedrock.data_volume, roll_value, axis=2)
            )
            bedrock_vadose = np.logical_and(bedrock,
                                            self.phreatic['vadose_zone'])
            bedrock_phreatic = np.logical_and(bedrock,
                                              self.phreatic['phreatic_zone'])
        
        # Sets the new domains
        self.bedrock_domains = {
            'bedrock': bedrock,
            'bedrock_vadose': bedrock_vadose,
            'bedrock_phreatic': bedrock_phreatic,
        }
        return None
    
    ###############
    ### METHODS ###
    ###############
    
    def get_subdomain(self, subdomain_name: str) -> np.ndarray:
        """Returns the numpy.ndarray modeling the requested subdomain.

        Parameters
        ----------
        subdomain_name : str
            Name of the requested subdomain: {}
            
        Returns
        -------
        out : np.ndarray
        """
        out = eval(self.subdomains_names[subdomain_name])
        return out

    def is_coordinate_in_domain(self, x: float, y: float, z: float) -> bool:
        """
        Returns true if a (x, y, z)-coordinate point is inside the domain,
        otherwise false.

        Parameters
        ----------
        x : float
            x-coordinate.
        y : float
            y-coordinate.
        z : float
            z-coordinate.

        Returns
        -------
        out : bool
        """
        if self.grid.is_inbox(x, y, z):
            i, j, k = self.grid.get_indices(x, y, z)
            out = bool(self.data_volume[i, j, k])
            return out
        else:
            out = False
            return out
        
    def is_coordinate_in_subdomain(self, subdomain: str, x: float, y: float,
                                   z: float) -> bool:
        """
        Returns true if a (x, y, z)-coordinate point is inside the subdomain,
        otherwise false.

        Parameters
        ----------
        subdomain : str
            Name of the subdomain to test: {}
        x : float
            x-coordinate.
        y : float
            y-coordinate.
        z : float
            z-coordinate.

        Returns
        -------
        out : bool
        """
        subdomain = eval(self.subdomains_names[subdomain])
        if self.grid.is_inbox(x, y, z):
            i, j, k = self.grid.get_indices(x, y, z)
            out = bool(subdomain[i, j, k])
            return out
        else:
            out = False
            return out
        
    def is_coordinate_2D_valid(self, x: float, y: float) -> bool:
        """
        Returns true if a z-coordinate exists for a (x, y)-coordinate point
        projected inside the domain, otherwise false.

        Parameters
        ----------
        x : float
            x-coordinate.
        y : float
            y-coordinate.

        Returns
        -------
        out : bool
        """
        point = Point(x, y)
        if self.grid.polygon.contains(point):
            i, j = self.grid.get_indices(x, y)
            out = bool(self.data_surfaces['z'][i, j])
        elif self.grid.polygon.touches(point):
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
        pts = np.column_stack((grid.X[row, col, 0].flatten(),
                               grid.Y[row, col, 0].flatten()))
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

subdomains_list = """
            - "domain"
            - "domain_surface"
            - "domain_bottom"
            - "domain_borders"
            - "domain_borders_sides"
            - "domain_borders_surface"
            - "domain_borders_bottom"
            - "vadose"
            - "vadose_borders"
            - "phreatic"
            - "phreatic_surface"
            - "phreatic_borders_surface"
"""

# Updates documentation
Domain.get_subdomain.__doc__ = (Domain.get_subdomain.__doc__
                                .format(subdomains_list))
Domain.is_coordinate_in_subdomain.__doc__ = (
    Domain.is_coordinate_in_subdomain.__doc__.format(subdomains_list)
)
