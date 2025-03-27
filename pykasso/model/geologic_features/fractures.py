import mpmath
import numpy as np
import pandas as pd

from ..fracturation import (
    voxelize_fractures,
    _calculate_normal,
)
from .geologicfeature import GeologicFeature
from ...core._namespaces import DEFAULT_FMM_COSTS

from ...core.grid import Grid
from numpy.random import Generator
from typing import Union


class Fractures(GeologicFeature):
    """
    Class modeling the fracturation model.
    """
    
    def __init__(
        self,
        grid: Grid,
        rng: Generator,
        *args,
        **kwargs,
    ) -> None:
        """
        Class modeling the fracturation model.

        Parameters
        ----------
        grid : Grid
            pyKasso's ``Grid`` of the model.
        rng : Generator
            Random Generator Number of numpy.
        """
        # Set super constructor parameters
        label = 'fractures'
        dim = 3
        super().__init__(grid, label, dim, *args, **kwargs)
        
        # Set class attributes
        self.rng = rng
        self.i = self.stats.index.max()
        self.families = pd.DataFrame()
        self.fractures = pd.DataFrame()
        self.fractures_voxelized = {}
        
    def set_names(
        self,
        names: dict[int, str],
        default_name: str = 'family {}',
    ) -> None:
        """
        Assign names to fracture families based on the provided ``names``
        dictionary, with an optional default naming pattern.

        Parameters
        ----------
        names : dict[int, str]
            A dictionary where the keys are fracture familiy indices (integers)
            and the values are the corresponding names (strings) to be
            assigned. This dictionary specifies which fracture families should
            receive custom names.
        default_name : str, default: 'family {}'
            A format string used to generate default fracture familiy names
            for items not explicitly named in the ``names`` dictionary.
            The format string should include a placeholder (e.g., '{}') that
            will be replaced by the item's index.
        
        Notes
        -----
        This function does not return a value. It rewrites the ``self.names``
        attribute with the new specified dictionary.
        """
        return super().set_names(names, default_name)

    def set_costs(
        self,
        costs: dict[int, str],
        default_cost: float = DEFAULT_FMM_COSTS['fractures'],
    ) -> None:
        """
        Assign costs to fracture families based on the provided dictionary,
        with an optional default cost.

        Parameters
        ----------
        costs : dict[int, str]
            A dictionary where the keys are fracture familiy indices (integers)
            and the values are the corresponding costs (floats) to be assigned.
            This dictionary specifies which fracture families should receive
            custom costs.
        default_cost : float, optional
            The default cost to be applied to fracture familes not explicitly
             listed in the ``costs`` dictionary. The default values are taken
             from the ``DEFAULT_FMM_COSTS['fractures']`` dictionary.
        
        Notes
        -----
        This function does not return a value. It rewrites the ``self.costs``
        attribute with the new specified dictionary.
        """
        return super().set_costs(costs, default_cost)
    
    def set_model(
        self,
        model: dict[int, bool],
        default_model: bool = True,
    ) -> None:
        """
        Indicate if a fracture family should be considered in the modelisation
        based on the provided dictionary, with an optional default setting.

        Parameters
        ----------
        model : dict[int, bool]
            A dictionary where the keys are fracture familiy indices (integers)
             and the values are booleans indicating if the fracture family is
            considered. This dictionary specifies which fracture families
             should be evaluated in the simulation.
            
        default_model : bool, default: True
            The default value to be applied to fracture families not
             explicitly listed in the ``model`` dictionary.
        
        Notes
        -----
        This function does not return a value. It rewrites the ``self.model``
        attribute with the new specified dictionary.
        """
        
        model.setdefault(0, False)
        return super().set_model(model, default_model)
        
    def generate_fracture_family(
        self,
        name: str,
        settings: dict,
        cost: float = DEFAULT_FMM_COSTS['fractures'],
    ) -> None:
        """
        Create a new fracture family.
        Populate the ``self.families`` and ``self.fractures`` dataframe attributes.

        Parameters
        ----------
        name : str
            Name of the fracture family.
        settings : dict
            Parameters for the fracture generation.
        cost : float, optional
            Travel cost of the fracture family. By default the default travel
            cost value for fractures.
        """
        # Update fracture families counter
        self.i = self.i + 1
        
        # Create a dataframe storing the new family settings
        frac_family = {
            'name': [name],
            'density': settings['density'],
            'cost': cost,
            # 'geology': [0,1,2] # TODO - add geology
        }
        frac_family = pd.DataFrame(data=frac_family, index=[self.i])
        
        # Generate fractures
        fractures = self.generate_fractures(**settings)
        fractures.insert(0, 'family_id', self.i)
        
        # Update tables
        self.fractures = pd.concat([self.fractures, fractures])
        self.families = pd.concat([self.families, frac_family])
        
        # Update dicts
        self.names |= self.families['name'].to_dict()
        self.costs |= self.families['cost'].to_dict()
        index = self.families.index.to_list()
        self.model |= {i: True for i in index}
        return None
     
    def compute_model(self) -> None:
        """
        Construct the array for the fracturation model by voxelizing the
        fractures.
        
        Voxelize each fracture family and populate the ``fractures_voxelized``
        dictionary attribute. Sort fracture families by travel cost and
        combine all the voxels according to their ranking.
        """
        # Voxelize the fractures
        fractures_voxelized = {}
        for i in self.families.index:
            fractures = self.fractures[self.fractures['family_id'] == i]
            if fractures.empty:
                fractures_voxelized[i] = np.zeros_like(self.grid.data_volume)
            else:
                fractures_voxelized[i] = voxelize_fractures(self.grid,
                                                            fractures)
                
        # Sort fracture families by cost
        self.families = self.families.sort_values(['cost'])
        fractures_family_ids = self.families.index
        
        # Construct the fracturation array
        frac_model = np.zeros_like(self.grid, dtype=np.float64)
        for family_id in fractures_family_ids:
            frac_model = np.where(
                fractures_voxelized[family_id] == 1,
                family_id,
                frac_model
            )
    
        # Update the attributes
        self.data_volume = frac_model.copy()
        self.fractures_voxelized = fractures_voxelized
                
        return None
                
    def generate_fractures(
        self,
        density: float,
        orientation: Union[list, float, int],
        dip: Union[list, float, int],
        length: Union[list, float, int],
        orientation_distribution: str = 'vonmises',
        dip_distribution: str = 'vonmises',
        length_distribution: str = 'power',
        **kwargs: dict,
    ) -> pd.DataFrame:

        ######################
        ### INITIALIZATION ###
        ######################
        
        ### Set optional parameters
        alpha = kwargs.get('alpha', 2)  # TODO : Default value ??
        
        ### Set angle parameters
        orientation_min = 0
        orientation_max = 360
        dip_min = 0
        dip_max = 360
        
        ### Define if parameter is a range or an unique value
        if isinstance(orientation, (list)):
            orientation_min, orientation_max = orientation
        else:
            orientation_distribution = 'unique'
        
        if isinstance(dip, (list)):
            dip_min, dip_max = dip
        else:
            dip_distribution = 'unique'
        
        if isinstance(length, (list)):
            length_min, length_max = min(length), max(length)
        else:
            length_distribution = 'unique'
            length_max = length

        ### Redefine fracturation domain
        Lx = self.grid.xmax - self.grid.xmin
        Ly = self.grid.ymax - self.grid.ymin
        Lz = self.grid.zmax - self.grid.zmin

        shift_x = min(Lx / 2, length_max / 2)
        shift_y = min(Ly / 2, length_max / 2)
        shift_z = min(Lz / 2, length_max / 2)

        Lex = 2 * shift_x + Lx
        Ley = 2 * shift_y + Ly
        Lez = 2 * shift_z + Lz

        area = Lex * Ley
        xmin = self.grid.xmin - shift_x
        ymin = self.grid.ymin - shift_y
        zmin = self.grid.zmin - shift_z

        ### Total numbers of fractures
        fractures_numbers = np.array(density) * area

        ### Generate poisson number for each fracture family
        real_frac_number = self.rng.poisson(fractures_numbers)
        
        ###########################
        ### Calculate fractures ###
        ###########################

        ##### FRACTURE CENTER LOCATION

        # Get fracture position from triple uniform distribution
        xm = Lex * self.rng.random(size=real_frac_number) + xmin
        ym = Ley * self.rng.random(size=real_frac_number) + ymin
        zm = Lez * self.rng.random(size=real_frac_number) + zmin

        ##### FRACTURE ORIENTATION
        
        if orientation_min > orientation_max:
            orientation_min = orientation_min - 360
        
        # No distribution case
        if orientation_distribution == 'unique':
            orientation_fractures = [orientation] * real_frac_number
        
        # Uniform distribution case
        elif orientation_distribution == 'uniform':
            orientation_fractures = self._uniform(orientation_min,
                                                  orientation_max,
                                                  real_frac_number)
        
        # Von Mises distribution case
        elif orientation_distribution == 'vonmises':
            orientation_fractures = self._vonmises(orientation_min,
                                                   orientation_max,
                                                   real_frac_number)
        
        ##### FRACTURE DIP
        if dip_min > dip_max:
            dip_min = dip_min - 360

        # No distribution case
        if dip_distribution == 'unique':
            dip_fractures = [dip] * real_frac_number
        
        # Uniform distribution case
        elif dip_distribution == 'uniform':
            dip_fractures = self._uniform(dip_min,
                                          dip_max,
                                          real_frac_number)
            
        # Von Mises distribution case
        elif dip_distribution == 'vonmises':
            dip_fractures = self._vonmises(dip_min,
                                           dip_max,
                                           real_frac_number)

        ##### FRACTURE LENGHT

        # No distribution case
        if length_distribution == 'unique':
            radius = [length / 2] * real_frac_number

        # Uniform distribution case
        elif length_distribution == 'uniform':
            radius = self._uniform(length_min,
                                   length_max,
                                   real_frac_number)

        # Trucated power law distribution case
        elif length_distribution == 'power':
            frac_length = self._power(length_min, length_max,
                                      alpha, real_frac_number)
            radius = frac_length / 2
            
        ##### FRACTURE NORMAL VECTOR
        normal = _calculate_normal(dip_fractures, orientation_fractures)

        #######################
        ### Returns results ###
        #######################
        
        data = {
            'x': xm,
            'y': ym,
            'z': zm,
            'radius': radius,
            'orientation': orientation_fractures,
            'dip': dip_fractures,
            'normal': normal
        }
        fractures = pd.DataFrame(data=data)

        return fractures
    
    def _uniform(
        self,
        value_min: float,
        value_max: float,
        n: int,
    ) -> np.ndarray:
        """
        Generate an array of ``n`` uniformly distributed random values between
        specified minimum and maximum bounds.
        """
        out = self.rng.uniform(low=value_min, high=value_max, size=n)
        return out
    
    @staticmethod
    def _vonmises_calculate_mu(
        theta_min: float,
        theta_max: float,
    ) -> float:
        """
        Calculate the mu parameter from the von Mises distribution.
        """
        mu = (theta_min + theta_max) / 2
        return mu
    
    @staticmethod
    def _vonmises_calculate_kappa(
        theta_max: float,
        mu: float,
    ) -> float:
        """
        Calculate the kappa parameter from the von Mises distribution.
        """
        std = (theta_max - mu) / 3
        kappa = Fractures._vonmises_solve_kappa(std)
        return kappa

    @staticmethod
    def _vonmises_solve_kappa(std: float):
        """
        Solve the variance equation for von Mises distribution by finding the
        corresponding kappa.
        """
        def func(kappa):
            a = std**2 - 1
            b = mpmath.besseli(1, kappa) / mpmath.besseli(0, kappa)
            return a + b
        kappa = mpmath.findroot(func, 0)
        return kappa
    
    def _vonmises(
        self,
        theta_min: float,
        theta_max: float,
        n: int,
    ) -> np.ndarray:
        """
        Generate an array of ``n`` random values following a von Mises
        distribution in an interval specified by minimum and maximum bounds.
        More details here : https://en.wikipedia.org/wiki/Von_Mises_distribution
        """
        theta_min_rad = np.radians(theta_min)
        theta_max_rad = np.radians(theta_max)
        mu = Fractures._vonmises_calculate_mu(theta_min_rad, theta_max_rad)
        kappa = Fractures._vonmises_calculate_kappa(theta_max_rad, mu)
        out = self.rng.vonmises(mu, kappa, size=n)
        out = np.degrees(out)
        return out
    
    def _power(
        self,
        value_min: float,
        value_max: float,
        alpha: float,
        n: int,
    ) -> np.ndarray:
        """
        Generate an array of ``n`` random values following a modified power-law
        distribution in an interval specified by minimum and maximum bounds.
        """
        palpha = (1 - alpha)
        invpalpha = 1 / palpha
        fmina = value_min**palpha
        frangea = value_max**palpha - fmina
        u = self.rng.random(size=n)
        out = (fmina + u * frangea)**invpalpha
        return out
