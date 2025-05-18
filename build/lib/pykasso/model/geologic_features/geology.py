from ...core._namespaces import DEFAULT_FMM_COSTS
from ...core.grid import Grid
from .geologicfeature import GeologicFeature


class Geology(GeologicFeature):
    """
    Class modeling the geologic model.
    """
    
    def __init__(
        self,
        grid: Grid,
        default_fmm_cost: float = DEFAULT_FMM_COSTS['geology'],
        *args,
        **kwargs,
    ) -> None:
        feature = 'geology'
        dim = 3
        super().__init__(grid, feature, dim, default_fmm_cost, *args, **kwargs)
        
    def set_names(
        self,
        names: dict[int, str],
        default_name: str = 'unit {}',
    ) -> None:
        """
        Assign names to geologic units based on the provided ``names``
        dictionary, with an optional default naming pattern.

        Parameters
        ----------
        names : dict[int, str]
            A dictionary where the keys are geologic unit indices (integers)
            and the values are the corresponding names (strings) to be
            assigned. This dictionary specifies which geologic unit should
            receive custom names.
        default_name : str, default: 'unit {}'
            A format string used to generate default geologic unit names for
            items not explicitly named in the ``names`` dictionary. The format
            string should include a placeholder (e.g., '{}') that will be
            replaced by the item's index.
        
        Notes
        -----
        This function does not return a value. It rewrites the ``self.names``
        attribute with the new specified dictionary.
        """
        return super().set_names(names, default_name)
    
    def set_model(
        self,
        model: dict[int, str],
        default_model: bool = True,
    ) -> None:
        """
        Indicate if a geologic unit should be considered in the modelisation
        based on the provided dictionary, with an optional default setting.

        Parameters
        ----------
        model : dict[int, bool]
            A dictionary where the keys are geologic unit indices (integers)
            and the values are booleans indicating if the item is considered or
            not.
        default_model : bool, default: True
            The default value to be applied to geologic units not explicitly
            listed in the ``model`` dictionary.
        
        Notes
        -----
        This function does not return a value. It rewrites the ``self.model``
        attribute with the new specified dictionary.
        """
        model.setdefault(0, True)
        return super().set_model(model, default_model)
