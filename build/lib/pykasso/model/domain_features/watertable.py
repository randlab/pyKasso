from ...core.grid import Grid
from ..geologic_features.surface import Surface


class WaterTable(Surface):
    """
    Class modeling the water level elevation, the phreatic/vadose limit of the
    study site.
    """
    def __init__(
        self,
        grid: Grid,
        *args,
        **kwargs,
    ) -> None:
        feature = 'water_table'
        super().__init__(grid, feature, *args, **kwargs)
