from ...core.grid import Grid
from ..geologic_features.surface import Surface


class Topography(Surface):
    """
    Class modeling the upper horizontal limit of the study site.
    """
    def __init__(
        self,
        grid: Grid,
        *args,
        **kwargs,
    ) -> None:
        feature = 'topography'
        super().__init__(grid, feature, *args, **kwargs)
