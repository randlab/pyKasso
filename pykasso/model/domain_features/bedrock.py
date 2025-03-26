from ...core.grid import Grid
from ..geologic_features.surface import Surface


class Bedrock(Surface):
    """
    Class modeling the lower horizontal limit of the study site.
    """
    def __init__(
        self,
        grid: Grid,
        *args,
        **kwargs,
    ) -> None:
        feature = 'bedrock'
        super().__init__(grid, feature, *args, **kwargs)
