"""Custom pyKasso's typing."""

# Typing
from typing import TypeVar, TYPE_CHECKING

from numpy import random
from pandas import core, io

if TYPE_CHECKING:
    import pykasso.core.project as pcp
    import pykasso.core.grid as pcg
    import pykasso.model.domain as pcd
    import pykasso.model.geologic_features as pcgf

### Custom internal types

# Core types
Project = TypeVar('Project', bound='pcp.Project')
Grid = TypeVar('Grid', bound='pcg.Grid')

# Model types
Domain = TypeVar('Domain', bound='pcd.Domain')
Delimitation = TypeVar('Delimitation', bound='pcd.Delimitation')
Topography = TypeVar('Topography', bound='pcd.Topography')
Bedrock = TypeVar('Bedrock', bound='pcd.Bedrock')
WaterLevel = TypeVar('WaterLevel', bound='pcd.WaterLevel')
Geology = TypeVar('Geology', bound='pcgf.Geology')

### Custom external types
# Numpy
RandomNumberGenerator = TypeVar('RandomNumberGenerator',
                                bound='random._generator.Generator')
# Pandas
Series = TypeVar('Series', bound='core.series.Series')
DataFrame = TypeVar('DataFrame', bound='core.frame.DataFrame')
Styler = TypeVar('Styler', bound='io.formats.style.Styler')
