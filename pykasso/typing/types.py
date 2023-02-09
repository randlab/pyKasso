# Typing
# from typing import TYPE_CHECKING, List, Any, Callable, TypeVar, Union
from typing import TypeVar

# if TYPE_CHECKING:  # https://mypy.readthedocs.io/en/latest/runtime_troubles.html
    # import pykasso as pk

# from pykasso.core.grid import Grid
# from pykasso.core.domain import Delimitation, Topography, Bedrock

import pykasso.core.grid              as pcg
import pykasso.core.domain            as pcd
import pykasso.core.geologic_features as pcgf

# Internal Types
Grid         = TypeVar('Grid'        , bound='pcg.Grid')
Domain       = TypeVar('Domain'      , bound='pcd.Domain')
Delimitation = TypeVar('Delimitation', bound='pcd.Delimitation')
Topography   = TypeVar('Topography'  , bound='pcd.Topography')
Bedrock      = TypeVar('Bedrock'     , bound='pcd.Bedrock')
WaterLevel   = TypeVar('WaterLevel'  , bound='pcd.WaterLevel')
Geology      = TypeVar('Geology'     , bound='pcgf.Geology')