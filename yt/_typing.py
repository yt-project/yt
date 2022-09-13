from typing import List, Optional, Tuple, Union

import unyt as un
from numpy import ndarray

FieldDescT = Tuple[str, Tuple[str, List[str], Optional[str]]]
KnownFieldsT = Tuple[FieldDescT, ...]

ParticleCoordinateTuple = Tuple[
    str,  # particle type
    Tuple[ndarray, ndarray, ndarray],  # xyz
    Union[float, ndarray],  # hsml
]


# types that can be converted to un.Unit
Unit = Union[un.Unit, str]

# types that can be converted to un.unyt_quantity
Quantity = Union[un.unyt_quantity, Tuple[float, Unit]]
