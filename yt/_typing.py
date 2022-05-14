from typing import List, Optional, Tuple, Union

from numpy import ndarray
from unyt import unyt_quantity

FieldDescT = Tuple[str, Tuple[str, List[str], Optional[str]]]
KnownFieldsT = Tuple[FieldDescT, ...]

ParticleCoordinateTuple = Tuple[
    str,  # particle type
    Tuple[ndarray, ndarray, ndarray],  # xyz
    Union[float, ndarray],  # hsml
]

# an intentionally restrictive list of types that can
# be passes to ds.quan (which is a proxy for unyt.unyt_quantity.__init__)
Quantity = Union[unyt_quantity, Tuple[float, str]]
