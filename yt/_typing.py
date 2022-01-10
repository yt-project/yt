from typing import DefaultDict, List, Optional, Tuple, Union

from numpy import ndarray
from numpy.typing import ArrayLike

FieldDescT = Tuple[str, Tuple[str, List[str], Optional[str]]]
KnownFieldsT = Tuple[FieldDescT, ...]

# the following custom types are related to particle io
ParticleCoordinateTuple = Tuple[
    str,  # particle type
    Tuple[ndarray, ndarray, ndarray],  # xyz
]

SPHParticleCoordinateTuple = Tuple[
    str,  # particle type
    Tuple[ndarray, ndarray, ndarray],  # xyz
    Union[float, ArrayLike, ndarray],  # hsml
]

ParticleTypeSizes = DefaultDict[str, int]
