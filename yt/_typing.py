from typing import List, Optional, Tuple, Union

from numpy import ndarray

FieldDescT = Tuple[str, Tuple[str, List[str], Optional[str]]]
KnownFieldsT = Tuple[FieldDescT, ...]

ParticleCoordinateTuple = Tuple[
    str,  # particle type
    Tuple[ndarray, ndarray, ndarray],  # xyz
    Union[float, ndarray],  # hsml
]
