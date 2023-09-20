from typing import Any, Optional, Union

import numpy as np
import unyt as un

FieldDescT = tuple[str, tuple[str, list[str], Optional[str]]]
KnownFieldsT = tuple[FieldDescT, ...]

ParticleType = str
FieldType = str
FieldName = str
FieldKey = tuple[FieldType, FieldName]
ImplicitFieldKey = FieldName
AnyFieldKey = Union[FieldKey, ImplicitFieldKey]
DomainDimensions = Union[tuple[int, ...], list[int], np.ndarray]

ParticleCoordinateTuple = tuple[
    str,  # particle type
    tuple[np.ndarray, np.ndarray, np.ndarray],  # xyz
    Union[float, np.ndarray],  # hsml
]

# Geometry specific types
AxisName = str
AxisOrder = tuple[AxisName, AxisName, AxisName]

# types that can be converted to un.Unit
Unit = Union[un.Unit, str]

# types that can be converted to un.unyt_quantity
Quantity = Union[un.unyt_quantity, tuple[float, Unit]]

# np.ndarray[...] syntax is runtime-valid from numpy 1.22, we quote it until our minimal
# runtime requirement is bumped to, or beyond this version

MaskT = Optional["np.ndarray[Any, np.dtype[np.bool_]]"]
AlphaT = Optional["np.ndarray[Any, np.dtype[np.float64]]"]
