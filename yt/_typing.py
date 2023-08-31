from typing import Optional, Union

import unyt as un
from numpy import ndarray

FieldDescT = tuple[str, tuple[str, list[str], Optional[str]]]
KnownFieldsT = tuple[FieldDescT, ...]

ParticleType = str
FieldType = str
FieldName = str
FieldKey = tuple[FieldType, FieldName]
ImplicitFieldKey = FieldName
AnyFieldKey = Union[FieldKey, ImplicitFieldKey]
DomainDimensions = Union[tuple[int, ...], list[int], ndarray]

ParticleCoordinateTuple = tuple[
    str,  # particle type
    tuple[ndarray, ndarray, ndarray],  # xyz
    Union[float, ndarray],  # hsml
]

# Geometry specific types
AxisName = str
AxisOrder = tuple[AxisName, AxisName, AxisName]

# types that can be converted to un.Unit
Unit = Union[un.Unit, str]

# types that can be converted to un.unyt_quantity
Quantity = Union[un.unyt_quantity, tuple[float, Unit]]
