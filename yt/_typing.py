import os
from typing import Tuple, TypeVar, Union

from matplotlib.colors import ListedColormap
from unyt.unit_object import Unit

# this defines a useful generic type to maintain type information
# across generic functions and parametrization
# It is close to the concept of templates in C++
# see https://docs.python.org/3/library/typing.html#generics
T = TypeVar("T")

FileLike = Union[str, bytes, "os.PathLike[str]", int]

# Note that Tuple is currently the only Sequence-like type
# that support specifying the size
# In most places where a field-like input is expected,
# List[str] should be acceptable as long as the size is exactly 2.
#
# here's a proto adapter function in python>=3.10
# def unpack_field(field:FieldLike) -> Tuple[str, str]:
#     match field:
#         case ftype, fname:
#             pass
#         case str():
#             ftype = "auto"
#             fname = field
#         case _:
#             raise TypeError
#     return ftype, fname

# and in python>=3.6 (only type hints are incompatible with earlier versions of Python 3)
# def unpack_field(field:FieldLike) -> Tuple[str, str]:
#     from yt.funcs import is_sequence
#     if is_sequence(field) and len(field) == 2:
#         return field
#     elif isinstance(field, str):
#         return "auto", field
#     raise TypeError
FieldLike = Union[Tuple[str, str], str]
UnitLike = Union[Unit, str]
CmapLike = Union[ListedColormap, str]
