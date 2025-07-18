import sys
from enum import auto

if sys.version_info >= (3, 11):
    from enum import StrEnum
else:
    from yt._maintenance.backports import StrEnum


# register all valid geometries
class Geometry(StrEnum):
    CARTESIAN = auto()
    CYLINDRICAL = auto()
    POLAR = auto()
    SPHERICAL = auto()
    GEOGRAPHIC = auto()
    INTERNAL_GEOGRAPHIC = auto()
    SPECTRAL_CUBE = auto()

    def __str__(self):
        # Implemented for backward compatibility.
        if sys.version_info >= (3, 11):
            return super().__str__()
        else:
            return self.name.lower()
