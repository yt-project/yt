from enum import StrEnum, auto


# register all valid geometries
class Geometry(StrEnum):
    CARTESIAN = auto()
    CYLINDRICAL = auto()
    POLAR = auto()
    SPHERICAL = auto()
    GEOGRAPHIC = auto()
    INTERNAL_GEOGRAPHIC = auto()
    SPECTRAL_CUBE = auto()
