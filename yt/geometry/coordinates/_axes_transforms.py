from enum import Enum, auto
from typing import Optional


class AxesTransform(Enum):
    DEFAULT = auto()
    GEOMETRY_NATIVE = auto()
    POLAR = auto()
    AITOFF_HAMMER = auto()


def parse_axes_transform(axes_transform: Optional[str]) -> AxesTransform:
    if axes_transform is None:
        # pass the responsability to ds.coordinates
        return AxesTransform.DEFAULT
    elif axes_transform == "geometry_native":
        return AxesTransform.GEOMETRY_NATIVE
    elif axes_transform == "polar":
        return AxesTransform.POLAR
    elif axes_transform == "aitoff_hammer":
        return AxesTransform.AITOFF_HAMMER
    else:
        raise ValueError(f"Unknown axes transform {axes_transform!r}")
