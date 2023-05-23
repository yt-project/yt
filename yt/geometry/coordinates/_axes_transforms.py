import sys
from enum import auto
from typing import Optional

if sys.version_info >= (3, 11):
    from enum import StrEnum
else:
    from yt._maintenance.backports import StrEnum


class AxesTransform(StrEnum):
    DEFAULT = auto()
    GEOMETRY_NATIVE = auto()
    POLAR = auto()
    AITOFF_HAMMER = auto()


def parse_axes_transform(axes_transform: Optional[str]) -> AxesTransform:
    if axes_transform is None:
        # pass the responsability to ds.coordinates
        axes_transform = "default"
    return AxesTransform(axes_transform)
