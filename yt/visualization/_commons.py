import os
import sys
import warnings
from typing import Optional, Type

import matplotlib
from packaging.version import Version

from ._mpl_imports import (
    FigureCanvasAgg,
    FigureCanvasBase,
    FigureCanvasPdf,
    FigureCanvasPS,
    FigureCanvasSVG,
)

MPL_VERSION = Version(matplotlib.__version__)

DEFAULT_FONT_PROPERTIES = {
    "family": "stixgeneral",
    "size": 18,
}

if MPL_VERSION >= Version("3.4"):
    DEFAULT_FONT_PROPERTIES["math_fontfamily"] = "cm"

SUPPORTED_FORMATS = frozenset(FigureCanvasBase.get_supported_filetypes().keys())
SUPPORTED_CANVAS_CLASSES = frozenset(
    (FigureCanvasAgg, FigureCanvasPdf, FigureCanvasPS, FigureCanvasSVG)
)


def get_canvas_class(suffix: str) -> Type[FigureCanvasBase]:
    s = normalize_extension_string(suffix)
    if s not in SUPPORTED_FORMATS:
        raise ValueError(f"Unsupported file format '{suffix}'.")
    for cls in SUPPORTED_CANVAS_CLASSES:
        if s in cls.get_supported_filetypes():
            return cls
    raise RuntimeError(
        "Something went terribly wrong. "
        f"File extension '{suffix}' is supposed to be supported "
        "but no compatible backend was found."
    )


def normalize_extension_string(s: str) -> str:
    if sys.version_info < (3, 9):
        if s.startswith("."):
            return s[1:]
        return s
    else:
        return s.removeprefix(".")


def validate_image_name(filename, suffix: Optional[str] = None) -> str:
    """
    Build a valid image filename with a specified extension (default to png).
    The suffix parameter is ignored if the input filename has a valid extension already.
    Otherwise, suffix is appended to the filename, replacing any existing extension.
    """
    name, psuffix = os.path.splitext(filename)
    psuffix = normalize_extension_string(psuffix)

    if suffix is not None:
        suffix = normalize_extension_string(suffix)

    if psuffix in SUPPORTED_FORMATS:
        if suffix in SUPPORTED_FORMATS and suffix != psuffix:
            warnings.warn(
                f"Received two valid image formats {psuffix!r} (from filename) "
                f"and {suffix!r} (from suffix). The former is ignored."
            )
            return f"{name}.{suffix}"
        return str(filename)

    if suffix is None:
        suffix = "png"

    if suffix not in SUPPORTED_FORMATS:
        raise ValueError(f"Unsupported file format {suffix!r}")

    return f"{filename}.{suffix}"


def get_canvas(figure, filename):

    name, suffix = os.path.splitext(filename)

    if not suffix:
        raise ValueError(
            f"Can not determine canvas class from filename '{filename}' "
            f"without an extension."
        )
    return get_canvas_class(suffix)(figure)
