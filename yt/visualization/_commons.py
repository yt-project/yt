import os
import sys
from typing import Optional, Type

from yt.utilities._version import MPL_VERSION
from yt.utilities.logger import ytLogger as mylog

from ._mpl_imports import (
    FigureCanvasAgg,
    FigureCanvasBase,
    FigureCanvasPdf,
    FigureCanvasPS,
    FigureCanvasSVG,
)

DEFAULT_FONT_PROPERTIES = {
    "family": "stixgeneral",
    "size": 18,
}

if MPL_VERSION >= (3, 4, 0):
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
    if normalize_extension_string(psuffix) in SUPPORTED_FORMATS:
        if suffix is not None:
            suffix = normalize_extension_string(suffix)
        if suffix in SUPPORTED_FORMATS and suffix != psuffix:
            mylog.warning(
                "Received two valid image formats '%s' (from `filename`) "
                "and '%s' (from `suffix`). The former is ignored.",
                psuffix,
                suffix,
            )
            return f"{name}.{suffix}"
        return str(filename)

    if suffix is None:
        suffix = ".png"

    suffix = normalize_extension_string(suffix)

    if suffix not in SUPPORTED_FORMATS:
        raise ValueError(f"Unsupported file format '{suffix}'.")

    return f"{filename}.{suffix}"


def get_canvas(figure, filename):

    name, suffix = os.path.splitext(filename)

    if not suffix:
        raise ValueError(
            f"Can not determine canvas class from filename '{filename}' "
            f"without an extension."
        )
    return get_canvas_class(suffix)(figure)
