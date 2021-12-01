import os
import sys
from typing import Optional, Tuple, Type, TypeVar

import matplotlib
from packaging.version import Version

from yt.utilities.logger import ytLogger as mylog

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
            mylog.warning(
                "Received two valid image formats '%s' (from `filename`) "
                "and '%s' (from `suffix`). The former is ignored.",
                psuffix,
                suffix,
            )
            return f"{name}.{suffix}"
        return str(filename)

    if suffix is None:
        suffix = "png"

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


T = TypeVar("T", tuple, list)


def _swap_axes_extents(extent: T) -> T:
    """
    swaps the x and y extent values, preserving type of extent

    Parameters
    ----------
    extent : sequence of four unyt quantities
        the current 4-element tuple or list of unyt quantities describing the
        plot extent. extent = (xmin, xmax, ymin, ymax).

    Returns
    -------
    tuple or list
        the extent axes swapped, now with (ymin, ymax, xmin, xmax).

    """
    extent_swapped = [extent[2], extent[3], extent[0], extent[1]]
    return type(extent)(extent_swapped)


def _swap_arg_pair_order(*args):
    """
    flips adjacent argument pairs, useful for swapping x-y plot arguments

    Parameters
    ----------
    *args
        argument pairs, must have an even number of *args

    Returns
    -------
    tuple
        args  with order of pairs switched, i.e,:

        _swap_arg_pair_order(x, y, px, py) returns:
            y, x, py, px

    """

    assert len(args) % 2 == 0
    n_pairs = len(args) // 2
    new_args = []
    for i in range(n_pairs):
        x_id = i * 2
        new_args.append(args[x_id + 1])
        new_args.append(args[x_id])
    return tuple(new_args)
