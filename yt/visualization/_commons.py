import os
import sys
import warnings
from functools import wraps
from typing import TYPE_CHECKING, Optional, Type, TypeVar

if sys.version_info >= (3, 8):
    from importlib.metadata import version
else:
    from importlib_metadata import version

from packaging.version import Version

if TYPE_CHECKING:
    from ._mpl_imports import FigureCanvasBase


DEFAULT_FONT_PROPERTIES = {
    "family": "stixgeneral",
    "size": 18,
}

MPL_VERSION = Version(version("matplotlib"))

if MPL_VERSION >= Version("3.4"):
    DEFAULT_FONT_PROPERTIES["math_fontfamily"] = "cm"


def _get_supported_image_file_formats():
    from ._mpl_imports import FigureCanvasBase

    return frozenset(FigureCanvasBase.get_supported_filetypes().keys())


def _get_supported_canvas_classes():
    from ._mpl_imports import (
        FigureCanvasAgg,
        FigureCanvasPdf,
        FigureCanvasPS,
        FigureCanvasSVG,
    )

    return frozenset(
        (FigureCanvasAgg, FigureCanvasPdf, FigureCanvasPS, FigureCanvasSVG)
    )


def get_canvas_class(suffix: str) -> Type["FigureCanvasBase"]:
    s = normalize_extension_string(suffix)
    if s not in _get_supported_image_file_formats():
        raise ValueError(f"Unsupported file format '{suffix}'.")
    for cls in _get_supported_canvas_classes():
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

    if psuffix in _get_supported_image_file_formats():
        if suffix in _get_supported_image_file_formats() and suffix != psuffix:
            warnings.warn(
                f"Received two valid image formats {psuffix!r} (from filename) "
                f"and {suffix!r} (from suffix). The former is ignored."
            )
            return f"{name}.{suffix}"
        return str(filename)

    if suffix is None:
        suffix = "png"

    if suffix not in _get_supported_image_file_formats():
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


def invalidate_plot(f):
    @wraps(f)
    def newfunc(self, *args, **kwargs):
        retv = f(self, *args, **kwargs)
        self._plot_valid = False
        return retv

    return newfunc


def invalidate_data(f):
    @wraps(f)
    def newfunc(self, *args, **kwargs):
        retv = f(self, *args, **kwargs)
        self._data_valid = False
        self._plot_valid = False
        return retv

    return newfunc


def invalidate_figure(f):
    @wraps(f)
    def newfunc(self, *args, **kwargs):
        retv = f(self, *args, **kwargs)
        for field in self.plots.keys():
            self.plots[field].figure = None
            self.plots[field].axes = None
            self.plots[field].cax = None
        self._setup_plots()
        return retv

    return newfunc


def validate_plot(f):
    @wraps(f)
    def newfunc(self, *args, **kwargs):
        # TODO: _profile_valid and _data_valid seem to play very similar roles,
        # there's probably room to abstract these into a common operation
        if hasattr(self, "_data_valid") and not self._data_valid:
            self._recreate_frb()
        if hasattr(self, "_profile_valid") and not self._profile_valid:
            self._recreate_profile()
        if not self._plot_valid:
            # it is the responsibility of _setup_plots to
            # call plot.run_callbacks()
            self._setup_plots()
        retv = f(self, *args, **kwargs)
        return retv

    return newfunc


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

    if len(args) % 2 != 0:
        raise TypeError("Number of arguments must be even.")
    n_pairs = len(args) // 2
    new_args = []
    for i in range(n_pairs):
        x_id = i * 2
        new_args.append(args[x_id + 1])
        new_args.append(args[x_id])
    return tuple(new_args)
