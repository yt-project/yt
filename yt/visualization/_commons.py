import os
import sys
import warnings
from functools import wraps
from typing import TYPE_CHECKING, Optional, Type, TypeVar

if sys.version_info >= (3, 8):
    from importlib.metadata import version
else:
    from importlib_metadata import version

import numpy as np
from more_itertools import always_iterable
from packaging.version import Version

from yt.config import ytcfg

if sys.version_info >= (3, 10):
    pass
else:
    from yt._maintenance.backports import zip

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


def get_log_minorticks(vmin: float, vmax: float) -> np.ndarray:
    """calculate positions of linear minorticks on a log colorbar

    Parameters
    ----------
    vmin : float
        the minimum value in the colorbar
    vmax : float
        the maximum value in the colorbar

    """
    expA = np.floor(np.log10(vmin))
    expB = np.floor(np.log10(vmax))
    cofA = np.ceil(vmin / 10**expA).astype("int64")
    cofB = np.floor(vmax / 10**expB).astype("int64")
    lmticks = np.empty(0)
    while cofA * 10**expA <= cofB * 10**expB:
        if expA < expB:
            lmticks = np.hstack((lmticks, np.linspace(cofA, 9, 10 - cofA) * 10**expA))
            cofA = 1
            expA += 1
        else:
            lmticks = np.hstack(
                (lmticks, np.linspace(cofA, cofB, cofB - cofA + 1) * 10**expA)
            )
            expA += 1
    return np.array(lmticks)


def get_symlog_minorticks(linthresh: float, vmin: float, vmax: float) -> np.ndarray:
    """calculate positions of linear minorticks on a symmetric log colorbar

    Parameters
    ----------
    linthresh : float
        the threshold for the linear region
    vmin : float
        the minimum value in the colorbar
    vmax : float
        the maximum value in the colorbar

    """
    if vmin > 0:
        return get_log_minorticks(vmin, vmax)
    elif vmax < 0 and vmin < 0:
        return -get_log_minorticks(-vmax, -vmin)
    elif vmin == 0:
        return np.hstack((0, get_log_minorticks(linthresh, vmax)))
    elif vmax == 0:
        return np.hstack((-get_log_minorticks(linthresh, -vmin)[::-1], 0))
    else:
        return np.hstack(
            (
                -get_log_minorticks(linthresh, -vmin)[::-1],
                0,
                get_log_minorticks(linthresh, vmax),
            )
        )


def get_symlog_majorticks(linthresh: float, vmin: float, vmax: float) -> np.ndarray:
    """calculate positions of major ticks on a log colorbar

    Parameters
    ----------
    linthresh : float
        the threshold for the linear region
    vmin : float
        the minimum value in the colorbar
    vmax : float
        the maximum value in the colorbar

    """
    if vmin >= 0.0:
        yticks = [vmin] + list(
            10
            ** np.arange(
                np.rint(np.log10(linthresh)),
                np.ceil(np.log10(1.1 * vmax)),
            )
        )
    elif vmax <= 0.0:
        if MPL_VERSION >= Version("3.5.0b"):
            offset = 0
        else:
            offset = 1

        yticks = list(
            -(
                10
                ** np.arange(
                    np.floor(np.log10(-vmin)),
                    np.rint(np.log10(linthresh)) - offset,
                    -1,
                )
            )
        ) + [vmax]
    else:
        yticks = (
            list(
                -(
                    10
                    ** np.arange(
                        np.floor(np.log10(-vmin)),
                        np.rint(np.log10(linthresh)) - 1,
                        -1,
                    )
                )
            )
            + [0]
            + list(
                10
                ** np.arange(
                    np.rint(np.log10(linthresh)),
                    np.ceil(np.log10(1.1 * vmax)),
                )
            )
        )
    if yticks[-1] > vmax:
        yticks.pop()
    return np.array(yticks)


def get_default_from_config(data_source, *, field, keys, defaults):
    _keys = list(always_iterable(keys))
    _defaults = list(always_iterable(defaults))

    ftype, fname = data_source._determine_fields(field)[0]
    ret = [
        ytcfg.get_most_specific("plot", ftype, fname, key, fallback=default)
        for key, default in zip(_keys, _defaults, strict=True)
    ]
    if len(ret) == 1:
        return ret[0]
    else:
        return ret
