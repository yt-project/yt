import os
import sys
import warnings
from functools import wraps
from typing import TYPE_CHECKING, Optional, TypeVar

import matplotlib as mpl
from matplotlib.ticker import SymmetricalLogLocator
from more_itertools import always_iterable

from yt.config import ytcfg

if sys.version_info >= (3, 10):
    pass
else:
    from yt._maintenance.backports import zip

if TYPE_CHECKING:
    from matplotlib.backend_bases import FigureCanvasBase


_DEFAULT_FONT_PROPERTIES = None


def get_default_font_properties():
    global _DEFAULT_FONT_PROPERTIES
    if _DEFAULT_FONT_PROPERTIES is None:
        import importlib.resources as importlib_resources

        _yt_style = mpl.rc_params_from_file(
            importlib_resources.files("yt") / "default.mplstyle",
            use_default_template=False,
        )
        _DEFAULT_FONT_PROPERTIES = {
            "family": _yt_style["font.family"][0],
            "math_fontfamily": _yt_style["mathtext.fontset"],
        }

    return _DEFAULT_FONT_PROPERTIES


def _get_supported_image_file_formats():
    from matplotlib.backend_bases import FigureCanvasBase

    return frozenset(FigureCanvasBase.get_supported_filetypes().keys())


def _get_supported_canvas_classes():
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.backends.backend_pdf import FigureCanvasPdf
    from matplotlib.backends.backend_ps import FigureCanvasPS
    from matplotlib.backends.backend_svg import FigureCanvasSVG

    return frozenset(
        (FigureCanvasAgg, FigureCanvasPdf, FigureCanvasPS, FigureCanvasSVG)
    )


def get_canvas_class(suffix: str) -> type["FigureCanvasBase"]:
    s = suffix.removeprefix(".")
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


def validate_image_name(filename, suffix: Optional[str] = None) -> str:
    """
    Build a valid image filename with a specified extension (default to png).
    The suffix parameter is ignored if the input filename has a valid extension already.
    Otherwise, suffix is appended to the filename, replacing any existing extension.
    """
    name, psuffix = os.path.splitext(filename)
    psuffix = psuffix.removeprefix(".")

    if suffix is not None:
        suffix = suffix.removeprefix(".")

    if psuffix in _get_supported_image_file_formats():
        if suffix in _get_supported_image_file_formats() and suffix != psuffix:
            warnings.warn(
                f"Received two valid image formats {psuffix!r} (from filename) "
                f"and {suffix!r} (from suffix). The former is ignored.",
                stacklevel=2,
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


class _MPL38_SymmetricalLogLocator(SymmetricalLogLocator):
    # Backporting behaviour from matplotlib 3.8 (in development at the time of writing)
    # see https://github.com/matplotlib/matplotlib/pull/25970

    def __init__(self, *args, **kwargs):
        if mpl.__version_info__ >= (3, 8):
            raise RuntimeError(
                "_MPL38_SymmetricalLogLocator is not needed with matplotlib>=3.8"
            )
        super().__init__(*args, **kwargs)

    def tick_values(self, vmin, vmax):
        linthresh = self._linthresh
        if vmax < vmin:
            vmin, vmax = vmax, vmin
        if -linthresh <= vmin < vmax <= linthresh:
            # only the linear range is present
            return sorted({vmin, 0, vmax})

        return super().tick_values(vmin, vmax)


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


def _get_units_label(units: str) -> str:
    if r"\frac" in units:
        return r"$\ \ \left(%s\right)$" % units
    elif units:
        return r"$\ \ (%s)$" % units
    else:
        return ""
