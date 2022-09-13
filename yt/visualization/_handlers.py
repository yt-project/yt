import sys
import weakref
from numbers import Real
from typing import Any, Dict, List, Optional, Tuple, Type, Union

import numpy as np
import unyt as un
from matplotlib.colors import Colormap, LogNorm, Normalize, SymLogNorm
from packaging.version import Version
from unyt import unyt_quantity

from yt._typing import Quantity, Unit
from yt.config import ytcfg
from yt.funcs import get_brewer_cmap, is_sequence, mylog
from yt.visualization._commons import MPL_VERSION
from yt.visualization.color_maps import _get_cmap

if sys.version_info >= (3, 8):
    from typing import Literal
else:
    from typing_extensions import Literal


class NormHandler:
    """
    A bookkeeper class that can hold a fully defined norm object, or dynamically
    build one on demand according to a set of constraints.

    If a fully defined norm object is added, any existing constraints are
    dropped, and vice versa. These rules are implemented with properties and
    watcher patterns.

    It also keeps track of display units so that vmin, vmax and linthresh can be
    updated with implicit units.
    """

    # using slots here to minimize the risk of introducing bugs
    # since attributes names are essential to this class's implementation
    __slots__ = (
        "data_source",
        "ds",
        "_display_units",
        "_vmin",
        "_vmax",
        "_dynamic_range",
        "_norm_type",
        "_linthresh",
        "_norm_type",
        "_norm",
    )
    _constraint_attrs: List[str] = [
        "vmin",
        "vmax",
        "dynamic_range",
        "norm_type",
        "linthresh",
    ]

    def __init__(
        self,
        data_source,
        *,
        display_units: un.Unit,
        vmin: Optional[un.unyt_quantity] = None,
        vmax: Optional[un.unyt_quantity] = None,
        dynamic_range: Optional[float] = None,
        norm_type: Optional[Type[Normalize]] = None,
        norm: Optional[Normalize] = None,
        linthresh: Optional[float] = None,
    ):
        self.data_source = weakref.proxy(data_source)
        self.ds = data_source.ds  # should already be a weakref proxy
        self._display_units = display_units

        self._norm = norm
        self._vmin = vmin
        self._vmax = vmax
        self._dynamic_range = dynamic_range
        self._norm_type = norm_type
        self._linthresh = linthresh

        if self.has_norm and self.has_constraints:
            raise TypeError(
                "NormHandler input is malformed. "
                "A norm cannot be passed along other constraints."
            )

    def _get_constraints(self) -> Dict[str, Any]:
        return {
            attr: getattr(self, attr)
            for attr in self.__class__._constraint_attrs
            if getattr(self, attr) is not None
        }

    @property
    def has_constraints(self) -> bool:
        return bool(self._get_constraints())

    def _reset_constraints(self) -> None:
        constraints = self._get_constraints()
        if not constraints:
            return

        msg = ", ".join([f"{name}={value}" for name, value in constraints.items()])
        mylog.warning("Dropping norm constraints (%s)", msg)
        for name in constraints.keys():
            setattr(self, name, None)

    @property
    def has_norm(self) -> bool:
        return self._norm is not None

    def _reset_norm(self):
        if not self.has_norm:
            return
        mylog.warning("Dropping norm (%s)", self.norm)
        self._norm = None

    def to_float(self, val: un.unyt_quantity) -> float:
        return float(val.to(self.display_units).d)

    def to_quan(self, val) -> un.unyt_quantity:
        if isinstance(val, un.unyt_quantity):
            return self.ds.quan(val)
        elif (
            is_sequence(val)
            and len(val) == 2
            and isinstance(val[0], Real)
            and isinstance(val[1], (str, un.Unit))
        ):
            return self.ds.quan(*val)
        elif isinstance(val, Real):
            return self.ds.quan(val, self.display_units)
        else:
            raise TypeError(f"Could not convert {val!r} to unyt_quantity")

    @property
    def display_units(self) -> un.Unit:
        return self._display_units

    @display_units.setter
    def display_units(self, newval: Unit) -> None:
        self._display_units = un.Unit(newval, registry=self.ds.unit_registry)

    def _set_quan_attr(
        self, attr: str, newval: Optional[Union[Quantity, float]]
    ) -> None:
        if newval is None:
            setattr(self, attr, None)
        else:
            try:
                quan = self.to_quan(newval)
            except TypeError as exc:
                raise TypeError(
                    "Expected None, a float, or a unyt_quantity, "
                    f"received {newval} with type {type(newval)}"
                ) from exc
            else:
                setattr(self, attr, quan)

    @property
    def vmin(self) -> Optional[Union[un.unyt_quantity, Literal["min"]]]:
        return self._vmin

    @vmin.setter
    def vmin(self, newval: Optional[Union[Quantity, float, Literal["min"]]]) -> None:
        self._reset_norm()
        if newval == "min":
            self._vmin = "min"
        else:
            self._set_quan_attr("_vmin", newval)

    @property
    def vmax(self) -> Optional[Union[un.unyt_quantity, Literal["max"]]]:
        return self._vmax

    @vmax.setter
    def vmax(self, newval: Optional[Union[Quantity, float, Literal["max"]]]) -> None:
        self._reset_norm()
        if newval == "max":
            self._vmax = "max"
        else:
            self._set_quan_attr("_vmax", newval)

    @property
    def dynamic_range(self) -> Optional[float]:
        return self._dynamic_range

    @dynamic_range.setter
    def dynamic_range(self, newval: Optional[float]) -> None:
        if newval is None:
            return

        try:
            newval = float(newval)
        except TypeError:
            raise TypeError(
                f"Expected a float. Received {newval} with type {type(newval)}"
            ) from None

        if newval == 0:
            raise ValueError("Dynamic range cannot be zero.")

        if newval == 1:
            raise ValueError("Dynamic range cannot be unity.")

        self._reset_norm()
        self._dynamic_range = newval

    def get_dynamic_range(
        self, dvmin: Optional[float], dvmax: Optional[float]
    ) -> Tuple[float, float]:
        if self.dynamic_range is None:
            raise RuntimeError(
                "Something went terribly wrong in setting up a dynamic range"
            )

        if self.vmax is None:
            if self.vmin is None:
                raise TypeError(
                    "Cannot set dynamic range with neither "
                    "vmin and vmax being constrained."
                )
            if dvmin is None:
                raise RuntimeError(
                    "Something went terribly wrong in setting up a dynamic range"
                )
            return dvmin, dvmin * self.dynamic_range
        elif self.vmin is None:
            if dvmax is None:
                raise RuntimeError(
                    "Something went terribly wrong in setting up a dynamic range"
                )
            return dvmax / self.dynamic_range, dvmax
        else:
            raise TypeError(
                "Cannot set dynamic range with both "
                "vmin and vmax already constrained."
            )

    @property
    def norm_type(self) -> Optional[Type[Normalize]]:
        return self._norm_type

    @norm_type.setter
    def norm_type(self, newval: Optional[Type[Normalize]]) -> None:
        if not (
            newval is None
            or (isinstance(newval, type) and issubclass(newval, Normalize))
        ):
            raise TypeError(
                "Expected a subclass of matplotlib.colors.Normalize, "
                f"received {newval} with type {type(newval)}"
            )
        self._reset_norm()
        if newval is not SymLogNorm:
            self.linthresh = None
        self._norm_type = newval

    @property
    def norm(self) -> Optional[Normalize]:
        return self._norm

    @norm.setter
    def norm(self, newval: Normalize) -> None:
        if not isinstance(newval, Normalize):
            raise TypeError(
                "Expected a matplotlib.colors.Normalize object, "
                f"received {newval} with type {type(newval)}"
            )
        self._reset_constraints()
        self._norm = newval

    @property
    def linthresh(self) -> Optional[float]:
        return self._linthresh

    @linthresh.setter
    def linthresh(self, newval: Optional[Union[Quantity, float]]) -> None:
        self._reset_norm()
        self._set_quan_attr("_linthresh", newval)
        if self._linthresh is not None and self._linthresh <= 0:
            raise ValueError(
                f"linthresh can only be set to strictly positive values, got {newval}"
            )
        if newval is not None:
            self.norm_type = SymLogNorm

    def get_norm(self, data: np.ndarray, *args, **kw) -> Normalize:
        if self.has_norm:
            return self.norm

        dvmin = dvmax = None

        finite_values_mask = np.isfinite(data)
        # FUTURE: when the minimal supported version of numpy reaches 1.16 or newer,
        # this complicated conditional can be simplified into
        # if self.vmin not in (None, "min"):
        if self.vmin is not None and not (
            isinstance(self.vmin, str) and self.vmin == "min"
        ):
            dvmin = self.to_float(self.vmin)
        elif np.any(finite_values_mask):
            dvmin = self.to_float(np.nanmin(data[finite_values_mask]))

        # FUTURE: see above
        if self.vmax is not None and not (
            isinstance(self.vmax, str) and self.vmax == "max"
        ):
            dvmax = self.to_float(self.vmax)
        elif np.any(finite_values_mask):
            dvmax = self.to_float(np.nanmax(data[finite_values_mask]))

        if self.dynamic_range is not None:
            dvmin, dvmax = self.get_dynamic_range(dvmin, dvmax)

        if dvmin is None:
            dvmin = 1 * getattr(data, "units", 1)
        kw.setdefault("vmin", dvmin)

        if dvmax is None:
            dvmax = 1 * getattr(data, "units", 1)
        kw.setdefault("vmax", dvmax)

        if data.ndim == 3:
            assert data.shape[-1] == 4
            # this is an RGBA array, only linear normalization makes sense here
            norm_type = Normalize
        elif self.norm_type is not None:
            # this is a convenience mechanism for backward compat,
            # allowing to toggle between lin and log scaling without detailed user input
            norm_type = self.norm_type
        else:
            if kw["vmin"] == kw["vmax"] or not np.any(finite_values_mask):
                norm_type = Normalize
            elif kw["vmin"] <= 0:
                # note: see issue 3944 (and PRs and issues linked therein) for a
                # discussion on when to switch to SymLog and related questions
                # of how to calculate a default linthresh value.
                norm_type = SymLogNorm
            else:
                norm_type = LogNorm

        if norm_type is SymLogNorm:
            if self.linthresh is not None:
                linthresh = self.to_float(self.linthresh)
            else:
                linthresh = self._guess_linthresh(data[finite_values_mask])

            kw.setdefault("linthresh", linthresh)
            if MPL_VERSION >= Version("3.2"):
                # note that this creates an inconsistency between mpl versions
                # since the default value previous to mpl 3.4.0 is np.e
                # but it is only exposed since 3.2.0
                kw.setdefault("base", 10)

        return norm_type(*args, **kw)

    def _guess_linthresh(self, finite_plot_data):
        # finite_plot_data is the ImageArray or ColorbarHandler data, already
        # filtered to be finite values

        # get the extrema for the negative and positive values separately
        # neg_min -> neg_max -> 0 -> pos_min -> pos_max
        def get_minmax(data):
            if len(data) > 0:
                return np.nanmin(data), np.nanmax(data)
            return None, None

        pos_min, pos_max = get_minmax(finite_plot_data[finite_plot_data > 0])
        neg_min, neg_max = get_minmax(finite_plot_data[finite_plot_data < 0])

        has_pos = pos_min is not None
        has_neg = neg_min is not None

        # the starting guess is the absolute value of the point closest to 0
        # (remember: neg_max is closer to 0 than neg_min)
        if has_pos and has_neg:
            linthresh = np.min((-neg_max, pos_min))
        elif has_pos:
            linthresh = pos_min
        elif has_neg:
            linthresh = -neg_max
        else:
            # this condition should be handled before here in get_norm
            raise RuntimeError("No finite data points.")

        log10_linthresh = np.log10(linthresh)

        # if either the pos or neg ranges exceed cutoff_sigdigs, then
        # linthresh is shifted to decrease the range to avoid floating point
        # precision errors in the normalization.
        cutoff_sigdigs = 15  # max allowable range in significant digits
        if has_pos and np.log10(pos_max) - log10_linthresh > cutoff_sigdigs:
            linthresh = pos_max / (10.0**cutoff_sigdigs)
            log10_linthresh = np.log10(linthresh)

        if has_neg and np.log10(-neg_min) - log10_linthresh > cutoff_sigdigs:
            linthresh = np.abs(neg_min) / (10.0**cutoff_sigdigs)

        if isinstance(linthresh, unyt_quantity):
            # if the original plot_data has units, linthresh will have units here
            return self.to_float(linthresh)

        return linthresh


class ColorbarHandler:
    __slots__ = ("_draw_cbar", "_draw_minorticks", "_cmap", "_background_color")

    def __init__(
        self,
        *,
        draw_cbar: bool = True,
        draw_minorticks: bool = True,
        cmap: Optional[Union[Colormap, str]] = None,
        background_color: Optional[str] = None,
    ):
        self._draw_cbar = draw_cbar
        self._draw_minorticks = draw_minorticks
        self._cmap: Optional[Colormap] = None
        self.cmap = cmap
        self._background_color = background_color

    @property
    def draw_cbar(self) -> bool:
        return self._draw_cbar

    @draw_cbar.setter
    def draw_cbar(self, newval) -> None:
        if not isinstance(newval, bool):
            raise TypeError(
                f"Excpected a boolean, got {newval} with type {type(newval)}"
            )
        self._draw_cbar = newval

    @property
    def draw_minorticks(self) -> bool:
        return self._draw_minorticks

    @draw_minorticks.setter
    def draw_minorticks(self, newval) -> None:
        if not isinstance(newval, bool):
            raise TypeError(
                f"Excpected a boolean, got {newval} with type {type(newval)}"
            )
        self._draw_minorticks = newval

    @property
    def cmap(self) -> Colormap:
        return self._cmap or _get_cmap(ytcfg.get("yt", "default_colormap"))

    @cmap.setter
    def cmap(self, newval) -> None:
        if isinstance(newval, Colormap) or newval is None:
            self._cmap = newval
        elif isinstance(newval, str):
            self._cmap = _get_cmap(newval)
        elif is_sequence(newval):
            # tuple colormaps are from palettable (or brewer2mpl)
            self._cmap = get_brewer_cmap(newval)
        else:
            raise TypeError(
                "Expected a colormap object or name, "
                f"got {newval} with type {type(newval)}"
            )

    @property
    def background_color(self) -> Any:
        return self._background_color or "white"

    @background_color.setter
    def background_color(self, newval: Any):
        # not attempting to constrain types here because
        # down the line it really depends on matplotlib.axes.Axes.set_faceolor
        # which is very type-flexibile
        if newval is None:
            self._background_color = self.cmap(0)
        else:
            self._background_color = newval

    @property
    def has_background_color(self) -> bool:
        return self._background_color is not None
