from unyt.array import (
    loadtxt,
    savetxt,
    unyt_array,
    unyt_quantity,
)
from unyt.unit_object import Unit, define_unit  # NOQA: F401
from unyt.unit_registry import UnitRegistry  # NOQA: Ffg401
from unyt.unit_systems import UnitSystem  # NOQA: F401

from yt.units.physical_constants import *
from yt.units.physical_constants import _ConstantContainer
from yt.units.unit_symbols import *
from yt.units.unit_symbols import _SymbolContainer
from yt.utilities.exceptions import YTArrayTooLargeToDisplay
from yt.units._numpy_wrapper_functions import (
    uconcatenate,
    ucross,
    udot,
    uhstack,
    uintersect1d,
    unorm,
    ustack,
    uunion1d,
    uvstack,
)
YTArray = unyt_array

YTQuantity = unyt_quantity


class UnitContainer:
    """A container for units and constants to associate with a dataset

    This object is usually accessed on a Dataset instance via ``ds.units``.

    Parameters
    ----------
    registry : UnitRegistry instance
        A unit registry to associate with units and constants accessed
        on this object.

    Example
    -------

    >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
    >>> code_mass = ds.units.code_mass
    >>> (12 * code_mass).to("Msun")
    unyt_quantity(4.89719136e+11, 'Msun')
    >>> code_mass.registry is ds.unit_registry
    True
    >>> ds.units.newtons_constant
    unyt_quantity(6.67384e-08, 'cm**3/(g*s**2)')

    """

    def __init__(self, registry):
        self.unit_symbols = _SymbolContainer(registry)
        self.physical_constants = _ConstantContainer(registry)

    def __dir__(self):
        all_dir = self.unit_symbols.__dir__() + self.physical_constants.__dir__()
        all_dir += object.__dir__(self)
        return list(set(all_dir))

    def __getattr__(self, item):
        pc = self.physical_constants
        us = self.unit_symbols
        ret = getattr(us, item, None) or getattr(pc, item, None)
        if not ret:
            raise AttributeError(item)
        return ret


def display_ytarray(arr):
    r"""
    Display a YTArray in a Jupyter widget that enables unit switching.

    The array returned by this function is read-only, and only works with
    arrays of size 3 or lower.

    Parameters
    ----------
    arr : YTArray
        The Array to display; must be of size 3 or lower.

    Examples
    --------
    >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
    >>> display_ytarray(ds.domain_width)
    """
    if arr.size > 3:
        raise YTArrayTooLargeToDisplay(arr.size, 3)
    import ipywidgets

    unit_registry = arr.units.registry
    equiv = unit_registry.list_same_dimensions(arr.units)
    dropdown = ipywidgets.Dropdown(options=sorted(equiv), value=str(arr.units))

    def arr_updater(arr, texts):
        def _value_updater(change):
            arr2 = arr.in_units(change["new"])
            if arr2.shape == ():
                arr2 = [arr2]
            for v, t in zip(arr2, texts):
                t.value = str(v.value)

        return _value_updater

    if arr.shape == ():
        arr_iter = [arr]
    else:
        arr_iter = arr
    texts = [ipywidgets.Text(value=str(_.value), disabled=True) for _ in arr_iter]
    dropdown.observe(arr_updater(arr, texts), names="value")
    return ipywidgets.HBox(texts + [dropdown])


def _wrap_display_ytarray(arr):
    from IPython.core.display import display

    display(display_ytarray(arr))


# monkeypatch __format__ method from unyt 2.9 (which requires Python >= 3.8)
# see https://github.com/yt-project/unyt/pull/188
# We should be able to require unyt >= 2.9 when we drop support for Python 3.7

from packaging.version import Version
from unyt import __version__
if Version(__version__) < Version("2.9"):
    def __mp_unyt_array_format(self, format_spec):
         return "{} {}".format(self.d.__format__(format_spec), self.units)
    unyt_array.__format__ = __mp_unyt_array_format
del __version__
del Version
