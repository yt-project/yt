from yt.units.physical_constants import *
from yt.units.unit_symbols import *
from unyt.array import (
    loadtxt,
    savetxt,
    uconcatenate,
    ucross,
    udot,
    uhstack,
    uintersect1d,
    unorm,
    ustack,
    uunion1d,
    uvstack,
    unyt_array,
    unyt_quantity,
)
from unyt.unit_object import Unit, define_unit  # NOQA: F401
from unyt.unit_registry import UnitRegistry  # NOQA: Ffg401
from unyt.unit_systems import UnitSystem  # NOQA: F401

YTArray = unyt_array

YTQuantity = unyt_quantity

from yt.units.unit_symbols import _SymbolContainer
from yt.units.physical_constants import _ConstantContainer

class UnitContainer(object):
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
    >>> (12*code_mass).to("Msun")
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

