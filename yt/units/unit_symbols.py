from unyt.unit_object import Unit
from unyt.unit_systems import add_symbols

from yt.units.unit_registry import default_unit_registry

add_symbols(globals(), registry=default_unit_registry)


class _SymbolContainer:
    """A container for units to associate with a dataset.

    This object is usually accessed on a Dataset instance via
    ``ds.units.unit_symbols``.

    Parameters
    ----------
    registry : UnitRegistry instance
        A unit registry to associate with units accessed on this object.

    Example
    -------

    >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
    >>> code_mass = ds.units.code_mass
    >>> (12 * code_mass).to("Msun")
    unyt_quantity(4.89719136e+11, 'Msun')
    >>> code_mass.registry is ds.unit_registry
    True
    """

    def __init__(self, registry):
        self._registry = registry
        self._cache = {}

    def __dir__(self):
        ret = [u for u in globals() if not u.startswith("_")]
        ret += list(self._registry.keys())
        ret += object.__dir__(self)
        return list(set(ret))

    def __getattr__(self, item):
        if item in self._cache:
            return self._cache[item]
        if hasattr(globals(), item):
            ret = Unit(globals()[item].expr, registry=self._registry)
        elif item in self._registry:
            ret = Unit(item, registry=self._registry)
        else:
            raise AttributeError(item)
        self._cache[item] = ret
        return ret
