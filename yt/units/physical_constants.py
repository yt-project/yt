from unyt.array import unyt_quantity
from unyt.unit_systems import add_constants

from yt.units.unit_registry import default_unit_registry

add_constants(globals(), registry=default_unit_registry)


class _ConstantContainer:
    """A container for physical constants to associate with a dataset.

    This object is usually accessed on a Dataset instance via
    ``ds.units.physical_constants``.

    Parameters
    ----------
    registry : UnitRegistry instance
        A unit registry to associate with units constants accessed on
        this object.

    Example
    -------

    >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
    >>> ds.units.physical_constants.newtons_constant
    unyt_quantity(6.67384e-08, 'cm**3/(g*s**2)')
    """

    def __init__(self, registry):
        self._registry = registry
        self._cache = {}

    def __dir__(self):
        ret = [p for p in globals() if not p.startswith("_")] + object.__dir__(self)
        return list(set(ret))

    def __getattr__(self, item):
        if item in self._cache:
            return self._cache[item]
        if item in globals():
            const = globals()[item].copy()
            const.units.registry = self._registry
            const.convert_to_base(self._registry.unit_system)
            const_v, const_unit = const.v, const.units
            ret = unyt_quantity(const_v, const_unit, registry=self._registry)
            self._cache[item] = ret
            return ret
        raise AttributeError(item)
