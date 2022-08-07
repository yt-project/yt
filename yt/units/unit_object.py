# re-export all of unyt.unit_object's public interface, except for `define_unit`
from unyt.unit_object import sympy_one, Unit, em_conversions, em_conversion_dims, NULL_UNIT


def define_unit(
    symbol, value, tex_repr=None, offset=None, prefixable=False, registry=None
):
    """
    yt.units.define_unit is a thin wrapper around unyt.define_unit
    The key difference is that yt.units.define_unit will register new units to
    yt's default unit registry instead of unyt's
    Units defined with this function are then available to import from the yt.units
    namespace.

    See unyt.define_unit for more details.
    """
    from unyt import define_unit as unyt_define_unit
    import yt.units

    if registry is None:
        registry = yt.units.default_unit_registry

    unyt_define_unit(
        symbol,
        value,
        tex_repr=tex_repr,
        offset=offset,
        prefixable=prefixable,
        registry=registry,
    )

    if registry is yt.units.default_unit_registry:
        u = Unit(symbol, registry=registry)
        setattr(yt.units, symbol, u)
