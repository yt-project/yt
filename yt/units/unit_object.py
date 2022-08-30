# re-export all of unyt.unit_object's public interface, except for `define_unit`
from unyt.unit_object import sympy_one, Unit, em_conversions, em_conversion_dims, NULL_UNIT


def define_unit(
    symbol, value, tex_repr=None, offset=None, prefixable=False, registry=None
):
    from unyt import define_unit as unyt_define_unit
    from yt._maintenance.deprecation import issue_deprecation_warning

    issue_deprecation_warning(
        "yt.units.unit_object.define_unit is deprecated. "
        "It is an alias for unyt.define_unit and does not actually make "
        "new units available to yt datasets. "
        "Use unyt.define_unit directly if it fits your needs. "
        "Otherwise, if you wish to add units to yt datasets, use the ds.define_unit "
        "method instead.",
        since="4.1"
    )
    unyt_define_unit(
        symbol,
        value,
        tex_repr=tex_repr,
        offset=offset,
        prefixable=prefixable,
        registry=registry,
    )
