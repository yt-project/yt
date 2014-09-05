"""
Cosmology related fields.




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.fields.field_plugin_registry import \
    register_field_plugin

@register_field_plugin
def setup_halo_analysis_fields(registry, ftype = "gas", slice_info = None):
    # slice_info would be the left, the right, and the factor.
    # For example, with the old Enzo-ZEUS fields, this would be:
    # slice(None, -2, None)
    # slice(1, -1, None)
    # 1.0
    # Otherwise, we default to a centered difference.
    if slice_info is None:
        sl_left = slice(None, -2, None)
        sl_right = slice(2, None, None)
        div_fac = 2.0
    else:
        sl_left, sl_right, div_fac = slice_info

    def _virial_radius(field, data):
        virial_radius = data.get_field_parameter("virial_radius")
        return data["radius"] / virial_radius

    registry.add_field(("index", "virial_radius"),
                       function=_virial_radius,
                       units="")
