"""
HaloCatalog-specific fields




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.funcs import mylog
from yt.fields.field_info_container import \
    FieldInfoContainer
from yt.units.yt_array import \
    YTArray

from yt.utilities.physical_constants import \
    mh, \
    mass_sun_cgs

m_units = "g"
p_units = "cm"
v_units = "cm / s"
r_units = "cm"

class HaloCatalogFieldInfo(FieldInfoContainer):
    known_other_fields = (
    )

    known_particle_fields = (
        ("particle_identifier", ("", [], None)),
        ("particle_position_x", (p_units, [], None)),
        ("particle_position_y", (p_units, [], None)),
        ("particle_position_z", (p_units, [], None)),
        ("particle_velocity_x", (v_units, [], None)),
        ("particle_velocity_y", (v_units, [], None)),
        ("particle_velocity_z", (v_units, [], None)),
        ("particle_mass", (m_units, [], "Virial Mass")),
        ("virial_radius", (r_units, [], "Virial Radius")),
)
