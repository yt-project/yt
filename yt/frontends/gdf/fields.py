"""
GDF-specific fields



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

# The nice thing about GDF is that for the most part, everything is in CGS,
# with potentially a scalar modification.

class GDFFieldInfo(FieldInfoContainer):
    known_other_fields = (
        ("density", ("g/cm**3", [], None)),
        ("specific_energy", ("erg/g", ["thermal_energy"], None)),
        ("pressure", ("erg/cm**3", [], None)),
        ("temperature", ("K", [], None)),
        ("velocity_x", ("cm/s", [], None)),
        ("velocity_y", ("cm/s", [], None)),
        ("velocity_z", ("cm/s", [], None)),
        ("mag_field_x", ("gauss", ["magnetic_field_x"], None)),
        ("mag_field_y", ("gauss", ["magnetic_field_y"], None)),
        ("mag_field_z", ("gauss", ["magnetic_field_z"], None)),
    )
    known_particle_fields = ()
