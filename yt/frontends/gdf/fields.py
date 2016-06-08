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

from yt.fields.field_info_container import \
    FieldInfoContainer

# The nice thing about GDF is that for the most part, everything is in CGS,
# with potentially a scalar modification.

class GDFFieldInfo(FieldInfoContainer):
    known_other_fields = (
        ("density", ("g/cm**3", ["density"], None)),
        ("specific_energy", ("erg/g", ["thermal_energy"], None)),
        ("pressure", ("erg/cm**3", ["pressure"], None)),
        ("temperature", ("K", ["temperature"], None)),
        ("velocity_x", ("cm/s", ["velocity_x"], None)),
        ("velocity_y", ("cm/s", ["velocity_y"], None)),
        ("velocity_z", ("cm/s", ["velocity_z"], None)),
        ("mag_field_x", ("gauss", ["magnetic_field_x"], None)),
        ("mag_field_y", ("gauss", ["magnetic_field_y"], None)),
        ("mag_field_z", ("gauss", ["magnetic_field_z"], None)),
    )
    known_particle_fields = ()

    def setup_fluid_fields(self):
        from yt.fields.magnetic_field import \
            setup_magnetic_field_aliases
        setup_magnetic_field_aliases(self, "gdf", ["magnetic_field_%s" % ax for ax in "xyz"])
