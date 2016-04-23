"""
Athena++-specific fields



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
from yt.utilities.physical_constants import \
    kboltz, mh

b_units = "code_magnetic"
pres_units = "code_mass/(code_length*code_time**2)"
rho_units = "code_mass / code_length**3"
vel_units = "code_length / code_time"

class AthenaPPFieldInfo(FieldInfoContainer):
    known_other_fields = (
        ("rho", (rho_units, ["density"], None)),
        ("pgas", (pres_units, ["pressure"], None)),
        ("vel1", (vel_units, ["velocity_x"], None)),
        ("vel2", (vel_units, ["velocity_y"], None)),
        ("vel3", (vel_units, ["velocity_z"], None)),
        ("B1", (b_units, [], None)),
        ("B2", (b_units, [], None)),
        ("B3", (b_units, [], None)),
    )

    def setup_fluid_fields(self):
        from yt.fields.magnetic_field import \
            setup_magnetic_field_aliases
        unit_system = self.ds.unit_system
        def _temperature(field, data):
            if data.has_field_parameter("mu"):
                mu = data.get_field_parameter("mu")
            else:
                mu = 0.6
            return mu*mh*data["gas","pressure"]/data["gas","density"]/kboltz
        self.add_field(("gas","temperature"), function=_temperature,
                       units=unit_system["temperature"])

        setup_magnetic_field_aliases(self, "athena++", ["B%d" % ax for ax in (1,2,3)])


