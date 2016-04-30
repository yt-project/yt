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

def velocity_field(j):
    def _velocity(field, data):
        return data["athena++", "mom%d" % j]/data["athena++","dens"]
    return _velocity

class AthenaPPFieldInfo(FieldInfoContainer):
    known_other_fields = (
        ("rho", (rho_units, ["density"], None)),
        ("dens", (rho_units, ["density"], None)),
        ("pgas", (pres_units, ["pressure"], None)),
        ("B1", (b_units, [], None)),
        ("B2", (b_units, [], None)),
        ("B3", (b_units, [], None)),
    )

    def setup_fluid_fields(self):
        from yt.fields.magnetic_field import \
            setup_magnetic_field_aliases
        unit_system = self.ds.unit_system
        # Add velocity fields
        vel_prefix = "velocity"
        if self.ds.geometry != "cartesian":
            vel_prefix += "_"+self.ds.geometry
        for i, comp in enumerate(self.ds.coordinates.axis_order):
            vel_field = ("athena++", "vel%d" % (i+1))
            mom_field = ("athena++", "mom%d" % (i+1))
            if comp == "r":
                comp = "radius"
            if vel_field in self.field_list:
                self.add_output_field(vel_field, units="code_length/code_time")
                self.alias(("gas","%s_%s" % (vel_prefix, comp)), vel_field,
                           units=unit_system["velocity"])
            elif mom_field in self.field_list:
                self.add_output_field(mom_field,
                                      units="code_mass/code_time/code_length**2")
                print(comp)
                self.add_field(("gas","%s_%s" % (vel_prefix, comp)),
                               function=velocity_field(i+1), units=unit_system["velocity"])
        # Add temperature field
        def _temperature(field, data):
            if data.has_field_parameter("mu"):
                mu = data.get_field_parameter("mu")
            else:
                mu = 0.6
            return mu*mh*data["gas","pressure"]/data["gas","density"]/kboltz
        self.add_field(("gas","temperature"), function=_temperature,
                       units=unit_system["temperature"])

        setup_magnetic_field_aliases(self, "athena++", ["B%d" % ax for ax in (1,2,3)])


