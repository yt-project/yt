"""
Athena-specific fields



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
pres_units = "code_pressure"
erg_units = "code_mass * (code_length/code_time)**2"
rho_units = "code_mass / code_length**3"

def velocity_field(comp):
    def _velocity(field, data):
        return data["athena", "momentum_%s" % comp]/data["athena","density"]
    return _velocity

class AthenaFieldInfo(FieldInfoContainer):
    known_other_fields = (
        ("density", ("code_mass/code_length**3", ["density"], None)),
        ("cell_centered_B_x", (b_units, [], None)),
        ("cell_centered_B_y", (b_units, [], None)),
        ("cell_centered_B_z", (b_units, [], None)),
        ("gravitational_potential", ("code_velocity**2", ["gravitational_potential"], None)),
    )

# In Athena, conservative or primitive variables may be written out.
# By default, yt concerns itself with primitive variables. The following
# field definitions allow for conversions to primitive variables in the
# case that the file contains the conservative ones.

    def setup_fluid_fields(self):
        from yt.fields.magnetic_field import \
            setup_magnetic_field_aliases
        unit_system = self.ds.unit_system
        # Add velocity fields
        for comp in "xyz":
            vel_field = ("athena", "velocity_%s" % comp)
            mom_field = ("athena", "momentum_%s" % comp)
            if vel_field in self.field_list:
                self.add_output_field(vel_field, units="code_length/code_time")
                self.alias(("gas","velocity_%s" % comp), vel_field,
                           units=unit_system["velocity"])
            elif mom_field in self.field_list:
                self.add_output_field(mom_field,
                                      units="code_mass/code_time/code_length**2")
                self.add_field(("gas","velocity_%s" % comp),
                               function=velocity_field(comp), units = unit_system["velocity"])
        # Add pressure, energy, and temperature fields
        def ekin1(data):
            return 0.5*(data["athena","momentum_x"]**2 +
                        data["athena","momentum_y"]**2 +
                        data["athena","momentum_z"]**2)/data["athena","density"]
        def ekin2(data):
            return 0.5*(data["athena","velocity_x"]**2 +
                        data["athena","velocity_y"]**2 +
                        data["athena","velocity_z"]**2)*data["athena","density"]
        def emag(data):
            return 0.5*(data["cell_centered_B_x"]**2 +
                        data["cell_centered_B_y"]**2 +
                        data["cell_centered_B_z"]**2)
        def eint_from_etot(data):
            eint = data["athena","total_energy"]
            eint -= ekin1(data)
            if ("athena","cell_centered_B_x") in self.field_list:
                eint -= emag(data)
            return eint
        def etot_from_pres(data):
            etot = data["athena","pressure"]/(data.ds.gamma-1.)
            etot += ekin2(data)
            if ("athena","cell_centered_B_x") in self.field_list:
                etot += emag(data)
            return etot
        if ("athena","pressure") in self.field_list:
            self.add_output_field(("athena","pressure"),
                                  units=pres_units)
            self.alias(("gas","pressure"),("athena","pressure"),
                       units=unit_system["pressure"])
            def _thermal_energy(field, data):
                return data["athena","pressure"] / \
                       (data.ds.gamma-1.)/data["athena","density"]
            self.add_field(("gas","thermal_energy"),
                           function=_thermal_energy,
                           units=unit_system["specific_energy"])
            def _total_energy(field, data):
                return etot_from_pres(data)/data["athena","density"]
            self.add_field(("gas","total_energy"),
                           function=_total_energy,
                           units=unit_system["specific_energy"])
        elif ("athena","total_energy") in self.field_list:
            self.add_output_field(("athena","total_energy"),
                                  units=pres_units)
            def _pressure(field, data):
                return eint_from_etot(data)*(data.ds.gamma-1.0)
            self.add_field(("gas","pressure"), function=_pressure,
                           units=unit_system["pressure"])
            def _thermal_energy(field, data):
                return eint_from_etot(data)/data["athena","density"]
            self.add_field(("gas","thermal_energy"),
                           function=_thermal_energy,
                           units=unit_system["specific_energy"])
            def _total_energy(field, data):
                return data["athena","total_energy"]/data["athena","density"]
            self.add_field(("gas","total_energy"),
                           function=_total_energy,
                           units=unit_system["specific_energy"])

        def _temperature(field, data):
            if data.has_field_parameter("mu"):
                mu = data.get_field_parameter("mu")
            else:
                mu = 0.6
            return mu*mh*data["gas","pressure"]/data["gas","density"]/kboltz
        self.add_field(("gas","temperature"), function=_temperature,
                       units=unit_system["temperature"])

        setup_magnetic_field_aliases(self, "athena", ["cell_centered_B_%s" % ax for ax in "xyz"])


