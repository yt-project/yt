"""
Cholla-specific fields



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2018, yt Development Team.
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
        return data["cholla", "momentum_%s" % comp]/data["cholla","density"]
    return _velocity

class ChollaFieldInfo(FieldInfoContainer):
    known_other_fields = (
        ("density", ("code_mass/code_length**3", ["density"], None)),
        ("cell_centered_B_x", (b_units, [], None)),
        ("cell_centered_B_y", (b_units, [], None)),
        ("cell_centered_B_z", (b_units, [], None)),
        ("gravitational_potential", ("code_velocity**2", ["gravitational_potential"], None)),
    )

# In Cholla, conservative or primitive variables may be written out.
# By default, yt concerns itself with primitive variables. The following
# field definitions allow for conversions to primitive variables in the
# case that the file contains the conservative ones.

    def setup_fluid_fields(self):
        from yt.fields.magnetic_field import \
            setup_magnetic_field_aliases
        unit_system = self.ds.unit_system
        # Add velocity fields
        for comp in "xyz":
            vel_field = ("cholla", "velocity_%s" % comp)
            mom_field = ("cholla", "momentum_%s" % comp)
            if vel_field in self.field_list:
                self.add_output_field(vel_field, sampling_type="cell",  units="code_length/code_time")
                self.alias(("gas","velocity_%s" % comp), vel_field,
                           units=unit_system["velocity"])
            elif mom_field in self.field_list:
                self.add_output_field(mom_field, sampling_type="cell",
                                      units="code_mass/code_time/code_length**2")
                self.add_field(("gas","velocity_%s" % comp), sampling_type="cell",
                               function=velocity_field(comp), units = unit_system["velocity"])
        # Add pressure, energy, and temperature fields
        def eint_from_etot(data):
            eint = data["cholla","total_energy"].copy()
            eint -= data["gas","kinetic_energy"]
            if ("cholla","cell_centered_B_x") in self.field_list:
                eint -= data["gas","magnetic_energy"]
            return eint
        def etot_from_pres(data):
            etot = data["cholla","pressure"].copy()
            etot /= data.ds.gamma-1.
            etot += data["gas", "kinetic_energy"]
            if ("cholla","cell_centered_B_x") in self.field_list:
                etot += data["gas", "magnetic_energy"]
            return etot
        if ("cholla","pressure") in self.field_list:
            self.add_output_field(("cholla","pressure"), sampling_type="cell",
                                  units=pres_units)
            self.alias(("gas","pressure"),("cholla","pressure"),
                       units=unit_system["pressure"])
            def _thermal_energy(field, data):
                return data["cholla","pressure"] / \
                       (data.ds.gamma-1.)/data["cholla","density"]
            self.add_field(("gas","thermal_energy"), sampling_type="cell",
                           function=_thermal_energy,
                           units=unit_system["specific_energy"])
            def _total_energy(field, data):
                return etot_from_pres(data)/data["cholla","density"]
            self.add_field(("gas","total_energy"), sampling_type="cell",
                           function=_total_energy,
                           units=unit_system["specific_energy"])
        elif ("cholla","total_energy") in self.field_list:
            self.add_output_field(("cholla","total_energy"), sampling_type="cell",
                                  units=pres_units)
            def _thermal_energy(field, data):
                return eint_from_etot(data)/data["cholla","density"]
            self.add_field(("gas","thermal_energy"), sampling_type="cell",
                           function=_thermal_energy,
                           units=unit_system["specific_energy"])
            def _total_energy(field, data):
                return data["cholla","total_energy"]/data["cholla","density"]
            self.add_field(("gas","total_energy"), sampling_type="cell",
                           function=_total_energy,
                           units=unit_system["specific_energy"])
        # Add temperature field
        def _temperature(field, data):
            if data.has_field_parameter("mu"):
                mu = data.get_field_parameter("mu")
            else:
                mu = 0.6
            return mu*mh*data["gas","pressure"]/data["gas","density"]/kboltz
        self.add_field(("gas","temperature"), sampling_type="cell",  function=_temperature,
                       units=unit_system["temperature"])

        setup_magnetic_field_aliases(self, "cholla", ["cell_centered_B_%s" % ax for ax in "xyz"])
