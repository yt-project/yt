"""
Chombo-specific fields



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
from yt.fields.field_info_container import \
    FieldInfoContainer
from yt.frontends.boxlib.fields import \
    rho_units, \
    mom_units, \
    eden_units, \
    _thermal_energy_density, \
    _thermal_energy, \
    _temperature

rho_units = "code_mass / code_length**3"
mom_units = "code_mass / (code_time * code_length**2)"
eden_units = "code_mass / (code_time**2 * code_length)" # erg / cm^3

# Chombo does not have any known fields by itself.
class ChomboFieldInfo(FieldInfoContainer):
    known_other_fields = ()
    known_particle_fields = ()

# Orion 2 Fields
# We duplicate everything here from Boxlib, because we want to be able to
# subclass it and that can be somewhat tricky.
class Orion2FieldInfo(ChomboFieldInfo):
    known_other_fields = (
        ("density", (rho_units, ["density"], None)),
        ("energy-density", (eden_units, ["energy_density"], None)),
        ("radiation-energy-density", (eden_units, ["radiation_energy_density"], None)),
        ("X-momentum", (mom_units, ["momentum_x"], None)),
        ("Y-momentum", (mom_units, ["momentum_y"], None)),
        ("Z-momentum", (mom_units, ["momentum_z"], None)),
        ("temperature", ("K", ["temperature"], None)),
        ("X-magnfield", ("gauss", ["magnetic_field_x"], None)),
        ("Y-magnfield", ("gauss", ["magnetic_field_y"], None)),
        ("Z-magnfield", ("gauss", ["magnetic_field_z"], None)),
    )

    known_particle_fields = (
        ("particle_mass", ("code_mass", [], None)),
        ("particle_position_x", ("code_length", [], None)),
        ("particle_position_y", ("code_length", [], None)),
        ("particle_position_z", ("code_length", [], None)),
        ("particle_momentum_x", (mom_units, [], None)),
        ("particle_momentum_y", (mom_units, [], None)),
        ("particle_momentum_z", (mom_units, [], None)),
        # Note that these are *internal* agmomen
        ("particle_angmomen_x", ("code_length**2/code_time", [], None)),
        ("particle_angmomen_y", ("code_length**2/code_time", [], None)),
        ("particle_angmomen_z", ("code_length**2/code_time", [], None)),
        ("particle_mlast", ("code_mass", [], None)),
        ("particle_r", ("code_length", [], None)),
        ("particle_mdeut", ("code_mass", [], None)),
        ("particle_n", ("", [], None)),
        ("particle_mdot", ("code_mass/code_time", [], None)),
        ("particle_burnstate", ("", [], None)),
        ("particle_luminosity", ("", [], None)),
        ("particle_id", ("", ["particle_index"], None)),
    )

    def setup_fluid_fields(self):
        def _get_vel(axis):
            def velocity(field, data):
                return data["momentum_%s" % ax]/data["density"]
            return velocity
        for ax in 'xyz':
            self.add_field("velocity_%s" % ax, function = _get_vel(ax),
                           units = "cm/s")
        self.add_field("thermal_energy",
                       function = _thermal_energy,
                       units = "erg/g")
        self.add_field("thermal_energy_density",
                       function = _thermal_energy_density,
                       units = "erg/cm**3")
        self.add_field("temperature", function=_temperature,
                       units="K")
