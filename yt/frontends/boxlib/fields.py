"""
Orion-specific fields



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.utilities.physical_constants import \
    mh, kboltz
from yt.fields.field_info_container import \
    FieldInfoContainer
import yt.fields.universal_fields

rho_units = "code_mass / code_length**3"
mom_units = "code_mass * code_length / code_time"
eden_units = "code_mass / (code_time**2 * code_length)" # erg / cm^3



def _thermal_energy_density(field, data):
    ke = 0.5 * ( data["momentum_x"]**2
               + data["momentum_y"]**2
               + data["momentum_z"]**2) / data["density"]
    return data["eden"] - ke

def _thermal_energy(field, data):
    return data["thermal_energy_density"] / data["density"]

def _temperature(field,data):
    mu = data.get_field_parameter("mu")
    return ( (data.pf.gamma-1.0) * mu * mh *
             data["thermal_energy"] / (kboltz * data["density"]) )


class BoxlibFieldInfo(FieldInfoContainer):
    known_other_fields = (
        ("density", (rho_units, ["density"], None)),
        ("eden", (eden_units, ["energy_density"], None)),
        ("xmom", (mom_units, ["momentum_x"], None)),
        ("ymom", (mom_units, ["momentum_y"], None)),
        ("zmom", (mom_units, ["momentum_z"], None)),
        ("temperature", ("K", ["temperature"], None)),
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
        ("particle_id", ("", ["particle_index"], None)),
        ("particle_mdot", ("code_mass/code_time", [], None)),
        # "mlast",
        # "r",
        # "mdeut",
        # "n",
        # "burnstate",
        # "luminosity",
    )

    def setup_fluid_fields(self):
        def _get_vel(axis):
            def velocity(field, data):
                return data["%smom" % ax]/data["density"]
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
