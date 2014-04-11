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
    mh, boltzmann_constant_cgs, amu_cgs
from yt.fields.field_info_container import \
    FieldInfoContainer

rho_units = "code_mass / code_length**3"
mom_units = "code_mass / (code_time * code_length**2)"
eden_units = "code_mass / (code_time**2 * code_length)" # erg / cm^3

def _thermal_energy_density(field, data):
    # What we've got here is UEINT:
    # u here is velocity
    # E is energy density from the file
    #   rho e = rho E - rho * u * u / 2
    ke = 0.5 * ( data["momentum_x"]**2
               + data["momentum_y"]**2
               + data["momentum_z"]**2) / data["density"]
    return data["eden"] - ke

def _thermal_energy(field, data):
    # This is little e, so we take thermal_energy_density and divide by density
    return data["thermal_energy_density"] / data["density"]

def _temperature(field,data):
    mu = data.pf.parameters["mu"]
    gamma = data.pf.parameters["gamma"]
    tr  = data["thermal_energy_density"] / data["density"]
    tr *= mu * amu_cgs / boltzmann_constant_cgs
    tr *= (gamma - 1.0)
    return tr

class BoxlibFieldInfo(FieldInfoContainer):
    known_other_fields = (
        ("density", (rho_units, ["density"], None)),
        ("eden", (eden_units, ["energy_density"], None)),
        ("xmom", (mom_units, ["momentum_x"], None)),
        ("ymom", (mom_units, ["momentum_y"], None)),
        ("zmom", (mom_units, ["momentum_z"], None)),
        ("temperature", ("K", ["temperature"], None)),
        ("Temp", ("K", ["temperature"], None)),
        ("x_velocity", ("cm/s", ["velocity_x"], None)),
        ("y_velocity", ("cm/s", ["velocity_y"], None)),
        ("z_velocity", ("cm/s", ["velocity_z"], None)),
        ("xvel", ("cm/s", ["velocity_x"], None)),
        ("yvel", ("cm/s", ["velocity_y"], None)),
        ("zvel", ("cm/s", ["velocity_z"], None)),
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
        # Now, let's figure out what fields are included.
        if any(f[1] == "xmom" for f in self.field_list):
            self.setup_momentum_to_velocity()
        self.add_field(("gas", "thermal_energy"),
                       function = _thermal_energy,
                       units = "erg/g")
        self.add_field(("gas", "thermal_energy_density"),
                       function = _thermal_energy_density,
                       units = "erg/cm**3")
        if ("gas", "temperature") not in self.field_aliases:
            self.add_field(("gas", "temperature"),
                           function=_temperature,
                           units="K")

    def setup_momentum_to_velocity(self):
        def _get_vel(axis):
            def velocity(field, data):
                return data["%smom" % axis]/data["density"]
        for ax in 'xyz':
            self.add_field(("gas", "velocity_%s" % ax),
                           function = _get_vel(ax),
                           units = "cm/s")

class MaestroFieldInfo(BoxlibFieldInfo):

    def __init__(self, *args, **kwargs):
        maestro_fields = (
            ("x_vel", ("cm/s", ["velocity_x"], None)),
            ("y_vel", ("cm/s", ["velocity_y"], None)),
            ("z_vel", ("cm/s", ["velocity_z"], None)),
            ("magvel", ("cm/s", ["velocity_magnitude"], None)),
        )
        self.known_other_fields += maestro_fields
        super(MaestroFieldInfo, self).__init__(*args, **kwargs)
