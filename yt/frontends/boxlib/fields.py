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
import string

from yt.utilities.physical_constants import \
    mh, boltzmann_constant_cgs, amu_cgs
from yt.fields.field_info_container import \
    FieldInfoContainer
from yt.fields.species_fields import \
    add_species_field_by_fraction
from yt.utilities.chemical_formulas import \
    ChemicalFormula

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

class MaestroFieldInfo(FieldInfoContainer):

    known_other_fields = (
        ("density", ("g/cm**3", ["density"], None)),
        ("x_vel", ("cm/s", ["velocity_x"], None)),
        ("y_vel", ("cm/s", ["velocity_y"], None)),
        ("z_vel", ("cm/s", ["velocity_z"], None)),
        ("magvel", ("cm/s", ["velocity_magnitude"], None)),
        ("tfromp", ("K", [], None)),
        ("tfromh", ("K", [], None)),
        ("Machnumber", ("", ["mach_number"], None)),
        ("S", ("1/s", [], None)),
        ("ad_excess", ("", [], "Adiabatic Excess")),
        ("deltaT", ("", [], None)),
        ("deltagamma", ("", [], None)),
        ("deltap", ("", [], None)),
        ("divw0", ("1/s", [], None)),
        # Specific entropy
        ("entropy", ("erg/(g*K)", ["entropy"], None)),
        ("entropypert", ("", [], None)),
        ("enucdot", ("ergs/(g*s)", [], None)),
        ("gpi_x", ("dyne/cm**3", [], None)), # Perturbational pressure grad
        ("gpi_y", ("dyne/cm**3", [], None)),
        ("gpi_z", ("dyne/cm**3", [], None)),
        ("h", ("erg/g", [], "Specific Enthalpy")),
        ("h0", ("erg/g", [], "Base State Specific Enthalpy")),
        # Momentum cannot be computed because we need to include base and
        # full state.
        ("momentum", ("g*cm/s", ["momentum_magnitude"], None)),
        ("p0", ("erg/cm**3", [], "p_0")),
        ("p0pluspi", ("erg/cm**3", [], "p_0 + \pi")),
        ("pi", ("erg/cm**3", [], None)),
        ("pioverp0", ("", [], "\pi/p_0")),
        # Base state density
        ("rho0", ("g/cm**3", [], "\\rho_0")),
        ("rhoh", ("erg/cm**3", ["enthalpy_density"], "(\\rho h)")),
        # Base state enthalpy density
        ("rhoh0", ("erg/cm**3", [], "(\\rho h)_0")),
        ("rhohpert", ("erg/cm**3", [], "(\\rho h)^\prime")),
        ("rhopert", ("g/cm**3", [], "\\rho^\prime")),
        ("soundspeed", ("cm/s", ["sound_speed"], None)),
        ("sponge", ("", [], None)),
        ("tpert", ("K", [], None)),
        # Again, base state -- so we can't compute ourselves.
        ("vort", ("1/s", ["vorticity_magnitude"], None)),
        # Base state
        ("w0_x", ("cm/s", [], None)),
        ("w0_y", ("cm/s", [], None)),
        ("w0_z", ("cm/s", [], None)),
    )

    def setup_fluid_fields(self):
        # Add omegadots, units of 1/s
        if self.pf.parameters["use_tfromp"]:
            self.alias(("gas", "temperature"), ("boxlib", "tfromp"),
                       units = "K")
        else:
            self.alias(("gas", "temperature"), ("boxlib", "tfromh"),
                       units = "K")

        for _, field in self.pf.field_list:
            if field.startswith("X("):
                # We have a fraction
                nice_name = field[2:-1]
                self.alias(("gas", "%s_fraction" % nice_name), ("boxlib", field),
                           units = "")
                def _create_density_func(field_name):
                    def _func(field, data):
                        return data[field_name] * data["gas", "density"]
                    return _func
                func = _create_density_func(("gas", "%s_fraction" % nice_name))
                self.add_field(name = ("gas", "%s_density" % nice_name),
                               function = func,
                               units = "g/cm**3")
                # We know this will either have one letter, or two.
                if field[3] in string.letters:
                    element, weight = field[2:4], field[4:-1]
                else:
                    element, weight = field[2:3], field[3:-1]
                weight = int(weight)
                # Here we can, later, add number density.
            if field.startswith("omegadot("):
                nice_name = field[9:-1]
                self.add_output_field(("boxlib", field), units = "1/s")
                self.alias(("gas", "%s_creation_rate" % nice_name),
                           ("boxlib", field), units = "1/s")
