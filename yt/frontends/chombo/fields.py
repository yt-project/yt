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
from yt.units.unit_object import Unit
from yt.fields.field_info_container import \
    FieldInfoContainer, \
    particle_deposition_functions, \
    particle_vector_functions, \
    standard_particle_fields

from yt.utilities.exceptions import YTFieldNotFound

rho_units = "code_mass / code_length**3"
mom_units = "code_mass / (code_time * code_length**2)"
eden_units = "code_mass / (code_time**2 * code_length)" # erg / cm^3
vel_units = "code_length / code_time"
b_units = "code_magnetic"


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
        ("X-magnfield", (b_units, [], None)),
        ("Y-magnfield", (b_units, [], None)),
        ("Z-magnfield", (b_units, [], None)),
        ("directrad-dedt-density", (eden_units, ["directrad-dedt-density"], None)),
        ("directrad-dpxdt-density", (mom_units, ["directrad-dpxdt-density"], None)),
        ("directrad-dpydt-density", (mom_units, ["directrad-dpydt-density"], None)),
        ("directrad-dpzdt-density", (mom_units, ["directrad-dpzdt-density"], None)),
        )
    known_particle_fields = (
        ("particle_mass", ("code_mass", [], None)),
        ("particle_position_x", ("code_length", [], None)),
        ("particle_position_y", ("code_length", [], None)),
        ("particle_position_z", ("code_length", [], None)),
        ("particle_momentum_x", ("code_mass*code_length/code_time", [], None)),
        ("particle_momentum_y", ("code_mass*code_length/code_time", [], None)),
        ("particle_momentum_z", ("code_mass*code_length/code_time", [], None)),
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

    def setup_particle_fields(self, ptype):

        def _get_vel(axis):
            def velocity(field, data):
                return data["particle_momentum_%s" % axis]/data["particle_mass"]
            return velocity

        for ax in 'xyz':
            self.add_field((ptype, "particle_velocity_%s" % ax), 
                           function=_get_vel(ax),
                           particle_type=True,
                           units="code_length/code_time")

        super(Orion2FieldInfo, self).setup_particle_fields(ptype)

    def setup_fluid_fields(self):
        from yt.fields.magnetic_field import \
            setup_magnetic_field_aliases
        unit_system = self.ds.unit_system
        def _thermal_energy_density(field, data):
            try:
                return data["energy-density"] - data["kinetic_energy"] - \
                    data["magnetic_energy"]
            except YTFieldNotFound:
                return data['energy-density'] - data["kinetic_energy"]

        def _thermal_energy(field, data):
            return data['thermal_energy_density']/data['density']

        def _magnetic_energy(field, data):
            ret = data["X-magnfield"]**2
            if data.ds.dimensionality > 1:
                ret = ret + data["Y-magnfield"]**2
            if data.ds.dimensionality > 2:
                ret = ret + data["Z-magnfield"]**2
            return ret/8.0/np.pi

        def _specific_magnetic_energy(field, data):
            return data['specific_magnetic_energy']/data['density']

        def _kinetic_energy(field, data):
            p2 = data['X-momentum']**2
            if data.ds.dimensionality > 1:
                p2 = p2 + data["Y-momentum"]**2
            if data.ds.dimensionality > 2:
                p2 = p2 + data["Z-momentum"]**2
            return 0.5 * p2/data['density']

        def _specific_kinetic_energy(field, data):
            return data['kinetic_energy']/data['density']

        def _temperature(field, data):
            c_v = data.ds.quan(data.ds.parameters['radiation.const_cv'], 
                               'erg/g/K')
            return (data["thermal_energy"]/c_v)

        def _get_vel(axis):
            def velocity(field, data):
                return data["momentum_%s" % axis]/data["density"]
            return velocity

        for ax in 'xyz':
            self.add_field(("gas", "velocity_%s" % ax), function = _get_vel(ax),
                           units = unit_system["velocity"])
        self.add_field(("gas", "thermal_energy"),
                       function = _thermal_energy,
                       units = unit_system["specific_energy"])
        self.add_field(("gas", "thermal_energy_density"),
                       function = _thermal_energy_density,
                       units = unit_system["pressure"])
        self.add_field(("gas", "kinetic_energy"),
                       function = _kinetic_energy,
                       units = unit_system["pressure"])
        self.add_field(("gas", "specific_kinetic_energy"),
                       function = _specific_kinetic_energy,
                       units = unit_system["specific_energy"])
        self.add_field(("gas", "magnetic_energy"),
                       function = _magnetic_energy,
                       units = unit_system["pressure"])
        self.add_field(("gas", "specific_magnetic_energy"),
                       function = _specific_magnetic_energy,
                       units = unit_system["specific_energy"])
        self.add_field(("gas", "temperature"), function=_temperature,
                       units=unit_system["temperature"])

        setup_magnetic_field_aliases(self, "chombo", ["%s-magnfield" % ax for ax in "XYZ"])


class ChomboPICFieldInfo3D(FieldInfoContainer):
    known_other_fields = (
        ("density", (rho_units, ["density", "Density"], None)),
        ("potential", ("code_length**2 / code_time**2", ["potential", "Potential"], None)),
        ("gravitational_field_x", ("code_length / code_time**2", [], None)),
        ("gravitational_field_y", ("code_length / code_time**2", [], None)),
        ("gravitational_field_z", ("code_length / code_time**2", [], None)),
    )
    known_particle_fields = (
        ("particle_mass", ("code_mass", [], None)),
        ("particle_position_x", ("code_length", [], None)),
        ("particle_position_y", ("code_length", [], None)),
        ("particle_position_z", ("code_length", [], None)),
        ("particle_velocity_x", ("code_length / code_time", [], None)),
        ("particle_velocity_y", ("code_length / code_time", [], None)),
        ("particle_velocity_z", ("code_length / code_time", [], None)),
    )

    # I am re-implementing this here to override a few default behaviors:
    # I don't want to skip output units for code_length and I do want
    # particle_fields to default to take_log = False.
    def setup_particle_fields(self, ptype, ftype='gas', num_neighbors=64 ):
        skip_output_units = ()
        for f, (units, aliases, dn) in sorted(self.known_particle_fields):
            units = self.ds.field_units.get((ptype, f), units)
            if (f in aliases or ptype not in self.ds.particle_types_raw) and \
                units not in skip_output_units:
                u = Unit(units, registry = self.ds.unit_registry)
                output_units = str(u.get_cgs_equivalent())
            else:
                output_units = units
            if (ptype, f) not in self.field_list:
                continue
            self.add_output_field((ptype, f),
                units = units, particle_type = True,
                display_name = dn, output_units = output_units, take_log=False)
            for alias in aliases:
                self.alias((ptype, alias), (ptype, f), units = output_units)

        ppos_fields = ["particle_position_%s" % ax for ax in 'xyz']
        pvel_fields = ["particle_velocity_%s" % ax for ax in 'xyz']
        particle_vector_functions(ptype, ppos_fields, pvel_fields, self)

        particle_deposition_functions(ptype, "particle_position",
            "particle_mass", self)
        standard_particle_fields(self, ptype)
        # Now we check for any leftover particle fields
        for field in sorted(self.field_list):
            if field in self: continue
            if not isinstance(field, tuple):
                raise RuntimeError
            if field[0] not in self.ds.particle_types:
                continue
            self.add_output_field(field, 
                                  units = self.ds.field_units.get(field, ""),
                                  particle_type = True)
        self.setup_smoothed_fields(ptype,
                                   num_neighbors=num_neighbors,
                                   ftype=ftype)


def _dummy_position(field, data):
    return 0.5*np.ones_like(data['particle_position_x'])


def _dummy_velocity(field, data):
    return np.zeros_like(data['particle_velocity_x'])


def _dummy_field(field, data):
    return 0.0 * data['gravitational_field_x']

fluid_field_types = ['chombo', 'gas']
particle_field_types = ['io', 'all']


class ChomboPICFieldInfo2D(ChomboPICFieldInfo3D):
    known_other_fields = (
        ("density", (rho_units, ["density", "Density"], None)),
        ("potential", ("code_length**2 / code_time**2", ["potential", "Potential"], None)),
        ("gravitational_field_x", ("code_length / code_time**2", [], None)),
        ("gravitational_field_y", ("code_length / code_time**2", [], None)),
    )
    known_particle_fields = (
        ("particle_mass", ("code_mass", [], None)),
        ("particle_position_x", ("code_length", [], None)),
        ("particle_position_y", ("code_length", [], None)),
        ("particle_velocity_x", ("code_length / code_time", [], None)),
        ("particle_velocity_y", ("code_length / code_time", [], None)),
    )

    def __init__(self, ds, field_list):
        super(ChomboPICFieldInfo2D, self).__init__(ds, field_list)

        for ftype in fluid_field_types:
            self.add_field((ftype, 'gravitational_field_z'), function = _dummy_field, 
                            units = "code_length / code_time**2")

        for ptype in particle_field_types:                
            self.add_field((ptype, "particle_position_z"), function = _dummy_position,
                           particle_type = True,
                           units = "code_length")

            self.add_field((ptype, "particle_velocity_z"), function = _dummy_velocity,
                           particle_type = True,
                           units = "code_length / code_time")


class ChomboPICFieldInfo1D(ChomboPICFieldInfo3D):
    known_other_fields = (
        ("density", (rho_units, ["density", "Density"], None)),
        ("potential", ("code_length**2 / code_time**2", ["potential", "Potential"], None)),
        ("gravitational_field_x", ("code_length / code_time**2", [], None)),
    )
    known_particle_fields = (
        ("particle_mass", ("code_mass", [], None)),
        ("particle_position_x", ("code_length", [], None)),
        ("particle_velocity_x", ("code_length / code_time", [], None)),
    )

    def __init__(self, ds, field_list):
        super(ChomboPICFieldInfo1D, self).__init__(ds, field_list)

        for ftype in fluid_field_types:
            self.add_field((ftype, 'gravitational_field_y'), function = _dummy_field, 
                            units = "code_length / code_time**2")

            self.add_field((ftype, 'gravitational_field_z'), function = _dummy_field, 
                    units = "code_length / code_time**2")

        for ptype in particle_field_types:
            self.add_field((ptype, "particle_position_y"), function = _dummy_position,
                           particle_type = True,
                           units = "code_length")
            self.add_field((ptype, "particle_position_z"), function = _dummy_position,
                           particle_type = True,
                           units = "code_length")
            self.add_field((ptype, "particle_velocity_y"), function = _dummy_velocity,
                           particle_type = True,
                           units = "code_length / code_time")
            self.add_field((ptype, "particle_velocity_z"), function = _dummy_velocity,
                           particle_type = True,
                           units = "code_length / code_time")


class PlutoFieldInfo(ChomboFieldInfo):
    known_other_fields = (
        ("rho", (rho_units, ["density"], None)),
        ("prs", ("code_mass / (code_length * code_time**2)", ["pressure"], None)),
        ("vx1", (vel_units, ["velocity_x"], None)),
        ("vx2", (vel_units, ["velocity_y"], None)),
        ("vx3", (vel_units, ["velocity_z"], None)),
        ("bx1", (b_units, [], None)),
        ("bx2", (b_units, [], None)),
        ("bx3", (b_units, [], None)),
    )

    known_particle_fields = ()

    def setup_fluid_fields(self):
        from yt.fields.magnetic_field import \
            setup_magnetic_field_aliases
        setup_magnetic_field_aliases(self, "chombo", ["bx%s" % ax for ax in [1,2,3]])

