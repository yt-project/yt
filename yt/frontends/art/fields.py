"""
ART-specific fields



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

b_units = "code_magnetic"
ra_units = "code_length / code_time**2"
rho_units = "code_mass / code_length**3"
vel_units = "code_velocity"
# NOTE: ARTIO uses momentum density.
mom_units = "code_mass / (code_length**2 * code_time)"
en_units = "code_mass*code_velocity**2/code_length**3"

class ARTFieldInfo(FieldInfoContainer):
    known_other_fields = (
        ("Density", (rho_units, ["density"], None)),
        ("TotalEnergy", (en_units, ["total_energy"], None)),
        ("XMomentumDensity", (mom_units, ["momentum_x"], None)),
        ("YMomentumDensity", (mom_units, ["momentum_y"], None)),
        ("ZMomentumDensity", (mom_units, ["momentum_z"], None)),
        ("Pressure", ("", ["pressure"], None)), # Unused
        ("Gamma", ("", ["gamma"], None)),
        ("GasEnergy", (en_units, ["thermal_energy"], None)),
        ("MetalDensitySNII", (rho_units, ["metal_ii_density"], None)),
        ("MetalDensitySNIa", (rho_units, ["metal_ia_density"], None)),
        ("PotentialNew", ("", ["potential"], None)),
        ("PotentialOld", ("", ["gas_potential"], None)),
    )

    known_particle_fields = (
        ("particle_position_x", ("code_length", [], None)),
        ("particle_position_y", ("code_length", [], None)),
        ("particle_position_z", ("code_length", [], None)),
        ("particle_velocity_x", (vel_units, [], None)),
        ("particle_velocity_y", (vel_units, [], None)),
        ("particle_velocity_z", (vel_units, [], None)),
        ("particle_mass", ("code_mass", [], None)),
        ("particle_index", ("", [], None)),
        ("particle_species", ("", ["particle_type"], None)),
        ("particle_creation_time", ("Gyr", [], None)),
        ("particle_mass_initial", ("code_mass", [], None)),
        ("particle_metallicity1", ("", [], None)),
        ("particle_metallicity2", ("", [], None)),
    )

    def setup_fluid_fields(self):
        unit_system = self.ds.unit_system
        def _temperature(field, data):
            r0 = data.ds.parameters['boxh'] / data.ds.parameters['ng']
            tr = data.ds.quan(3.03e5 * r0**2, 'K/code_velocity**2')
            tr *= data.ds.parameters['wmu'] * data.ds.parameters['Om0']
            tr *= (data.ds.parameters['gamma'] - 1.)
            tr /= data.ds.parameters['aexpn']**2
            return tr * data['art', 'GasEnergy'] / data['art', 'Density']
        self.add_field(('gas', 'temperature'),
                       function=_temperature, 
                       units=unit_system["temperature"])

        def _get_vel(axis):
            def velocity(field, data):
                return (data[('gas','momentum_%s' % axis)] /
                        data[('gas','density')])
            return velocity
        for ax in 'xyz':
            self.add_field(('gas','velocity_%s' % ax),
                           function = _get_vel(ax),
                           units=unit_system["velocity"])

        def _momentum_magnitude(field, data):
            tr = (data['gas','momentum_x']**2 +
                  data['gas','momentum_y']**2 +
                  data['gas','momentum_z']**2)**0.5
            tr *= data['index','cell_volume'].in_units('cm**3')
            return tr
        self.add_field(('gas', 'momentum_magnitude'),
                       function=_momentum_magnitude,
                       units=unit_system["momentum"])

        def _velocity_magnitude(field, data):
            tr = data['gas','momentum_magnitude']
            tr /= data['gas','cell_mass']
            return tr
        self.add_field(('gas', 'velocity_magnitude'),
                       function=_velocity_magnitude,
                       units=unit_system["velocity"])

        def _metal_density(field, data):
            tr = data['gas','metal_ia_density']
            tr += data['gas','metal_ii_density']
            return tr
        self.add_field(('gas','metal_density'),
                       function=_metal_density,
                       units=unit_system["density"])

        def _metal_mass_fraction(field, data):
            tr = data['gas','metal_density']
            tr /= data['gas','density']
            return tr
        self.add_field(('gas', 'metal_mass_fraction'),
                       function=_metal_mass_fraction,
                       units='')

        def _H_mass_fraction(field, data):
            tr = (1. - data.ds.parameters['Y_p'] - 
                  data['gas', 'metal_mass_fraction'])
            return tr
        self.add_field(('gas', 'H_mass_fraction'),
                       function=_H_mass_fraction,
                       units='')

        def _metallicity(field, data):
            tr = data['gas','metal_mass_fraction']
            tr /= data['gas','H_mass_fraction']
            return tr
        self.add_field(('gas', 'metallicity'),
                       function=_metallicity,
                       units='')
