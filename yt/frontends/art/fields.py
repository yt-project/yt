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
from yt.utilities.physical_constants import mh

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
        self.add_field(('gas', 'temperature'), sampling_type="cell", 
                       function=_temperature, 
                       units=unit_system["temperature"])

        def _get_vel(axis):
            def velocity(field, data):
                return (data[('gas','momentum_%s' % axis)] /
                        data[('gas','density')])
            return velocity
        for ax in 'xyz':
            self.add_field(('gas','velocity_%s' % ax), sampling_type="cell", 
                           function = _get_vel(ax),
                           units=unit_system["velocity"])

        def _momentum_magnitude(field, data):
            tr = (data['gas','momentum_x']**2 +
                  data['gas','momentum_y']**2 +
                  data['gas','momentum_z']**2)**0.5
            tr *= data['index','cell_volume'].in_units('cm**3')
            return tr
        self.add_field(('gas', 'momentum_magnitude'), sampling_type="cell", 
                       function=_momentum_magnitude,
                       units=unit_system["momentum"])

        def _velocity_magnitude(field, data):
            tr = data['gas','momentum_magnitude']
            tr /= data['gas','cell_mass']
            return tr
        self.add_field(('gas', 'velocity_magnitude'), sampling_type="cell", 
                       function=_velocity_magnitude,
                       units=unit_system["velocity"])

        def _metal_density(field, data):
            tr = data['gas','metal_ia_density']
            tr += data['gas','metal_ii_density']
            return tr
        self.add_field(('gas','metal_density'), sampling_type="cell", 
                       function=_metal_density,
                       units=unit_system["density"])

        def _metal_mass_fraction(field, data):
            tr = data['gas','metal_density']
            tr /= data['gas','density']
            return tr
        self.add_field(('gas', 'metal_mass_fraction'), sampling_type="cell", 
                       function=_metal_mass_fraction,
                       units='')

        def _H_mass_fraction(field, data):
            tr = (1. - data.ds.parameters['Y_p'] - 
                  data['gas', 'metal_mass_fraction'])
            return tr
        self.add_field(('gas', 'H_mass_fraction'), sampling_type="cell", 
                       function=_H_mass_fraction,
                       units='')

        def _metallicity(field, data):
            tr = data['gas','metal_mass_fraction']
            tr /= data['gas','H_mass_fraction']
            return tr
        self.add_field(('gas', 'metallicity'), sampling_type="cell", 
                       function=_metallicity,
                       units='')

        atoms = ['C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', \
        'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', \
        'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']
        for atom in atoms:
            def _specific_metal_density(field, data):
                nucleus_densityIa = data[('gas','metal_ia_density')].in_units("g / cm**3")*\
                                    data.ds.quan(SNIa_abundance[atom],"1 / g")*atomic_mass[atom]*mh
                nucleus_densityII = data[('gas','metal_ii_density')].in_units("g / cm**3")*\
                                    data.ds.quan(SNIa_abundance[atom],"1 / g")*atomic_mass[atom]*mh
                return nucleus_densityIa + nucleus_densityII
            self.add_field(('gas','%s_nuclei_mass_density'%atom),sampling_type="cell",
                            function=_specific_metal_density,
                            units=unit_system["density"])



# based on Iwamoto et al 1999
# number of atoms per gram of SNIa metal
SNIa_abundance = {
    'H'  : 0.00E+00, 'He' : 0.00E+00, 'C'  : 1.75E+21, 
    'N'  : 3.61E+16, 'O'  : 3.90E+21, 'F'  : 1.30E+13, 
    'Ne' : 9.76E+19, 'Na' : 1.20E+18, 'Mg' : 1.54E+20, 
    'Al' : 1.59E+19, 'Si' : 2.43E+21, 'P'  : 5.02E+18, 
    'S'  : 1.18E+21, 'Cl' : 2.14E+18, 'Ar' : 1.71E+20, 
    'K'  : 8.74E+17, 'Ca' : 1.30E+20, 'Sc' : 2.14E+15, 
    'Ti' : 3.12E+18, 'V'  : 6.41E+17, 'Cr' : 7.11E+19, 
    'Mn' : 7.04E+19, 'Fe' : 5.85E+21, 'Co' : 7.69E+18, 
    'Ni' : 9.34E+20, 'Cu' : 2.06E+16, 'Zn' : 1.88E+17}

# number of atoms per gram of SNII metal
SNII_abundance = {
    'H'  : 0.00E+00, 'He' : 0.00E+00, 'C'  : 1.55E+21, 
    'N'  : 2.62E+19, 'O'  : 2.66E+22, 'F'  : 1.44E+13, 
    'Ne' : 2.70E+21, 'Na' : 6.67E+19, 'Mg' : 1.19E+21, 
    'Al' : 1.29E+20, 'Si' : 1.02E+21, 'P'  : 9.20E+18, 
    'S'  : 3.02E+20, 'Cl' : 7.95E+17, 'Ar' : 4.71E+19, 
    'K'  : 4.06E+17, 'Ca' : 3.45E+19, 'Sc' : 1.20E+15, 
    'Ti' : 6.47E+17, 'V'  : 4.62E+16, 'Cr' : 5.96E+18, 
    'Mn' : 1.65E+18, 'Fe' : 3.82E+20, 'Co' : 2.90E+17, 
    'Ni' : 2.40E+19, 'Cu' : 4.61E+15, 'Zn' : 6.81E+16}

# taken from TRIDENT
atomic_mass = {
    'H' : 1.00794,   'He': 4.002602,  'Li': 6.941,
    'Be': 9.012182,  'B' : 10.811,    'C' : 12.0107,
    'N' : 14.0067,   'O' : 15.9994,   'F' : 18.9984032,
    'Ne': 20.1797,   'Na': 22.989770, 'Mg': 24.3050,
    'Al': 26.981538, 'Si': 28.0855,   'P' : 30.973761,
    'S' : 32.065,    'Cl': 35.453,    'Ar': 39.948,
    'K' : 39.0983,   'Ca': 40.078,    'Sc': 44.955910,
    'Ti': 47.867,    'V' : 50.9415,   'Cr': 51.9961,
    'Mn': 54.938049, 'Fe': 55.845,    'Co': 58.933200,
    'Ni': 58.6934,   'Cu': 63.546,    'Zn': 65.409}
