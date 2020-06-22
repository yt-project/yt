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
                       sampling_type="cell",
                       function=_temperature, 
                       units=unit_system["temperature"])

        def _get_vel(axis):
            def velocity(field, data):
                return (data[('gas','momentum_%s' % axis)] /
                        data[('gas','density')])
            return velocity
        for ax in 'xyz':
            self.add_field(('gas','velocity_%s' % ax),
                           sampling_type="cell",
                           function = _get_vel(ax),
                           units=unit_system["velocity"])

        def _momentum_magnitude(field, data):
            tr = (data['gas','momentum_x']**2 +
                  data['gas','momentum_y']**2 +
                  data['gas','momentum_z']**2)**0.5
            tr *= data['index','cell_volume'].in_units('cm**3')
            return tr
        self.add_field(('gas', 'momentum_magnitude'),
                       sampling_type="cell",
                       function=_momentum_magnitude,
                       units=unit_system["momentum"])

        def _velocity_magnitude(field, data):
            tr = data['gas','momentum_magnitude']
            tr /= data['gas','cell_mass']
            return tr
        self.add_field(('gas', 'velocity_magnitude'),
                       sampling_type="cell",
                       function=_velocity_magnitude,
                       units=unit_system["velocity"])

        def _metal_density(field, data):
            tr = data['gas','metal_ia_density']
            tr += data['gas','metal_ii_density']
            return tr
        self.add_field(('gas','metal_density'),
                       sampling_type="cell",
                       function=_metal_density,
                       units=unit_system["density"])

        def _metal_mass_fraction(field, data):
            tr = data['gas','metal_density']
            tr /= data['gas','density']
            return tr
        self.add_field(('gas', 'metal_mass_fraction'),
                       sampling_type="cell",
                       function=_metal_mass_fraction,
                       units='')

        def _H_mass_fraction(field, data):
            tr = (1. - data.ds.parameters['Y_p'] - 
                  data['gas', 'metal_mass_fraction'])
            return tr
        self.add_field(('gas', 'H_mass_fraction'),
                       sampling_type="cell",
                       function=_H_mass_fraction,
                       units='')

        def _metallicity(field, data):
            tr = data['gas','metal_mass_fraction']
            tr /= data['gas','H_mass_fraction']
            return tr
        self.add_field(('gas', 'metallicity'),
                       sampling_type="cell",
                       function=_metallicity,
                       units='')

        atoms = ['C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', \
        'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', \
        'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']
        def _specific_metal_density_function(atom):
            def _specific_metal_density(field, data):
                nucleus_densityIa = data['gas','metal_ia_density']*\
                                    SNIa_abundance[atom]
                nucleus_densityII = data['gas','metal_ii_density']*\
                                    SNII_abundance[atom]
                return nucleus_densityIa + nucleus_densityII
            return _specific_metal_density
        for atom in atoms:
            self.add_field(('gas','%s_nuclei_mass_density'%atom),
                            sampling_type="cell",
                            function=_specific_metal_density_function(atom),
                            units=unit_system["density"])


# based on Iwamoto et al 1999
# mass fraction of each atom in SNIa metal 
SNIa_abundance = {
    'H'  : 0.00E+00, 'He' : 0.00E+00, 'C'  : 3.52E-02, 
    'N'  : 8.47E-07, 'O'  : 1.04E-01, 'F'  : 4.14E-10, 
    'Ne' : 3.30E-03, 'Na' : 4.61E-05, 'Mg' : 6.25E-03, 
    'Al' : 7.19E-04, 'Si' : 1.14E-01, 'P'  : 2.60E-04, 
    'S'  : 6.35E-02, 'Cl' : 1.27E-04, 'Ar' : 1.14E-02, 
    'K'  : 5.72E-05, 'Ca' : 8.71E-03, 'Sc' : 1.61E-07, 
    'Ti' : 2.50E-04, 'V'  : 5.46E-05, 'Cr' : 6.19E-03, 
    'Mn' : 6.47E-03, 'Fe' : 5.46E-01, 'Co' : 7.59E-04, 
    'Ni' : 9.17E-02, 'Cu' : 2.19E-06, 'Zn' : 2.06E-05}

# mass fraction of each atom in SNII metal 
SNII_abundance = {
    'H'  : 0.00E+00, 'He' : 0.00E+00, 'C'  : 3.12E-02, 
    'N'  : 6.15E-04, 'O'  : 7.11E-01, 'F'  : 4.57E-10, 
    'Ne' : 9.12E-02, 'Na' : 2.56E-03, 'Mg' : 4.84E-02, 
    'Al' : 5.83E-03, 'Si' : 4.81E-02, 'P'  : 4.77E-04, 
    'S'  : 1.62E-02, 'Cl' : 4.72E-05, 'Ar' : 3.15E-03, 
    'K'  : 2.65E-05, 'Ca' : 2.31E-03, 'Sc' : 9.02E-08, 
    'Ti' : 5.18E-05, 'V'  : 3.94E-06, 'Cr' : 5.18E-04, 
    'Mn' : 1.52E-04, 'Fe' : 3.58E-02, 'Co' : 2.86E-05, 
    'Ni' : 2.35E-03, 'Cu' : 4.90E-07, 'Zn' : 7.46E-06}
