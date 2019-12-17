"""
AMRVAC-specific fields

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
from yt.fields.magnetic_field import setup_magnetic_field_aliases
from yt.units import dimensions
from yt import mylog

# We need to specify which fields we might have in our dataset.  The field info
# container subclass here will define which fields it knows about.  There are
# optionally methods on it that get called which can be subclassed.

direction_aliases = {
    "cartesian": ("x", "y", "z"),
    "polar": ("r", "theta", "z"),
    "cylindrical": ("r", "z", "theta"),
    "spherical": ("r", "theta", "phi")
}

class AMRVACFieldInfo(FieldInfoContainer):
    code_density = "code_mass / code_length**3"
    code_moment = "code_mass / code_length**2 / code_time"
    code_pressure = "code_mass / code_length / code_time**2"

    # format: (native(?) field, (units, [aliases], display_name))
    # note: aliases will correspond to "gas" typed fields, whereas the native ones are "amrvac" typed
    known_other_fields = (
        ("rho", (code_density, ["density"], None)),
        ("m1", (code_moment, ["moment_1"], None)),
        ("m2", (code_moment, ["moment_2"], None)),
        ("m3", (code_moment, ["moment_3"], None)),
        ("e", (code_pressure, ["energy_density"], None)),
        ("b1", ("code_magnetic", ["magnetic_1"], None)),
        ("b2", ("code_magnetic", ["magnetic_2"], None)),
        ("b3", ("code_magnetic", ["magnetic_3"], None))
    )

    known_particle_fields = ()

    def setup_fluid_fields(self):
        def _v1(field, data):
            return data["gas", "moment_1"] / data["gas", "density"]
        def _v2(field, data):
            return data["gas", "moment_2"] / data["gas", "density"]
        def _v3(field, data):
            return data["gas", "moment_3"] / data["gas", "density"]

        us = self.ds.unit_system
        aliases = direction_aliases[self.ds.geometry]
        for idir, alias, func in zip("123", aliases, (_v1, _v2, _v3)):
            if not ("amrvac", "m%s" % idir) in self.field_list:
                break
            self.add_field(("gas", "velocity_%s" % alias), function=func,
                            units=us['velocity'],
                            dimensions=dimensions.velocity,
                            sampling_type="cell")
            self.alias(("gas", "velocity_%s" % idir), ("gas", "velocity_%s" % alias),
                        units=us["velocity"])
            self.alias(("gas", "moment_%s" % alias), ("gas", "moment_%s" % idir),
                        units=us["density"]*us["velocity"])

        setup_magnetic_field_aliases(self, "amrvac", ["mag%s" % ax for ax in "xyz"])


        # fields with nested dependencies are defined thereafter by increasing level of complexity


        # kinetic pressure is given by 0.5 * rho * v**2
        def _kinetic_energy_density(field, data):
            # devnote : have a look at issue 1301
            return 0.5 * data['gas', 'density'] * data['gas', 'velocity_magnitude']**2

        self.add_field(("gas", "kinetic_energy_density"), function=_kinetic_energy_density,
                        units=us["density"]*us["velocity"]**2,
                        dimensions=dimensions.density*dimensions.velocity**2,
                        sampling_type="cell")


        # magnetic energy density
        if ('amrvac', 'b1') in self.field_list:
            def _magnetic_energy_density(field, data):
                emag = 0.5 * data['gas', 'magnetic_1']**2
                for idim in '23':
                    if not ('amrvac', 'b%s' % idim) in self.field_list:
                        break
                    emag += 0.5 * data['gas', 'magnetic_%s' % idim]**2
                # important note: in AMRVAC the magnetic field is defined in units where mu0 = 1, such that
                # Emag = 0.5*B**2 instead of Emag = 0.5*B**2 / mu0
                # To correctly transform the dimensionality from gauss**2 -> rho*v**2, we have to take mu0 into account.
                # If we divide here, units when adding the field should be us["density"]*us["velocity"]**2
                # If not, they should be us["magnetic_field"]**2 and division should happen elsewhere.
                emag /= 4 * np.pi
                # divided by mu0 = 4pi in cgs, yt handles 'mks' and 'code' unit systems internally.
                return emag

            self.add_field(('gas', 'magnetic_energy_density'), function=_magnetic_energy_density,
                           units=us["density"] * us["velocity"] ** 2,
                           dimensions=dimensions.density * dimensions.velocity ** 2,
                           sampling_type='cell')

        # Adding the thermal pressure field.
        # In AMRVAC we have multiple physics possibilities:
        # - if HD/MHD + energy equation, pressure is (gamma-1)*(e - ekin (- emag)) for (M)HD
        # - if HD/MHD but solve_internal_e is true in parfile, pressure is (gamma-1)*e for both
        # - if (m)hd_energy is false in parfile (isothermal), pressure is c_adiab * rho**gamma


        def _full_thermal_pressure_HD(field, data):
            # important note : energy density and pressure are actually expressed in the same unit
            pthermal = (data.ds.gamma - 1) * (data['gas', 'energy_density'] - data['gas', 'kinetic_energy_density'])
            return pthermal

        def _full_thermal_pressure_MHD(field, data):
            pthermal = _full_thermal_pressure_HD(field, data) \
                       - (data.ds.gamma - 1) * data["gas", "magnetic_energy_density"]
            return pthermal

        def _polytropic_thermal_pressure(field, data):
            return (data.ds.gamma - 1) * data['gas', 'energy_density']

        def _adiabatic_thermal_pressure(field, data):
            return data.ds._c_adiab * data["gas", "density"]**data.ds.gamma

        pressure_recipe = None
        if ("amrvac", "e") in self.field_list:
            if self.ds._e_is_internal:
                pressure_recipe = _polytropic_thermal_pressure
                mylog.info('Using polytropic EoS for thermal pressure.')
            elif ('amrvac', 'b1') in self.field_list:
                pressure_recipe = _full_thermal_pressure_MHD
                mylog.info('Using full MHD energy for thermal pressure.')
            else:
                pressure_recipe = _full_thermal_pressure_HD
                mylog.info('Using full HD energy for thermal pressure.')
        elif self.ds._c_adiab is not None:
            pressure_recipe = _adiabatic_thermal_pressure
            mylog.info('Using adiabatic EoS for thermal pressure (isothermal).')
            mylog.warning('If you used usr_set_pthermal you should redefine the thermal_pressure field.')

        if pressure_recipe is not None:
            self.add_field(('gas', 'thermal_pressure'), function=pressure_recipe,
                           units=us['density']*us['velocity']**2,
                           dimensions=dimensions.density*dimensions.velocity**2,
                           sampling_type='cell')

            # sound speed and temperature depend on thermal pressure
            def _sound_speed(field, data):
                return np.sqrt(data.ds.gamma * data["gas", "thermal_pressure"] / data["gas", "density"])

            self.add_field(("gas", "sound_speed"), function=_sound_speed,
                            units=us["velocity"],
                            dimensions=dimensions.velocity,
                            sampling_type="cell")
        else:
            mylog.warning("e not found and no parfile passed, can not set thermal_pressure.")
