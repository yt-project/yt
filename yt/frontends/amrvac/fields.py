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

from numpy import sqrt
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


        # kinetic energy density depends on velocities
        def _kinetic_energy_density(field, data):
            # devnote : have a look at issue 1301
            return 0.5 * data['gas', 'density'] * data['gas', 'velocity_magnitude']**2

        self.add_field(("gas", "kinetic_energy_density"), function=_kinetic_energy_density,
                        units=us["density"]*us["velocity"]**2,
                        dimensions=dimensions.density*dimensions.velocity**2,
                        sampling_type="cell")


        # magnetic energy
        def _magnetic_energy_density(field, data):
            emag = 0
            for idim in '123':
                if not ('amrvac', 'b%s' % idim) in self.field_list:
                    break
                emag += 0.5 * data['gas', 'magnetic_%s' % idim]**2
            return emag

        if ('amrvac', 'b1') in self.field_list:
            self.add_field(('gas', 'magnetic_energy_density'), function=_magnetic_energy_density,
                           units=us['density']*us['velocity']**2,
                           dimensions=dimensions.density*dimensions.velocity**2,
                           sampling_type='cell')


        # Adding the thermal pressure field.
        # In AMRVAC we have multiple physics possibilities:
        # - if HD/MHD + energy equation, pressure is (gamma-1)*(e - ekin (- emag)) for (M)HD
        # - if HD/MHD but solve_internal_e is true in parfile, pressure is (gamma-1)*e for both
        # - if (m)hd_energy is false in parfile (isothermal), pressure is c_adiab * rho**gamma
        def _full_thermal_pressure(field, data):
            # important note : energy density and pressure are actually expressed in the same unit
            pthermal = (data.ds.gamma - 1) * (data['gas', 'energy_density'] - data['gas', 'kinetic_energy_density'])
            if ('amrvac', 'b1') in self.field_list:
                pthermal -= (data.ds.gamma - 1) * data['gas', 'magnetic_energy_density']
            return pthermal

        def _polytropic_thermal_pressure(field, data):
            return (data.ds.gamma - 1) * data['gas', 'energy_density']

        def _adiabatic_thermal_pressure(field, data):
            return data.ds._c_adiab * data["gas", "density"]**data.ds.gamma

        pressure_is_defined = False
        pressure_equation = None
        if ("amrvac", "e") in self.field_list:
            if self.ds.namelist and self.ds.namelist['solve_internal_e']:
                pressure_equation = _polytropic_thermal_pressure
                mylog.info('Using polytropic EoS for thermal pressure.')
            else:
                pressure_equation = _full_thermal_pressure
                mylog.info('Using full (M)HD energy for thermal pressure.')
        elif self.ds._c_adiab is not None:
            pressure_equation = _adiabatic_thermal_pressure
            mylog.info('Using adiabatic EoS for thermal pressure (isothermal).')
            mylog.warning('If you used usr_set_pthermal you should redefine the thermal_pressure field.')

        if pressure_equation:
            self.add_field(('gas', 'thermal_pressure'), function=pressure_equation,
                           units=us['density']*us['velocity']**2,
                           dimensions=dimensions.density*dimensions.velocity**2,
                           sampling_type='cell')
            pressure_is_defined = True
        else:
            mylog.warning("e not found and no parfile passed, can not set thermal_pressure.")


        # fields that depend on thermal pressure
        if pressure_is_defined:
            def _sound_speed(field, data):
                return sqrt(data.ds.gamma * data["thermal_pressure"] / data["gas", "density"])

            self.add_field(("gas", "sound_speed"), function=_sound_speed,
                            units=us["velocity"],
                            dimensions=dimensions.velocity,
                            sampling_type="cell")

            # temperature depends on thermal pressure through T = p/rho
            def _temperature(field, data):
                return data['thermal_pressure'] / data['density']

            self.add_field(('gas', 'temperature'), function=_temperature,
                           units=us['temperature'], dimensions=dimensions.temperature, sampling_type='cell')
