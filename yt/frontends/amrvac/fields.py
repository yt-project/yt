"""
AMRVAC-specific fields

"""


import functools

import numpy as np

from yt.fields.field_info_container import FieldInfoContainer
from yt.fields.magnetic_field import setup_magnetic_field_aliases
from yt.units import dimensions
from yt.utilities.logger import ytLogger as mylog

# We need to specify which fields we might have in our dataset.  The field info
# container subclass here will define which fields it knows about.  There are
# optionally methods on it that get called which can be subclassed.

direction_aliases = {
    "cartesian": ("x", "y", "z"),
    "polar": ("r", "theta", "z"),
    "cylindrical": ("r", "z", "theta"),
    "spherical": ("r", "theta", "phi"),
}


def _velocity(field, data, idir, prefix=None):
    """Velocity = linear momentum / density"""
    # This is meant to be used with functools.partial to produce
    # functions with only 2 arguments (field, data)
    # idir : int
    #    the direction index (1, 2 or 3)
    # prefix : str
    #    used to generalize to dust fields
    if prefix is None:
        prefix = ""
    moment = data["gas", "%smoment_%d" % (prefix, idir)]
    rho = data["gas", f"{prefix}density"]

    mask1 = rho == 0
    if mask1.any():
        mylog.info(
            "zeros found in %sdensity, "
            "patching them to compute corresponding velocity field.",
            prefix,
        )
        mask2 = moment == 0
        if not ((mask1 & mask2) == mask1).all():
            raise RuntimeError
        rho[mask1] = 1
    return moment / rho


code_density = "code_mass / code_length**3"
code_moment = "code_mass / code_length**2 / code_time"
code_pressure = "code_mass / code_length / code_time**2"


class AMRVACFieldInfo(FieldInfoContainer):
    # for now, define a finite family of dust fields (up to 100 species)
    MAXN_DUST_SPECIES = 100
    known_dust_fields = [
        ("rhod%d" % idust, (code_density, ["dust%d_density" % idust], None))
        for idust in range(1, MAXN_DUST_SPECIES + 1)
    ] + [
        (
            "m%dd%d" % (idir, idust),
            (code_moment, ["dust%d_moment_%d" % (idust, idir)], None),
        )
        for idust in range(1, MAXN_DUST_SPECIES + 1)
        for idir in (1, 2, 3)
    ]
    # format: (native(?) field, (units, [aliases], display_name))
    # note: aliases will correspond to "gas" typed fields
    # whereas the native ones are "amrvac" typed
    known_other_fields = (
        ("rho", (code_density, ["density"], None)),
        ("m1", (code_moment, ["moment_1"], None)),
        ("m2", (code_moment, ["moment_2"], None)),
        ("m3", (code_moment, ["moment_3"], None)),
        ("e", (code_pressure, ["energy_density"], None)),
        ("b1", ("code_magnetic", ["magnetic_1"], None)),
        ("b2", ("code_magnetic", ["magnetic_2"], None)),
        ("b3", ("code_magnetic", ["magnetic_3"], None)),
        ("Te", ("code_temperature", ["temperature"], None)),
        *known_dust_fields,
    )

    known_particle_fields = ()

    def _setup_velocity_fields(self, idust=None):
        if idust is None:
            dust_flag = dust_label = ""
        else:
            dust_flag = "d%d" % idust
            dust_label = "dust%d_" % idust

        us = self.ds.unit_system
        for idir, alias in enumerate(direction_aliases[self.ds.geometry], start=1):
            if not ("amrvac", "m%d%s" % (idir, dust_flag)) in self.field_list:
                break
            velocity_fn = functools.partial(_velocity, idir=idir, prefix=dust_label)
            self.add_field(
                ("gas", f"{dust_label}velocity_{alias}"),
                function=velocity_fn,
                units=us["velocity"],
                dimensions=dimensions.velocity,
                sampling_type="cell",
            )
            self.alias(
                ("gas", "%svelocity_%d" % (dust_label, idir)),
                ("gas", f"{dust_label}velocity_{alias}"),
                units=us["velocity"],
            )
            self.alias(
                ("gas", f"{dust_label}moment_{alias}"),
                ("gas", "%smoment_%d" % (dust_label, idir)),
                units=us["density"] * us["velocity"],
            )

    def _setup_dust_fields(self):
        idust = 1
        imax = self.__class__.MAXN_DUST_SPECIES
        while ("amrvac", "rhod%d" % idust) in self.field_list:
            if idust > imax:
                mylog.error(
                    "Only the first %d dust species are currently read by yt. "
                    "If you read this, please consider issuing a ticket. ",
                    imax,
                )
                break
            self._setup_velocity_fields(idust)
            idust += 1
        n_dust_found = idust - 1

        us = self.ds.unit_system
        if n_dust_found > 0:

            def _total_dust_density(field, data):
                tot = np.zeros_like(data[("gas", "density")])
                for idust in range(1, n_dust_found + 1):
                    tot += data["dust%d_density" % idust]
                return tot

            self.add_field(
                ("gas", "total_dust_density"),
                function=_total_dust_density,
                dimensions=dimensions.density,
                units=us["density"],
                sampling_type="cell",
            )

            def dust_to_gas_ratio(field, data):
                return data[("gas", "total_dust_density")] / data[("gas", "density")]

            self.add_field(
                ("gas", "dust_to_gas_ratio"),
                function=dust_to_gas_ratio,
                dimensions=dimensions.dimensionless,
                sampling_type="cell",
            )

    def setup_fluid_fields(self):

        setup_magnetic_field_aliases(self, "amrvac", [f"mag{ax}" for ax in "xyz"])
        self._setup_velocity_fields()  # gas velocities
        self._setup_dust_fields()  # dust derived fields (including velocities)

        # fields with nested dependencies are defined thereafter
        # by increasing level of complexity
        us = self.ds.unit_system

        def _kinetic_energy_density(field, data):
            # devnote : have a look at issue 1301
            return 0.5 * data["gas", "density"] * data["gas", "velocity_magnitude"] ** 2

        self.add_field(
            ("gas", "kinetic_energy_density"),
            function=_kinetic_energy_density,
            units=us["density"] * us["velocity"] ** 2,
            dimensions=dimensions.density * dimensions.velocity**2,
            sampling_type="cell",
        )

        # magnetic energy density
        if ("amrvac", "b1") in self.field_list:

            def _magnetic_energy_density(field, data):
                emag = 0.5 * data["gas", "magnetic_1"] ** 2
                for idim in "23":
                    if not ("amrvac", f"b{idim}") in self.field_list:
                        break
                    emag += 0.5 * data["gas", f"magnetic_{idim}"] ** 2
                # in AMRVAC the magnetic field is defined in units where mu0 = 1,
                # such that
                # Emag = 0.5*B**2 instead of Emag = 0.5*B**2 / mu0
                # To correctly transform the dimensionality from gauss**2 -> rho*v**2,
                # we have to take mu0 into account. If we divide here, units when adding
                # the field should be us["density"]*us["velocity"]**2.
                # If not, they should be us["magnetic_field"]**2 and division should
                # happen elsewhere.
                emag /= 4 * np.pi
                # divided by mu0 = 4pi in cgs,
                # yt handles 'mks' and 'code' unit systems internally.
                return emag

            self.add_field(
                ("gas", "magnetic_energy_density"),
                function=_magnetic_energy_density,
                units=us["density"] * us["velocity"] ** 2,
                dimensions=dimensions.density * dimensions.velocity**2,
                sampling_type="cell",
            )

        # Adding the thermal pressure field.
        # In AMRVAC we have multiple physics possibilities:
        # - if HD/MHD + energy equation P = (gamma-1)*(e - ekin (- emag)) for (M)HD
        # - if HD/MHD but solve_internal_e is true in parfile, P = (gamma-1)*e for both
        # - if (m)hd_energy is false in parfile (isothermal), P = c_adiab * rho**gamma

        def _full_thermal_pressure_HD(field, data):
            # energy density and pressure are actually expressed in the same unit
            pthermal = (data.ds.gamma - 1) * (
                data["gas", "energy_density"] - data["gas", "kinetic_energy_density"]
            )
            return pthermal

        def _full_thermal_pressure_MHD(field, data):
            pthermal = (
                _full_thermal_pressure_HD(field, data)
                - (data.ds.gamma - 1) * data["gas", "magnetic_energy_density"]
            )
            return pthermal

        def _polytropic_thermal_pressure(field, data):
            return (data.ds.gamma - 1) * data["gas", "energy_density"]

        def _adiabatic_thermal_pressure(field, data):
            return data.ds._c_adiab * data["gas", "density"] ** data.ds.gamma

        pressure_recipe = None
        if ("amrvac", "e") in self.field_list:
            if self.ds._e_is_internal:
                pressure_recipe = _polytropic_thermal_pressure
                mylog.info("Using polytropic EoS for thermal pressure.")
            elif ("amrvac", "b1") in self.field_list:
                pressure_recipe = _full_thermal_pressure_MHD
                mylog.info("Using full MHD energy for thermal pressure.")
            else:
                pressure_recipe = _full_thermal_pressure_HD
                mylog.info("Using full HD energy for thermal pressure.")
        elif self.ds._c_adiab is not None:
            pressure_recipe = _adiabatic_thermal_pressure
            mylog.info("Using adiabatic EoS for thermal pressure (isothermal).")
            mylog.warning(
                "If you used usr_set_pthermal you should "
                "redefine the thermal_pressure field."
            )

        if pressure_recipe is not None:
            self.add_field(
                ("gas", "thermal_pressure"),
                function=pressure_recipe,
                units=us["density"] * us["velocity"] ** 2,
                dimensions=dimensions.density * dimensions.velocity**2,
                sampling_type="cell",
            )

            # sound speed and temperature depend on thermal pressure
            def _sound_speed(field, data):
                return np.sqrt(
                    data.ds.gamma
                    * data["gas", "thermal_pressure"]
                    / data["gas", "density"]
                )

            self.add_field(
                ("gas", "sound_speed"),
                function=_sound_speed,
                units=us["velocity"],
                dimensions=dimensions.velocity,
                sampling_type="cell",
            )
        else:
            mylog.warning(
                "e not found and no parfile passed, can not set thermal_pressure."
            )
