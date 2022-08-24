import numpy as np

from yt._typing import KnownFieldsT
from yt.fields.field_info_container import FieldInfoContainer
from yt.units.yt_array import YTArray
from yt.utilities.physical_constants import amu_cgs, boltzmann_constant_cgs

b_units = "code_magnetic"
ra_units = "code_length / code_time**2"
rho_units = "code_mass / code_length**3"
vel_units = "code_velocity"
# NOTE: ARTIO uses momentum density.
mom_units = "code_mass / (code_length**2 * code_time)"
en_units = "code_mass*code_velocity**2/code_length**3"
p_units = "code_mass / (code_length * code_time**2)"


class ARTIOFieldInfo(FieldInfoContainer):
    known_other_fields: KnownFieldsT = (
        ("HVAR_GAS_DENSITY", (rho_units, ["density"], None)),
        ("HVAR_GAS_ENERGY", (en_units, ["total_energy_density"], None)),
        ("HVAR_INTERNAL_ENERGY", (en_units, ["thermal_energy_density"], None)),
        ("HVAR_PRESSURE", (p_units, ["pressure"], None)),
        ("HVAR_MOMENTUM_X", (mom_units, ["momentum_density_x"], None)),
        ("HVAR_MOMENTUM_Y", (mom_units, ["momentum_density_y"], None)),
        ("HVAR_MOMENTUM_Z", (mom_units, ["momentum_density_z"], None)),
        ("HVAR_GAMMA", ("", ["gamma"], None)),
        ("HVAR_METAL_DENSITY_Ia", (rho_units, ["metal_ia_density"], None)),
        ("HVAR_METAL_DENSITY_II", (rho_units, ["metal_ii_density"], None)),
        ("VAR_POTENTIAL", ("", ["potential"], None)),
        ("VAR_POTENTIAL_HYDRO", ("", ["gas_potential"], None)),
        ("RT_HVAR_HI", (rho_units, ["H_density"], None)),
        ("RT_HVAR_HII", (rho_units, ["H_p1_density"], None)),
        ("RT_HVAR_H2", (rho_units, ["H2_density"], None)),
        ("RT_HVAR_HeI", (rho_units, ["He_density"], None)),
        ("RT_HVAR_HeII", (rho_units, ["He_p1_density"], None)),
        ("RT_HVAR_HeIII", (rho_units, ["He_p2_density"], None)),
    )

    known_particle_fields: KnownFieldsT = (
        ("POSITION_X", ("code_length", ["particle_position_x"], None)),
        ("POSITION_Y", ("code_length", ["particle_position_y"], None)),
        ("POSITION_Z", ("code_length", ["particle_position_z"], None)),
        ("VELOCITY_X", (vel_units, ["particle_velocity_x"], None)),
        ("VELOCITY_Y", (vel_units, ["particle_velocity_y"], None)),
        ("VELOCITY_Z", (vel_units, ["particle_velocity_z"], None)),
        ("MASS", ("code_mass", ["particle_mass"], None)),
        ("PID", ("", ["particle_index"], None)),
        ("SPECIES", ("", ["particle_type"], None)),
        ("BIRTH_TIME", ("", [], None)),  # code-units defined as dimensionless to
        # avoid incorrect conversion
        ("INITIAL_MASS", ("code_mass", ["initial_mass"], None)),
        ("METALLICITY_SNIa", ("", ["metallicity_snia"], None)),
        ("METALLICITY_SNII", ("", ["metallicity_snii"], None)),
    )

    def setup_fluid_fields(self):
        unit_system = self.ds.unit_system

        def _get_vel(axis):
            def velocity(field, data):
                return (
                    data[("gas", f"momentum_density_{axis}")] / data[("gas", "density")]
                )

            return velocity

        for ax in "xyz":
            self.add_field(
                ("gas", f"velocity_{ax}"),
                sampling_type="cell",
                function=_get_vel(ax),
                units=unit_system["velocity"],
            )

        def _temperature(field, data):
            tr = data[("gas", "thermal_energy_density")] / data[("gas", "density")]
            # We want this to match *exactly* what ARTIO would compute
            # internally.  We therefore use the exact values that are internal
            # to ARTIO, rather than yt's own internal constants.
            mH = 1.007825 * amu_cgs
            mHe = 4.002602 * amu_cgs
            Yp = 0.24
            XH = 1.0 - Yp
            XHe = 0.25 * Yp
            mb = XH * mH + XHe * mHe
            wmu = 4.0 / (8.0 - 5.0 * Yp)
            # Note that we have gamma = 5.0/3.0 here
            tr *= data[("gas", "gamma")] - 1.0
            tr *= wmu
            tr *= mb / boltzmann_constant_cgs
            return tr

        # TODO: The conversion factor here needs to be addressed, as previously
        # it was set as:
        # unit_T = unit_v**2.0*mb / constants.k
        self.add_field(
            ("gas", "temperature"),
            sampling_type="cell",
            function=_temperature,
            units=unit_system["temperature"],
        )

        # Create a metal_density field as sum of existing metal fields.
        flag1 = ("artio", "HVAR_METAL_DENSITY_Ia") in self.field_list
        flag2 = ("artio", "HVAR_METAL_DENSITY_II") in self.field_list
        if flag1 or flag2:
            if flag1 and flag2:

                def _metal_density(field, data):
                    tr = data[("gas", "metal_ia_density")].copy()
                    np.add(tr, data[("gas", "metal_ii_density")], out=tr)
                    return tr

            elif flag1 and not flag2:

                def _metal_density(field, data):
                    tr = data["metal_ia_density"]
                    return tr

            else:

                def _metal_density(field, data):
                    tr = data["metal_ii_density"]
                    return tr

            self.add_field(
                ("gas", "metal_density"),
                sampling_type="cell",
                function=_metal_density,
                units=unit_system["density"],
                take_log=True,
            )

    def setup_particle_fields(self, ptype):
        if ptype == "STAR":

            def _creation_time(field, data):
                return YTArray(
                    data.ds._handle.tphys_from_tcode_array(data["STAR", "BIRTH_TIME"]),
                    "yr",
                )

            def _age(field, data):
                return data.ds.current_time - data["STAR", "creation_time"]

            self.add_field(
                (ptype, "creation_time"),
                sampling_type="particle",
                function=_creation_time,
                units="yr",
            )
            self.add_field(
                (ptype, "age"), sampling_type="particle", function=_age, units="yr"
            )

            if self.ds.cosmological_simulation:

                def _creation_redshift(field, data):
                    return (
                        1.0
                        / data.ds._handle.auni_from_tcode_array(
                            data["STAR", "BIRTH_TIME"]
                        )
                        - 1.0
                    )

                self.add_field(
                    (ptype, "creation_redshift"),
                    sampling_type="particle",
                    function=_creation_redshift,
                )

        super().setup_particle_fields(ptype)
