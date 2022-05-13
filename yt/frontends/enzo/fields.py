import numpy as np

from yt._typing import KnownFieldsT
from yt.fields.field_info_container import FieldInfoContainer
from yt.utilities.physical_constants import me, mp

b_units = "code_magnetic"
e_units = "code_magnetic/c"
ra_units = "code_length / code_time**2"
rho_units = "code_mass / code_length**3"
vel_units = "code_velocity"

known_species_names = {
    "HI": "H_p0",
    "HII": "H_p1",
    "HeI": "He_p0",
    "HeII": "He_p1",
    "HeIII": "He_p2",
    "H2I": "H2_p0",
    "H2II": "H2_p1",
    "HM": "H_m1",
    "HeH": "HeH_p0",
    "DI": "D_p0",
    "DII": "D_p1",
    "HDI": "HD_p0",
    "Electron": "El",
    "OI": "O_p0",
    "OII": "O_p1",
    "OIII": "O_p2",
    "OIV": "O_p3",
    "OV": "O_p4",
    "OVI": "O_p5",
    "OVII": "O_p6",
    "OVIII": "O_p7",
    "OIX": "O_p8",
}

NODAL_FLAGS = {
    "BxF": [1, 0, 0],
    "ByF": [0, 1, 0],
    "BzF": [0, 0, 1],
    "Ex": [0, 1, 1],
    "Ey": [1, 0, 1],
    "Ez": [1, 1, 0],
    "AvgElec0": [0, 1, 1],
    "AvgElec1": [1, 0, 1],
    "AvgElec2": [1, 1, 0],
}


class EnzoFieldInfo(FieldInfoContainer):
    known_other_fields: KnownFieldsT = (
        ("Cooling_Time", ("s", ["cooling_time"], None)),
        ("Dengo_Cooling_Rate", ("erg/g/s", [], None)),
        ("Grackle_Cooling_Rate", ("erg/s/cm**3", [], None)),
        ("HI_kph", ("1/code_time", ["H_p0_ionization_rate"], None)),
        ("HeI_kph", ("1/code_time", ["He_p0_ionization_rate"], None)),
        ("HeII_kph", ("1/code_time", ["He_p1_ionization_rate"], None)),
        ("H2I_kdiss", ("1/code_time", ["H2_p0_dissociation_rate"], None)),
        ("HM_kph", ("1/code_time", ["H_m1_ionization_rate"], None)),
        ("H2II_kdiss", ("1/code_time", ["H2_p1_dissociation_rate"], None)),
        ("Bx", (b_units, [], None)),
        ("By", (b_units, [], None)),
        ("Bz", (b_units, [], None)),
        ("BxF", (b_units, [], None)),
        ("ByF", (b_units, [], None)),
        ("BzF", (b_units, [], None)),
        ("Ex", (e_units, [], None)),
        ("Ey", (e_units, [], None)),
        ("Ez", (e_units, [], None)),
        ("AvgElec0", (e_units, [], None)),
        ("AvgElec1", (e_units, [], None)),
        ("AvgElec2", (e_units, [], None)),
        ("RadAccel1", (ra_units, ["radiation_acceleration_x"], None)),
        ("RadAccel2", (ra_units, ["radiation_acceleration_y"], None)),
        ("RadAccel3", (ra_units, ["radiation_acceleration_z"], None)),
        ("Dark_Matter_Density", (rho_units, ["dark_matter_density"], None)),
        ("Temperature", ("K", ["temperature"], None)),
        ("Dust_Temperature", ("K", ["dust_temperature"], None)),
        ("x-velocity", (vel_units, ["velocity_x"], None)),
        ("y-velocity", (vel_units, ["velocity_y"], None)),
        ("z-velocity", (vel_units, ["velocity_z"], None)),
        ("RaySegments", ("", ["ray_segments"], None)),
        ("PhotoGamma", ("eV/code_time", ["photo_gamma"], None)),
        ("PotentialField", ("code_velocity**2", ["gravitational_potential"], None)),
        ("Density", (rho_units, ["density"], None)),
        ("Metal_Density", (rho_units, ["metal_density"], None)),
        ("SN_Colour", (rho_units, [], None)),
        # Note: we do not alias Electron_Density to anything
        ("Electron_Density", (rho_units, [], None)),
    )

    known_particle_fields: KnownFieldsT = (
        ("particle_position_x", ("code_length", [], None)),
        ("particle_position_y", ("code_length", [], None)),
        ("particle_position_z", ("code_length", [], None)),
        ("particle_velocity_x", (vel_units, [], None)),
        ("particle_velocity_y", (vel_units, [], None)),
        ("particle_velocity_z", (vel_units, [], None)),
        ("creation_time", ("code_time", [], None)),
        ("dynamical_time", ("code_time", [], None)),
        ("metallicity_fraction", ("code_metallicity", [], None)),
        ("metallicity", ("", [], None)),
        ("particle_type", ("", [], None)),
        ("particle_index", ("", [], None)),
        ("particle_mass", ("code_mass", [], None)),
        ("GridID", ("", [], None)),
        ("identifier", ("", ["particle_index"], None)),
        ("level", ("", [], None)),
        ("AccretionRate", ("code_mass/code_time", [], None)),
        ("AccretionRateTime", ("code_time", [], None)),
        ("AccretionRadius", ("code_length", [], None)),
        ("RadiationLifetime", ("code_time", [], None)),
    )

    def __init__(self, ds, field_list):
        hydro_method = ds.parameters.get("HydroMethod", None)
        if hydro_method is None:
            hydro_method = ds.parameters["Physics"]["Hydro"]["HydroMethod"]
        if hydro_method == 2:
            sl_left = slice(None, -2, None)
            sl_right = slice(1, -1, None)
            div_fac = 1.0
        else:
            sl_left = slice(None, -2, None)
            sl_right = slice(2, None, None)
            div_fac = 2.0
        slice_info = (sl_left, sl_right, div_fac)
        super().__init__(ds, field_list, slice_info)

        # setup nodal flag information
        for field in NODAL_FLAGS:
            if ("enzo", field) in self:
                finfo = self["enzo", field]
                finfo.nodal_flag = np.array(NODAL_FLAGS[field])

    def add_species_field(self, species):
        # This is currently specific to Enzo.  Hopefully in the future we will
        # have deeper integration with other systems, such as Dengo, to provide
        # better understanding of ionization and molecular states.
        #
        # We have several fields to add based on a given species field.  First
        # off, we add the species field itself.  Then we'll add a few more
        # items...
        #
        self.add_output_field(
            ("enzo", f"{species}_Density"),
            sampling_type="cell",
            take_log=True,
            units="code_mass/code_length**3",
        )
        yt_name = known_species_names[species]
        # don't alias electron density since mass is wrong
        if species != "Electron":
            self.alias(("gas", f"{yt_name}_density"), ("enzo", f"{species}_Density"))

    def setup_species_fields(self):
        species_names = [
            fn.rsplit("_Density")[0]
            for ft, fn in self.field_list
            if fn.endswith("_Density")
        ]
        species_names = [sp for sp in species_names if sp in known_species_names]

        def _electron_density(field, data):
            return data[("enzo", "Electron_Density")] * (me / mp)

        self.add_field(
            ("gas", "El_density"),
            sampling_type="cell",
            function=_electron_density,
            units=self.ds.unit_system["density"],
        )
        for sp in species_names:
            self.add_species_field(sp)
            self.species_names.append(known_species_names[sp])
        self.species_names.sort()  # bb #1059

    def setup_fluid_fields(self):
        from yt.fields.magnetic_field import setup_magnetic_field_aliases

        # Now we conditionally load a few other things.
        params = self.ds.parameters
        multi_species = params.get("MultiSpecies", None)
        dengo = params.get("DengoChemistryModel", 0)
        if multi_species is None:
            multi_species = params["Physics"]["AtomicPhysics"]["MultiSpecies"]
        if multi_species > 0 or dengo == 1:
            self.setup_species_fields()
        self.setup_energy_field()
        setup_magnetic_field_aliases(self, "enzo", [f"B{ax}" for ax in "xyz"])

    def setup_energy_field(self):
        unit_system = self.ds.unit_system
        # We check which type of field we need, and then we add it.
        ge_name = None
        te_name = None
        params = self.ds.parameters
        multi_species = params.get("MultiSpecies", None)
        if multi_species is None:
            multi_species = params["Physics"]["AtomicPhysics"]["MultiSpecies"]
        hydro_method = params.get("HydroMethod", None)
        if hydro_method is None:
            hydro_method = params["Physics"]["Hydro"]["HydroMethod"]
        dual_energy = params.get("DualEnergyFormalism", None)
        if dual_energy is None:
            dual_energy = params["Physics"]["Hydro"]["DualEnergyFormalism"]
        if ("enzo", "Gas_Energy") in self.field_list:
            ge_name = "Gas_Energy"
        elif ("enzo", "GasEnergy") in self.field_list:
            ge_name = "GasEnergy"
        if ("enzo", "Total_Energy") in self.field_list:
            te_name = "Total_Energy"
        elif ("enzo", "TotalEnergy") in self.field_list:
            te_name = "TotalEnergy"

        if hydro_method == 2 and te_name is not None:
            self.add_output_field(
                ("enzo", te_name), sampling_type="cell", units="code_velocity**2"
            )
            self.alias(("gas", "specific_thermal_energy"), ("enzo", te_name))

            def _ge_plus_kin(field, data):
                ret = data[("enzo", te_name)] + 0.5 * data[("gas", "velocity_x")] ** 2.0
                if data.ds.dimensionality > 1:
                    ret += 0.5 * data[("gas", "velocity_y")] ** 2.0
                if data.ds.dimensionality > 2:
                    ret += 0.5 * data[("gas", "velocity_z")] ** 2.0
                return ret

            self.add_field(
                ("gas", "specific_total_energy"),
                sampling_type="cell",
                function=_ge_plus_kin,
                units=unit_system["specific_energy"],
            )
        elif dual_energy == 1:
            if te_name is not None:
                self.add_output_field(
                    ("enzo", te_name), sampling_type="cell", units="code_velocity**2"
                )
                self.alias(
                    ("gas", "specific_total_energy"),
                    ("enzo", te_name),
                    units=unit_system["specific_energy"],
                )
            if ge_name is not None:
                self.add_output_field(
                    ("enzo", ge_name), sampling_type="cell", units="code_velocity**2"
                )
                self.alias(
                    ("gas", "specific_thermal_energy"),
                    ("enzo", ge_name),
                    units=unit_system["specific_energy"],
                )
        elif hydro_method in (4, 6) and te_name is not None:
            self.add_output_field(
                ("enzo", te_name), sampling_type="cell", units="code_velocity**2"
            )

            # Subtract off B-field energy
            def _sub_b(field, data):
                ret = data[("enzo", te_name)] - 0.5 * data[("gas", "velocity_x")] ** 2.0
                if data.ds.dimensionality > 1:
                    ret -= 0.5 * data[("gas", "velocity_y")] ** 2.0
                if data.ds.dimensionality > 2:
                    ret -= 0.5 * data[("gas", "velocity_z")] ** 2.0
                ret -= (
                    data[("gas", "magnetic_energy_density")] / data[("gas", "density")]
                )
                return ret

            self.add_field(
                ("gas", "specific_thermal_energy"),
                sampling_type="cell",
                function=_sub_b,
                units=unit_system["specific_energy"],
            )
        elif te_name is not None:  # Otherwise, we assume TotalEnergy is kinetic+thermal
            self.add_output_field(
                ("enzo", te_name), sampling_type="cell", units="code_velocity**2"
            )
            self.alias(
                ("gas", "specific_total_energy"),
                ("enzo", te_name),
                units=unit_system["specific_energy"],
            )

            def _tot_minus_kin(field, data):
                ret = data[("enzo", te_name)] - 0.5 * data[("gas", "velocity_x")] ** 2.0
                if data.ds.dimensionality > 1:
                    ret -= 0.5 * data[("gas", "velocity_y")] ** 2.0
                if data.ds.dimensionality > 2:
                    ret -= 0.5 * data[("gas", "velocity_z")] ** 2.0
                return ret

            self.add_field(
                ("gas", "specific_thermal_energy"),
                sampling_type="cell",
                function=_tot_minus_kin,
                units=unit_system["specific_energy"],
            )
        if multi_species == 0 and "Mu" in params:

            def _mean_molecular_weight(field, data):
                return params["Mu"] * data["index", "ones"]

            self.add_field(
                ("gas", "mean_molecular_weight"),
                sampling_type="cell",
                function=_mean_molecular_weight,
                units="",
            )

            def _number_density(field, data):
                return data["gas", "density"] / (mp * params["Mu"])

            self.add_field(
                ("gas", "number_density"),
                sampling_type="cell",
                function=_number_density,
                units=unit_system["number_density"],
            )

    def setup_particle_fields(self, ptype):
        def _age(field, data):
            return data.ds.current_time - data[("all", "creation_time")]

        self.add_field(
            (ptype, "age"), sampling_type="particle", function=_age, units="yr"
        )

        super().setup_particle_fields(ptype)
