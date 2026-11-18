import numpy as np

from yt._typing import KnownFieldsT
from yt.fields.field_info_container import FieldInfoContainer
from yt.units import dimensions
from yt.utilities.physical_constants import kboltz, mh

b_units = "code_magnetic"
pres_units = "code_mass / (code_length * code_time**2)"
rho_units = "code_mass / code_length**3"
vel_units = "code_length / code_time"
mom_units = "code_mass / (code_length**2 * code_time)"
pot_units = "code_length**2 / code_time**2"


def velocity_field(j):
    def _velocity(data):
        return data["pasta", f"fluid_mom_{j}"] / data["gas", "density"]

    return _velocity


class PastaFieldInfo(FieldInfoContainer):
    known_other_fields: KnownFieldsT = (
        ("fluid_den_prim", (rho_units, ["density"], None)),
        ("fluid_den_cons", (rho_units, ["density"], None)),
        ("fluid_vel_1", (vel_units, [], None)),
        ("fluid_vel_2", (vel_units, [], None)),
        ("fluid_vel_3", (vel_units, [], None)),
        ("fluid_mom_1", (mom_units, [], None)),
        ("fluid_mom_2", (mom_units, [], None)),
        ("fluid_mom_3", (mom_units, [], None)),
        ("fluid_pressure", (pres_units, ["pressure"], None)),
        ("fluid_internale", (pres_units, ["thermal_energy_density"], None)),
        ("fluid_internal_energy", (pres_units, ["thermal_energy_density"], None)),
        (
            "fluid_specific_internal_energy",
            (pot_units, ["specific_thermal_energy"], None),
        ),
        ("fluid_energy", (pres_units, ["total_energy_per_volume"], None)),
        ("temperature", ("K", [], None)),
        ("field_bcc_1", (b_units, [], None)),
        ("field_bcc_2", (b_units, [], None)),
        ("field_bcc_3", (b_units, [], None)),
        ("field_bf_1", (b_units, [], None)),
        ("field_bf_2", (b_units, [], None)),
        ("field_bf_3", (b_units, [], None)),
        ("field_emf_1", (f"{b_units} * {vel_units}", [], None)),
        ("field_emf_2", (f"{b_units} * {vel_units}", [], None)),
        ("field_emf_3", (f"{b_units} * {vel_units}", [], None)),
        ("gravity_phi", (pot_units, ["gravitational_potential"], None)),
        ("gravity_phi_old", (pot_units, [], None)),
        ("particle_rho", (rho_units, ["dark_matter_density"], None)),
        ("scalar_cons_0", (rho_units, [], None)),
        ("scalar_cons_1", (rho_units, [], None)),
        ("scalar_cons_2", (rho_units, [], None)),
        ("scalar_cons_3", (rho_units, [], None)),
        ("scalar_cons_4", (rho_units, [], None)),
        ("scalar_cons_5", (rho_units, [], None)),
        ("scalar_cons_6", (rho_units, [], None)),
        ("scalar_cons_7", (rho_units, [], None)),
        ("scalar_cons_8", (rho_units, [], None)),
        ("scalar_cons_9", (rho_units, [], None)),
        ("scalar_cons_10", (rho_units, [], None)),
        ("scalar_cons_11", (rho_units, [], None)),
        ("scalar_cons_12", (rho_units, [], None)),
        ("scalar_prim_0", ("", [], None)),
        ("scalar_prim_1", ("", [], None)),
        ("scalar_prim_2", ("", [], None)),
        ("scalar_prim_3", ("", [], None)),
        ("scalar_prim_4", ("", [], None)),
        ("scalar_prim_5", ("", [], None)),
        ("scalar_prim_6", ("", [], None)),
        ("scalar_prim_7", ("", [], None)),
        ("scalar_prim_8", ("", [], None)),
        ("scalar_prim_9", ("", [], None)),
        ("scalar_prim_10", ("", [], None)),
        ("scalar_prim_11", ("", [], None)),
        ("scalar_prim_12", ("", [], None)),
    )

    known_particle_fields: KnownFieldsT = (
        ("particle_position_x", ("code_length", [], None)),
        ("particle_position_y", ("code_length", [], None)),
        ("particle_position_z", ("code_length", [], None)),
        ("particle_velocity_x", (vel_units, [], None)),
        ("particle_velocity_y", (vel_units, [], None)),
        ("particle_velocity_z", (vel_units, [], None)),
        ("particle_mass", ("code_mass", [], None)),
        ("particle_index", ("", [], None)),
        ("x1", ("code_length", [], None)),
        ("x2", ("code_length", [], None)),
        ("x3", ("code_length", [], None)),
        ("v1", (vel_units, [], None)),
        ("v2", (vel_units, [], None)),
        ("v3", (vel_units, [], None)),
        ("mass", ("code_mass", [], None)),
        ("pid", ("", [], None)),
        ("blk", ("", [], None)),
        ("must_refine", ("", [], None)),
    )

    def setup_fluid_fields(self):
        from yt.fields.magnetic_field import setup_magnetic_field_aliases

        unit_system = self.ds.unit_system

        # Velocity and momentum aliasing based on active axis ordering
        for i, comp in enumerate(self.ds.coordinates.axis_order):
            vel_field = ("pasta", f"fluid_vel_{i + 1}")
            mom_field = ("pasta", f"fluid_mom_{i + 1}")
            if vel_field in self.field_list:
                self.add_output_field(vel_field, sampling_type="cell", units=vel_units)
                self.alias(
                    ("gas", f"velocity_{comp}"),
                    vel_field,
                    units=unit_system["velocity"],
                )
            elif mom_field in self.field_list:
                self.add_output_field(mom_field, sampling_type="cell", units=mom_units)
                self.add_field(
                    ("gas", f"velocity_{comp}"),
                    sampling_type="cell",
                    function=velocity_field(i + 1),
                    units=unit_system["velocity"],
                )

            if mom_field in self.field_list:
                self.add_output_field(mom_field, sampling_type="cell", units=mom_units)
                self.alias(
                    ("gas", f"momentum_density_{comp}"),
                    mom_field,
                    units=unit_system["density"] * unit_system["velocity"],
                )

        # Thermal energy and pressure
        if ("pasta", "fluid_pressure") in self.field_list:

            def _specific_thermal_energy(data):
                return (
                    data["gas", "pressure"]
                    / (data.ds.gamma - 1.0)
                    / data["gas", "density"]
                )

            self.add_field(
                ("gas", "specific_thermal_energy"),
                sampling_type="cell",
                function=_specific_thermal_energy,
                units=unit_system["specific_energy"],
            )
        elif ("pasta", "fluid_specific_internal_energy") in self.field_list:

            def _pressure(data):
                return (
                    (data.ds.gamma - 1.0)
                    * data["gas", "density"]
                    * data["gas", "specific_thermal_energy"]
                )

            self.add_field(
                ("gas", "pressure"),
                sampling_type="cell",
                function=_pressure,
                units=unit_system["pressure"],
            )
        elif ("pasta", "fluid_internale") in self.field_list or (
            "pasta",
            "fluid_internal_energy",
        ) in self.field_list:
            int_e_field = (
                ("pasta", "fluid_internal_energy")
                if ("pasta", "fluid_internal_energy") in self.field_list
                else ("pasta", "fluid_internale")
            )

            def _specific_thermal_energy(data):
                return data[int_e_field] / data["gas", "density"]

            self.add_field(
                ("gas", "specific_thermal_energy"),
                sampling_type="cell",
                function=_specific_thermal_energy,
                units=unit_system["specific_energy"],
            )

            def _pressure(data):
                return (data.ds.gamma - 1.0) * data[int_e_field]

            self.add_field(
                ("gas", "pressure"),
                sampling_type="cell",
                function=_pressure,
                units=unit_system["pressure"],
            )
        elif ("pasta", "fluid_energy") in self.field_list:

            def eint_from_etot(data):
                eint = (
                    data["pasta", "fluid_energy"]
                    - data["gas", "kinetic_energy_density"]
                )
                if ("gas", "magnetic_energy_density") in self:
                    eint = eint - data["gas", "magnetic_energy_density"]
                return np.maximum(eint, 0.0)

            def _specific_thermal_energy(data):
                return eint_from_etot(data) / data["gas", "density"]

            self.add_field(
                ("gas", "specific_thermal_energy"),
                sampling_type="cell",
                function=_specific_thermal_energy,
                units=unit_system["specific_energy"],
            )

            if ("pasta", "fluid_pressure") not in self.field_list:

                def _pressure(data):
                    return (data.ds.gamma - 1.0) * eint_from_etot(data)

                self.add_field(
                    ("gas", "pressure"),
                    sampling_type="cell",
                    function=_pressure,
                    units=unit_system["pressure"],
                )

        has_thermal = any(
            ("pasta", f) in self.field_list
            for f in (
                "fluid_pressure",
                "fluid_internale",
                "fluid_internal_energy",
                "fluid_specific_internal_energy",
                "fluid_energy",
            )
        )

        # Sound speed
        if has_thermal and ("gas", "sound_speed") not in self:

            def _sound_speed(data):
                return np.sqrt(
                    data.ds.gamma * data["gas", "pressure"] / data["gas", "density"]
                )

            self.add_field(
                ("gas", "sound_speed"),
                sampling_type="cell",
                function=_sound_speed,
                units=unit_system["velocity"],
            )

        # Temperature
        if ("pasta", "temperature") in self.field_list:
            self.alias(
                ("gas", "temperature"),
                ("pasta", "temperature"),
                units=unit_system["temperature"],
            )
        elif has_thermal:

            def _temperature(data):
                return (
                    (data["gas", "pressure"] / data["gas", "density"])
                    * data.ds.mu
                    * mh
                    / kboltz
                )

            self.add_field(
                ("gas", "temperature"),
                sampling_type="cell",
                function=_temperature,
                units=unit_system["temperature"],
            )

        # Magnetic fields
        b_fields = [f"field_bcc_{i + 1}" for i in range(3)]
        bf_fields = [f"field_bf_{i + 1}" for i in range(3)]
        if any(("pasta", bf) in self.field_list for bf in b_fields):
            setup_magnetic_field_aliases(self, "pasta", b_fields)
        elif any(("pasta", bf) in self.field_list for bf in bf_fields):
            setup_magnetic_field_aliases(self, "pasta", bf_fields)

        for i, comp in enumerate(self.ds.coordinates.axis_order):
            if ("pasta", f"field_bf_{i + 1}") in self.field_list:
                b_target_units = (
                    unit_system["magnetic_field_cgs"]
                    if unit_system.base_units[dimensions.current_mks] is None
                    else unit_system["magnetic_field_mks"]
                )
                self.alias(
                    ("gas", f"magnetic_field_face_{comp}"),
                    ("pasta", f"field_bf_{i + 1}"),
                    units=b_target_units,
                )

        # Species fields
        prim_chem = self.ds.parameters.get("chemistry/primordial_chemistry")
        if prim_chem is None:
            prim_chem = self.ds.parameters.get("primordial_chemistry")

        scalar_cons_indices = set()
        scalar_prim_indices = set()
        for ft, fn in self.field_list:
            if ft == "pasta":
                if fn.startswith("scalar_cons_"):
                    try:
                        scalar_cons_indices.add(int(fn.split("scalar_cons_")[1]))
                    except ValueError:
                        pass
                elif fn.startswith("scalar_prim_"):
                    try:
                        scalar_prim_indices.add(int(fn.split("scalar_prim_")[1]))
                    except ValueError:
                        pass

        all_scalar_indices = scalar_cons_indices | scalar_prim_indices

        if prim_chem is None and all_scalar_indices:
            max_idx = max(all_scalar_indices)
            if max_idx >= 11:
                prim_chem = 3
            elif max_idx >= 8:
                prim_chem = 2
            elif max_idx >= 5 or 0 in all_scalar_indices:
                prim_chem = 1
            else:
                prim_chem = 0

        self.ds.primordial_chemistry = prim_chem

        species_map = []
        metal_idx = None
        if prim_chem == 1:
            species_map = [
                (0, "H_p0"),
                (1, "H_p1"),
                (2, "He_p0"),
                (3, "He_p1"),
                (4, "He_p2"),
                (5, "El"),
            ]
            metal_idx = 6
        elif prim_chem == 2:
            species_map = [
                (0, "H_p0"),
                (1, "H_p1"),
                (2, "He_p0"),
                (3, "He_p1"),
                (4, "He_p2"),
                (5, "El"),
                (6, "H_m1"),
                (7, "H2_p0"),
                (8, "H2_p1"),
            ]
            metal_idx = 9
        elif prim_chem == 3:
            species_map = [
                (0, "H_p0"),
                (1, "H_p1"),
                (2, "He_p0"),
                (3, "He_p1"),
                (4, "He_p2"),
                (5, "El"),
                (6, "H_m1"),
                (7, "H2_p0"),
                (8, "H2_p1"),
                (9, "D_p0"),
                (10, "D_p1"),
                (11, "HD_p0"),
            ]
            metal_idx = 12

        species_found = []
        for idx, sp in species_map:
            has_cons = idx in scalar_cons_indices
            has_prim = idx in scalar_prim_indices
            if not (has_cons or has_prim):
                continue

            species_found.append(sp)

            if sp == "El":
                from yt.utilities.physical_constants import me, mp

                if has_cons:

                    def _el_density(data, i=idx):
                        return data["pasta", f"scalar_cons_{i}"] * (me / mp)

                    self.add_field(
                        ("gas", "El_density"),
                        sampling_type="cell",
                        function=_el_density,
                        units=unit_system["density"],
                    )
                elif has_prim:

                    def _el_density(data, i=idx):
                        return (
                            data["pasta", f"scalar_prim_{i}"]
                            * data["gas", "density"]
                            * (me / mp)
                        )

                    self.add_field(
                        ("gas", "El_density"),
                        sampling_type="cell",
                        function=_el_density,
                        units=unit_system["density"],
                    )
            else:
                if has_cons:
                    self.alias(
                        ("gas", f"{sp}_density"),
                        ("pasta", f"scalar_cons_{idx}"),
                        units=unit_system["density"],
                    )
                if has_prim:
                    self.alias(
                        ("gas", f"{sp}_fraction"),
                        ("pasta", f"scalar_prim_{idx}"),
                        units="",
                    )

        self.species_names = species_found

        if metal_idx is not None:
            has_cons_metal = metal_idx in scalar_cons_indices
            has_prim_metal = metal_idx in scalar_prim_indices
            if has_cons_metal:
                self.alias(
                    ("gas", "metal_density"),
                    ("pasta", f"scalar_cons_{metal_idx}"),
                    units=unit_system["density"],
                )
            if has_prim_metal:
                self.alias(
                    ("gas", "metallicity"),
                    ("pasta", f"scalar_prim_{metal_idx}"),
                    units="",
                )
                self.alias(
                    ("gas", "metal_fraction"),
                    ("pasta", f"scalar_prim_{metal_idx}"),
                    units="",
                )
            if has_cons_metal and not has_prim_metal:

                def _metallicity(data):
                    return data["gas", "metal_density"] / data["gas", "density"]

                self.add_field(
                    ("gas", "metallicity"),
                    sampling_type="cell",
                    function=_metallicity,
                    units="",
                )
                self.alias(("gas", "metal_fraction"), ("gas", "metallicity"))
            elif has_prim_metal and not has_cons_metal:

                def _metal_density(data):
                    return data["gas", "metallicity"] * data["gas", "density"]

                self.add_field(
                    ("gas", "metal_density"),
                    sampling_type="cell",
                    function=_metal_density,
                    units=unit_system["density"],
                )

    def setup_particle_fields(self, ptype):
        super().setup_particle_fields(ptype)
        unit_system = self.ds.unit_system

        if (ptype, "mass") in self.field_list and (ptype, "particle_mass") not in self:
            self.alias(
                (ptype, "particle_mass"), (ptype, "mass"), units=unit_system["mass"]
            )
        if (ptype, "pid") in self.field_list and (ptype, "particle_index") not in self:
            self.alias((ptype, "particle_index"), (ptype, "pid"), units="")

        if self.ds.geometry.name == "cylindrical":
            for ax, raw in [("r", "x1"), ("theta", "x2"), ("z", "x3")]:
                if (ptype, raw) in self.field_list and (
                    ptype,
                    f"particle_position_{ax}",
                ) not in self:
                    self.alias(
                        (ptype, f"particle_position_{ax}"),
                        (ptype, raw),
                        units=unit_system["length"],
                    )
            for ax, raw in [("r", "v1"), ("theta", "v2"), ("z", "v3")]:
                if (ptype, raw) in self.field_list and (
                    ptype,
                    f"particle_velocity_{ax}",
                ) not in self:
                    self.alias(
                        (ptype, f"particle_velocity_{ax}"),
                        (ptype, raw),
                        units=unit_system["velocity"],
                    )
        elif self.ds.geometry.name == "spherical":
            for ax, raw in [("r", "x1"), ("theta", "x2"), ("phi", "x3")]:
                if (ptype, raw) in self.field_list and (
                    ptype,
                    f"particle_position_{ax}",
                ) not in self:
                    self.alias(
                        (ptype, f"particle_position_{ax}"),
                        (ptype, raw),
                        units=unit_system["length"],
                    )
            for ax, raw in [("r", "v1"), ("theta", "v2"), ("phi", "v3")]:
                if (ptype, raw) in self.field_list and (
                    ptype,
                    f"particle_velocity_{ax}",
                ) not in self:
                    self.alias(
                        (ptype, f"particle_velocity_{ax}"),
                        (ptype, raw),
                        units=unit_system["velocity"],
                    )
