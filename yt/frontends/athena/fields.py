from yt._typing import KnownFieldsT
from yt.fields.field_info_container import FieldInfoContainer
from yt.utilities.physical_constants import kboltz, mh

b_units = "code_magnetic"
pres_units = "code_pressure"
erg_units = "code_mass * (code_length/code_time)**2"
rho_units = "code_mass / code_length**3"


def velocity_field(comp):
    def _velocity(field, data):
        return data["athena", f"momentum_{comp}"] / data["athena", "density"]

    return _velocity


class AthenaFieldInfo(FieldInfoContainer):
    known_other_fields: KnownFieldsT = (
        ("density", ("code_mass/code_length**3", ["density"], None)),
        ("cell_centered_B_x", (b_units, [], None)),
        ("cell_centered_B_y", (b_units, [], None)),
        ("cell_centered_B_z", (b_units, [], None)),
        ("total_energy", ("code_pressure", ["total_energy_density"], None)),
        (
            "gravitational_potential",
            ("code_velocity**2", ["gravitational_potential"], None),
        ),
    )

    # In Athena, conservative or primitive variables may be written out.
    # By default, yt concerns itself with primitive variables. The following
    # field definitions allow for conversions to primitive variables in the
    # case that the file contains the conservative ones.

    def setup_fluid_fields(self):
        from yt.fields.magnetic_field import setup_magnetic_field_aliases

        unit_system = self.ds.unit_system
        # Add velocity fields
        for comp in "xyz":
            vel_field = ("athena", f"velocity_{comp}")
            mom_field = ("athena", f"momentum_{comp}")
            if vel_field in self.field_list:
                self.add_output_field(
                    vel_field, sampling_type="cell", units="code_length/code_time"
                )
                self.alias(
                    ("gas", f"velocity_{comp}"),
                    vel_field,
                    units=unit_system["velocity"],
                )
            elif mom_field in self.field_list:
                self.add_output_field(
                    mom_field,
                    sampling_type="cell",
                    units="code_mass/code_time/code_length**2",
                )
                self.alias(
                    ("gas", f"momentum_density_{comp}"),
                    mom_field,
                    units=unit_system["density"] * unit_system["velocity"],
                )
                self.add_field(
                    ("gas", f"velocity_{comp}"),
                    sampling_type="cell",
                    function=velocity_field(comp),
                    units=unit_system["velocity"],
                )

        # Add pressure, energy, and temperature fields
        def eint_from_etot(data):
            eint = (
                data["athena", "total_energy"] - data["gas", "kinetic_energy_density"]
            )
            if ("athena", "cell_centered_B_x") in self.field_list:
                eint -= data["gas", "magnetic_energy_density"]
            return eint

        def etot_from_pres(data):
            etot = data["athena", "pressure"] / (data.ds.gamma - 1.0)
            etot += data["gas", "kinetic_energy_density"]
            if ("athena", "cell_centered_B_x") in self.field_list:
                etot += data["gas", "magnetic_energy_density"]
            return etot

        if ("athena", "pressure") in self.field_list:
            self.add_output_field(
                ("athena", "pressure"), sampling_type="cell", units=pres_units
            )
            self.alias(
                ("gas", "pressure"),
                ("athena", "pressure"),
                units=unit_system["pressure"],
            )

            def _specific_thermal_energy(field, data):
                return (
                    data["athena", "pressure"]
                    / (data.ds.gamma - 1.0)
                    / data["athena", "density"]
                )

            self.add_field(
                ("gas", "specific_thermal_energy"),
                sampling_type="cell",
                function=_specific_thermal_energy,
                units=unit_system["specific_energy"],
            )

            def _specific_total_energy(field, data):
                return etot_from_pres(data) / data["athena", "density"]

            self.add_field(
                ("gas", "specific_total_energy"),
                sampling_type="cell",
                function=_specific_total_energy,
                units=unit_system["specific_energy"],
            )
        elif ("athena", "total_energy") in self.field_list:
            self.add_output_field(
                ("athena", "total_energy"), sampling_type="cell", units=pres_units
            )

            def _specific_thermal_energy(field, data):
                return eint_from_etot(data) / data["athena", "density"]

            self.add_field(
                ("gas", "specific_thermal_energy"),
                sampling_type="cell",
                function=_specific_thermal_energy,
                units=unit_system["specific_energy"],
            )

            def _specific_total_energy(field, data):
                return data["athena", "total_energy"] / data["athena", "density"]

            self.add_field(
                ("gas", "specific_total_energy"),
                sampling_type="cell",
                function=_specific_total_energy,
                units=unit_system["specific_energy"],
            )

        # Add temperature field
        def _temperature(field, data):
            return (
                data.ds.mu
                * data["gas", "pressure"]
                / data["gas", "density"]
                * mh
                / kboltz
            )

        self.add_field(
            ("gas", "temperature"),
            sampling_type="cell",
            function=_temperature,
            units=unit_system["temperature"],
        )

        setup_magnetic_field_aliases(
            self, "athena", [f"cell_centered_B_{ax}" for ax in "xyz"]
        )
