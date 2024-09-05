from yt._typing import KnownFieldsT
from yt.fields.field_info_container import FieldInfoContainer
from yt.utilities.physical_constants import kboltz, mh

b_units = "code_magnetic"
pres_units = "code_mass/(code_length*code_time**2)"
rho_units = "code_mass / code_length**3"
vel_units = "code_length / code_time"


def velocity_field(j):
    def _velocity(field, data):
        return data["athena_pp", f"mom{j}"] / data["athena_pp", "dens"]

    return _velocity


class AthenaPPFieldInfo(FieldInfoContainer):
    known_other_fields: KnownFieldsT = (
        ("rho", (rho_units, ["density"], None)),
        ("dens", (rho_units, ["density"], None)),
        ("Bcc1", (b_units, [], None)),
        ("Bcc2", (b_units, [], None)),
        ("Bcc3", (b_units, [], None)),
    )

    def setup_fluid_fields(self):
        from yt.fields.magnetic_field import setup_magnetic_field_aliases

        unit_system = self.ds.unit_system
        # Add velocity fields
        vel_prefix = "velocity"
        for i, comp in enumerate(self.ds.coordinates.axis_order):
            vel_field = ("athena_pp", "vel%d" % (i + 1))
            mom_field = ("athena_pp", "mom%d" % (i + 1))
            if vel_field in self.field_list:
                self.add_output_field(
                    vel_field, sampling_type="cell", units="code_length/code_time"
                )
                self.alias(
                    ("gas", f"{vel_prefix}_{comp}"),
                    vel_field,
                    units=unit_system["velocity"],
                )
            elif mom_field in self.field_list:
                self.add_output_field(
                    mom_field,
                    sampling_type="cell",
                    units="code_mass/code_time/code_length**2",
                )
                self.add_field(
                    ("gas", f"{vel_prefix}_{comp}"),
                    sampling_type="cell",
                    function=velocity_field(i + 1),
                    units=unit_system["velocity"],
                )
        # Figure out thermal energy field
        if ("athena_pp", "press") in self.field_list:
            self.add_output_field(
                ("athena_pp", "press"), sampling_type="cell", units=pres_units
            )
            self.alias(
                ("gas", "pressure"),
                ("athena_pp", "press"),
                units=unit_system["pressure"],
            )

            def _specific_thermal_energy(field, data):
                return (
                    data["athena_pp", "press"]
                    / (data.ds.gamma - 1.0)
                    / data["athena_pp", "rho"]
                )

            self.add_field(
                ("gas", "specific_thermal_energy"),
                sampling_type="cell",
                function=_specific_thermal_energy,
                units=unit_system["specific_energy"],
            )
        elif ("athena_pp", "Etot") in self.field_list:
            self.add_output_field(
                ("athena_pp", "Etot"), sampling_type="cell", units=pres_units
            )

            def _specific_thermal_energy(field, data):
                eint = data["athena_pp", "Etot"] - data["gas", "kinetic_energy_density"]
                if ("athena_pp", "B1") in self.field_list:
                    eint -= data["gas", "magnetic_energy_density"]
                return eint / data["athena_pp", "dens"]

            self.add_field(
                ("gas", "specific_thermal_energy"),
                sampling_type="cell",
                function=_specific_thermal_energy,
                units=unit_system["specific_energy"],
            )

        # Add temperature field
        def _temperature(field, data):
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

        setup_magnetic_field_aliases(
            self, "athena_pp", ["Bcc%d" % ax for ax in (1, 2, 3)]
        )
