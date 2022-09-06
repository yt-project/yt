import numpy as np
from unyt import Zsun

from yt._typing import KnownFieldsT
from yt.fields.field_info_container import FieldInfoContainer
from yt.utilities.physical_constants import kboltz, mh

# Copied from Athena frontend
pres_units = "code_pressure"
erg_units = "code_mass * (code_length/code_time)**2"
rho_units = "code_mass / code_length**3"
mom_units = "code_mass / code_length**2 / code_time"


def velocity_field(comp):
    def _velocity(field, data):
        return data["cholla", f"momentum_{comp}"] / data["cholla", "density"]

    return _velocity


class ChollaFieldInfo(FieldInfoContainer):
    known_other_fields: KnownFieldsT = (
        # Each entry here is of the form
        # ( "name", ("units", ["fields", "to", "alias"], # "display_name")),
        ("density", (rho_units, ["density"], None)),
        ("momentum_x", (mom_units, ["momentum_x"], None)),
        ("momentum_y", (mom_units, ["momentum_y"], None)),
        ("momentum_z", (mom_units, ["momentum_z"], None)),
        ("Energy", ("code_pressure", ["total_energy_density"], None)),
        ("scalar0", (rho_units, [], None)),
        ("metal_density", (rho_units, ["metal_density"], None)),
    )

    known_particle_fields = ()

    # In Cholla, conservative variables are written out.

    def setup_fluid_fields(self):

        unit_system = self.ds.unit_system

        # Add velocity fields
        for comp in "xyz":
            self.add_field(
                ("gas", f"velocity_{comp}"),
                sampling_type="cell",
                function=velocity_field(comp),
                units=unit_system["velocity"],
            )

        # Add pressure field
        if ("cholla", "GasEnergy") in self.field_list:
            self.add_output_field(
                ("cholla", "GasEnergy"), sampling_type="cell", units=pres_units
            )
            self.alias(
                ("gas", "thermal_energy"),
                ("cholla", "GasEnergy"),
                units=unit_system["pressure"],
            )

            def _pressure(field, data):
                return (data.ds.gamma - 1.0) * data["cholla", "GasEnergy"]

        else:

            def _pressure(field, data):
                return (data.ds.gamma - 1.0) * (
                    data["cholla", "Energy"] - data["gas", "kinetic_energy_density"]
                )

        self.add_field(
            ("gas", "pressure"),
            sampling_type="cell",
            function=_pressure,
            units=unit_system["pressure"],
        )

        def _specific_total_energy(field, data):
            return data["cholla", "Energy"] / data["cholla", "density"]

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

        # Add color field if present (scalar0 / density)
        if ("cholla", "scalar0") in self.field_list:
            self.add_output_field(
                ("cholla", "scalar0"),
                sampling_type="cell",
                units=rho_units,
            )

            def _color(field, data):
                return data["cholla", "scalar0"] / data["cholla", "density"]

            self.add_field(
                ("cholla", "color"),
                sampling_type="cell",
                function=_color,
                units="",
            )

            self.alias(
                ("gas", "color"),
                ("cholla", "color"),
                units="",
            )

            # Using color field to define metallicity field, where a color of 1
            # indicates solar metallicity

            def _metallicity(field, data):
                # Ensuring that there are no negative metallicities
                return np.clip(data[("cholla", "color")], 0, np.inf) * Zsun

            self.add_field(
                ("cholla", "metallicity"),
                sampling_type="cell",
                function=_metallicity,
                units="Zsun",
            )

            self.alias(
                ("gas", "metallicity"),
                ("cholla", "metallicity"),
                units="Zsun",
            )

    def setup_particle_fields(self, ptype):
        super().setup_particle_fields(ptype)
        # This will get called for every particle type.
