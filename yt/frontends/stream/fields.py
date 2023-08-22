import re

from yt._typing import KnownFieldsT
from yt.fields.field_info_container import FieldInfoContainer


class StreamFieldInfo(FieldInfoContainer):
    known_other_fields: KnownFieldsT = (
        ("density", ("code_mass/code_length**3", ["density"], None)),
        (
            "dark_matter_density",
            ("code_mass/code_length**3", ["dark_matter_density"], None),
        ),
        ("number_density", ("1/code_length**3", ["number_density"], None)),
        ("pressure", ("dyne/code_length**2", ["pressure"], None)),
        ("specific_thermal_energy", ("erg / g", ["specific_thermal_energy"], None)),
        ("temperature", ("K", ["temperature"], None)),
        ("velocity_x", ("code_length/code_time", ["velocity_x"], None)),
        ("velocity_y", ("code_length/code_time", ["velocity_y"], None)),
        ("velocity_z", ("code_length/code_time", ["velocity_z"], None)),
        ("magnetic_field_x", ("gauss", [], None)),
        ("magnetic_field_y", ("gauss", [], None)),
        ("magnetic_field_z", ("gauss", [], None)),
        ("velocity_r", ("code_length/code_time", ["velocity_r"], None)),
        ("velocity_theta", ("code_length/code_time", ["velocity_theta"], None)),
        ("velocity_phi", ("code_length/code_time", ["velocity_phi"], None)),
        ("magnetic_field_r", ("gauss", [], None)),
        ("magnetic_field_theta", ("gauss", [], None)),
        ("magnetic_field_phi", ("gauss", [], None)),
        (
            "radiation_acceleration_x",
            ("code_length/code_time**2", ["radiation_acceleration_x"], None),
        ),
        (
            "radiation_acceleration_y",
            ("code_length/code_time**2", ["radiation_acceleration_y"], None),
        ),
        (
            "radiation_acceleration_z",
            ("code_length/code_time**2", ["radiation_acceleration_z"], None),
        ),
        ("metallicity", ("Zsun", ["metallicity"], None)),
        # We need to have a bunch of species fields here, too
        ("metal_density", ("code_mass/code_length**3", ["metal_density"], None)),
        ("hi_density", ("code_mass/code_length**3", ["hi_density"], None)),
        ("hii_density", ("code_mass/code_length**3", ["hii_density"], None)),
        ("h2i_density", ("code_mass/code_length**3", ["h2i_density"], None)),
        ("h2ii_density", ("code_mass/code_length**3", ["h2ii_density"], None)),
        ("h2m_density", ("code_mass/code_length**3", ["h2m_density"], None)),
        ("hei_density", ("code_mass/code_length**3", ["hei_density"], None)),
        ("heii_density", ("code_mass/code_length**3", ["heii_density"], None)),
        ("heiii_density", ("code_mass/code_length**3", ["heiii_density"], None)),
        ("hdi_density", ("code_mass/code_length**3", ["hdi_density"], None)),
        ("di_density", ("code_mass/code_length**3", ["di_density"], None)),
        ("dii_density", ("code_mass/code_length**3", ["dii_density"], None)),
    )

    known_particle_fields: KnownFieldsT = (
        ("particle_position", ("code_length", ["particle_position"], None)),
        ("particle_position_x", ("code_length", ["particle_position_x"], None)),
        ("particle_position_y", ("code_length", ["particle_position_y"], None)),
        ("particle_position_z", ("code_length", ["particle_position_z"], None)),
        ("particle_velocity", ("code_length/code_time", ["particle_velocity"], None)),
        (
            "particle_velocity_x",
            ("code_length/code_time", ["particle_velocity_x"], None),
        ),
        (
            "particle_velocity_y",
            ("code_length/code_time", ["particle_velocity_y"], None),
        ),
        (
            "particle_velocity_z",
            ("code_length/code_time", ["particle_velocity_z"], None),
        ),
        ("particle_index", ("", ["particle_index"], None)),
        (
            "particle_gas_density",
            ("code_mass/code_length**3", ["particle_gas_density"], None),
        ),
        ("particle_gas_temperature", ("K", ["particle_gas_temperature"], None)),
        ("particle_mass", ("code_mass", ["particle_mass"], None)),
        ("smoothing_length", ("code_length", ["smoothing_length"], None)),
        ("density", ("code_mass/code_length**3", ["density"], None)),
        ("temperature", ("code_temperature", ["temperature"], None)),
        ("creation_time", ("code_time", ["creation_time"], None)),
        ("age", ("code_time", [], None)),
    )

    def setup_fluid_fields(self):
        from yt.fields.magnetic_field import setup_magnetic_field_aliases
        from yt.fields.species_fields import setup_species_fields
        from yt.utilities.periodic_table import periodic_table

        # First grab all the element symbols from the periodic table
        # (this includes the electron and deuterium)
        symbols = list(periodic_table.elements_by_symbol)
        # Now add some common molecules
        symbols += ["H2", "CO"]
        species_names = []
        for field in self.ds.stream_handler.field_units:
            if field[0] in self.ds.particle_types:
                continue
            units = self.ds.stream_handler.field_units[field]
            if units != "":
                self.add_output_field(field, sampling_type="cell", units=units)
                # Check to see if this could be a species fraction field
                if field[1].endswith("_fraction"):
                    sp = field[1].rsplit("_fraction")[0]
                    parts = sp.split("_")
                    # parts is now either an element or molecule symbol
                    # by itself:
                    valid = parts[0] in symbols
                    # or it may also have an ionization state after it
                    if valid and len(parts) > 1 and parts[0] != "El":
                        # Note that this doesn't catch invalid ionization states,
                        # which would indicate more electron states empty than actually
                        # exist, but we'll leave that to the user to do correctly.
                        valid &= re.match("^[pm](0|[1-9][0-9]*)$", parts[1]) is not None
                    if valid:
                        # Add the species name to the list
                        species_names.append(sp)
                        # Alias the field
                        self.alias(("gas", field[1]), ("stream", field[1]))
        self.species_names = sorted(species_names)
        setup_magnetic_field_aliases(
            self,
            "stream",
            [f"magnetic_field_{ax}" for ax in self.ds.coordinates.axis_order],
        )
        setup_species_fields(self)

    def add_output_field(self, name, sampling_type, **kwargs):
        if name in self.ds.stream_handler.field_units:
            kwargs["units"] = self.ds.stream_handler.field_units[name]
        super().add_output_field(name, sampling_type, **kwargs)
