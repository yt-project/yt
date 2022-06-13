"""
CF-radial-specific fields



"""
from typing import List

from yt.fields.field_info_container import FieldInfoContainer


class CFRadialFieldInfo(FieldInfoContainer):
    known_other_fields = ()  # fields are set dynamically
    known_particle_fields = ()
    units_to_ignore = ("dBz", "dBZ", "ratio")  # set as nondimensional if found
    field_units_ignored: List[str] = []  # fields for which units have been ignored

    # (find, replace) pairs for sanitizing:
    unit_subs = (("degrees", "degree"), ("meters", "m"), ("_per_", "/"))

    def setup_fluid_fields(self):
        # Here we dynamically add fields available in our netcdf file for to the
        # FieldInfoContainer with sanitized units.

        for field in self.field_list:  # field here is ('fluid_type', 'field') tuple
            units = self.ds.field_units.get(field, "")

            # sanitization of the units
            if units in self.units_to_ignore:
                self.field_units_ignored.append(field)
                units = ""

            for findstr, repstr in self.unit_subs:
                units = units.replace(findstr, repstr)

            self.add_output_field(field, "cell", units=units)
