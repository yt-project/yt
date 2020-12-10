from yt.fields.field_info_container import FieldInfoContainer
from yt.fields.magnetic_field import setup_magnetic_field_aliases


class IdefixFieldInfo(FieldInfoContainer):
    known_other_fields = (
        ("Vc-RHO", ("code_mass / code_length**3", ["density"], None)),
        ("Vc-VX1", ("code_length / code_time", ["velocity_x"], None)),
        ("Vc-VX2", ("code_length / code_time", ["velocity_y"], None)),
        ("Vc-VX3", ("code_length / code_time", ["velocity_z"], None)),
        ("Vc-BX1", ("code_magnetic", [], None)),
        ("Vc-BX2", ("code_magnetic", [], None)),
        ("Vc-BX3", ("code_magnetic", [], None)),
        ("Vc-PRS", ("code_pressure", ["pressure"], None)),
    )
    # note that velocity '_x', '_y' and '_z' aliases are meant to be
    # overwriten according to geometry in self.setup_fluid_aliases

    known_particle_fields = ()

    def setup_fluid_fields(self):
        setup_magnetic_field_aliases(self, "idefix", [f"Vc-BX{idir}" for idir in "123"])
