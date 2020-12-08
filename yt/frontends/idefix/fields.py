from yt.fields.field_info_container import FieldInfoContainer

# We need to specify which fields we might have in our dataset.  The field info
# container subclass here will define which fields it knows about.  There are
# optionally methods on it that get called which can be subclassed.


class IdefixFieldInfo(FieldInfoContainer):
    known_other_fields = (
        ("Vc-RHO", ("code_mass / code_length**3", ["density"], None)),
        ("Vc-VX1", ("code_length / code_time", ["velocity_x"], None)),
        ("Vc-VX2", ("code_length / code_time", ["velocity_y"], None)),
        ("Vc-VX3", ("code_length / code_time", ["velocity_z"], None)),
    )
    # note that velocity '_x', '_y' and '_z' aliases are meant to be
    # overwriten according to geometry in self.setup_fluid_aliases

    known_particle_fields = (
        # Identical form to above
        # ( "name", ("units", ["fields", "to", "alias"], "display_name")),
    )

    def __init__(self, ds, field_list):
        super(IdefixFieldInfo, self).__init__(ds, field_list)
        # If you want, you can check self.field_list

    def setup_fluid_fields(self):
        # Here we do anything that might need info about the dataset.
        # You can use self.alias, self.add_output_field (for on-disk fields)
        # and self.add_field (for derived fields).
        pass

    def setup_particle_fields(self, ptype):
        super(IdefixFieldInfo, self).setup_particle_fields(ptype)
        # This will get called for every particle type.
