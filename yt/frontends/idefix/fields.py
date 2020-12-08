from yt.fields.field_info_container import FieldInfoContainer

# We need to specify which fields we might have in our dataset.  The field info
# container subclass here will define which fields it knows about.  There are
# optionally methods on it that get called which can be subclassed.


# this is copied from frontends/amrvac/fields.py
direction_aliases = {
    "cartesian": ("x", "y", "z"),
    "polar": ("r", "theta", "z"),
    "cylindrical": ("r", "z", "theta"),
    "spherical": ("r", "theta", "phi"),
}


class IdefixFieldInfo(FieldInfoContainer):
    known_other_fields = (
        # Each entry here is of the form
        # ( "name", ("units", ["fields", "to", "alias"], "display_name")),
        ("Vc-RHO", ("code_mass / code_length**3", ["density"], None)),
        ("Vc-VX1", ("code_length / code_time", ["velocity_1"], None)),
        ("Vc-VX2", ("code_length / code_time", ["velocity_2"], None)),
        ("Vc-VX3", ("code_length / code_time", ["velocity_3"], None)),
    )

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

        # setup velocity fields
        ax_aliases = direction_aliases[self.ds.geometry]
        for idir, alias in enumerate(ax_aliases, start=1):
            self.alias(
                ("gas", f"velocity_{alias}"),
                ("gas", f"velocity_{idir}"),
                units=self.ds.unit_system["velocity"],
            )

    def setup_particle_fields(self, ptype):
        super(IdefixFieldInfo, self).setup_particle_fields(ptype)
        # This will get called for every particle type.
