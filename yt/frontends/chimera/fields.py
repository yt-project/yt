"""
Chimera-specific fields



"""
from yt._typing import KnownFieldsT
from yt.fields.field_info_container import FieldInfoContainer

# We need to specify which fields we might have in our dataset.  The field info
# container subclass here will define which fields it knows about.  There are
# optionally methods on it that get called which can be subclassed.


class ChimeraFieldInfo(FieldInfoContainer):
    known_other_fields: KnownFieldsT = (
        ("e_int", ("erg", ["Internal Energy"], "Internal Energy")),
        ("entropy", ("", ["Entropy"], None)),
        ("rho_c", ("g/cm**3", ["density", "Density"], "Density")),
        ("dudt_nu", ("erg/s", [], None)),
        ("dudt_nuc", ("erg/s", [], None)),
        ("grav_x_c", ("cm/s**2", [], None)),
        ("grav_y_c", ("cm/s**2", [], None)),
        ("grav_z_c", ("cm/s**2", [], None)),
        ("press", ("erg/cm**3", ["pressure"], "Pressure")),
        ("t_c", ("K", ["temperature"], "Temperature")),
        ("u_c", ("cm/s", ["v_radial"], "Radial Velocity")),
        ("v_c", ("cm/s", ["v_theta"], "Theta Velocity")),
        ("v_csound", ("", [], None)),
        ("wBVMD", ("1/s", [], "BruntViasala_freq")),
        ("w_c", ("cm/s", ["v_phi"], "Phi Velocity")),
        ("ye_c", ("", [], None)),
        ("ylep", ("", [], None)),
        ("a_nuc_rep_c", ("", [], None)),
        ("be_nuc_rep_c", ("", [], None)),
        ("e_book", ("", [], None)),
        ("nse_c", ("", [], None)),
        ("z_nuc_rep_c", ("", [], None)),
    )
    # Each entry here is of the form
    # ( "name", ("units", ["fields", "to", "alias"], # "display_name")),

    known_particle_fields = (
        # Identical form to above
    )

    def __init__(self, ds, field_list):
        super().__init__(ds, field_list)
        # If you want, you can check self.field_list

    def setup_fluid_fields(self):
        # Here we do anything that might need info about the dataset.
        # You can use self.alias, self.add_output_field (for on-disk fields)
        # and self.add_field (for derived fields).

        def _test(field, data):

            return data["chimera", "rho_c"]

        self.add_field(
            ("chimera", "test"), sampling_type="cell", function=_test, units="g/cm**3"
        )

    def setup_particle_fields(self, ptype):
        super().setup_particle_fields(ptype)
        # This will get called for every particle type.
