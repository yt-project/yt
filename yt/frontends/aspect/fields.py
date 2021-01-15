from yt.fields.field_info_container import FieldInfoContainer

# We need to specify which fields we might have in our dataset.  The field info
# container subclass here will define which fields it knows about.  There are
# optionally methods on it that get called which can be subclassed.


class ASPECTFieldInfo(FieldInfoContainer):

    known_other_fields = (
        ("T", ("K", [], None)),
        ("p", ("K", [], None)),
        ("viscosity", ("Pa", [], None)),
        ("velocity_x", ("m/s", [], None)),
        ("velocity_y", ("m/s", [], None)),
        ("velocity_z", ("m/s", [], None)),
        ("strain_rate", ("1/s", [], None)),
        ("thermal_diffusivity", ("m**2/s", [], None)),
        ("thermal_conductivity", ("W/m/K", [], None)),
        ("density", ("kg/m/m/m", [], None)),
        ("elastic_shear_modulus", ("Pa", [], None)),
        ("noninitial_plastic_strain", ("", [], None)),
        ("plastic_strain", ("", [], None)),
        ("plastic_yielding", ("", [], None)),
        ("artificial_viscosity", ("", [], None)),
        ("current_cohesions", ("", [], None)),
        ("current_friction_angles", ("", [], None)),
    )

    known_particle_fields = (
        # Identical form to above
        # ( "name", ("units", ["fields", "to", "alias"], # "display_name")),
    )

    log_fields = ["viscosity", "strain_rate"]

    def __init__(self, ds, field_list):

        # add the stress tensors to the known fields
        for xx in ["xx", "xy", "xz", "yy", "yz", "zz"]:
            self.known_other_fields += (("shear_stress_" + xx, ("Pa", [], None)),)
            self.known_other_fields += (("stress_" + xx, ("Pa", [], None)),)

        super(ASPECTFieldInfo, self).__init__(ds, field_list)
        for name in self:
            if name not in self.log_fields:
                self[name].take_log = False

    def setup_fluid_fields(self):
        # Here we do anything that might need info about the dataset.
        # You can use self.alias, self.add_output_field and self.add_field .
        pass

    def setup_particle_fields(self, ptype):
        # This will get called for every particle type.
        pass
