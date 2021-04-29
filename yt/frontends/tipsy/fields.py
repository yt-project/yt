from yt.frontends.sph.fields import SPHFieldInfo


class TipsyFieldInfo(SPHFieldInfo):
    known_particle_fields = SPHFieldInfo.known_particle_fields + (
        ("smoothing_length", ("code_length", [], None)),
    )
    aux_particle_fields = {
        "uDotFB": ("uDotFB", ("code_mass * code_velocity**2", [""], None)),
        "uDotAV": ("uDotAV", ("code_mass * code_velocity**2", [""], None)),
        "uDotPdV": ("uDotPdV", ("code_mass * code_velocity**2", [""], None)),
        "uDotHydro": ("uDotHydro", ("code_mass * code_velocity**2", [""], None)),
        "uDotDiff": ("uDotDiff", ("code_mass * code_velocity**2", [""], None)),
        "uDot": ("uDot", ("code_mass * code_velocity**2", [""], None)),
        "coolontime": ("coolontime", ("code_time", [""], None)),
        "timeform": ("timeform", ("code_time", [""], None)),
        "massform": ("massform", ("code_mass", [""], None)),
        "HI": ("HI", ("dimensionless", ["H_fraction"], None)),
        "HII": ("HII", ("dimensionless", ["H_p1_fraction"], None)),
        "HeI": ("HeI", ("dimensionless", ["He_fraction"], None)),
        "HeII": ("HeII", ("dimensionless", ["He_p2_fraction"], None)),
        "OxMassFrac": ("OxMassFrac", ("dimensionless", ["O_fraction"], None)),
        "FeMassFrac": ("FeMassFrac", ("dimensionless", ["Fe_fraction"], None)),
        "c": ("c", ("code_velocity", [""], None)),
        "acc": ("acc", ("code_velocity / code_time", [""], None)),
        "accg": ("accg", ("code_velocity / code_time", [""], None)),
        "smoothlength": ("smoothlength", ("code_length", ["smoothing_length"], None)),
    }

    def __init__(self, ds, field_list, slice_info=None):
        for field in field_list:
            if (
                field[1] in self.aux_particle_fields.keys()
                and self.aux_particle_fields[field[1]] not in self.known_particle_fields
            ):
                self.known_particle_fields += (self.aux_particle_fields[field[1]],)
        super().__init__(ds, field_list, slice_info)
