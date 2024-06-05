from yt.frontends.sph.fields import SPHFieldInfo


class SwiftFieldInfo(SPHFieldInfo):
    def __init__(self, ds, field_list, slice_info=None):
        self.known_particle_fields += (
            (
                "InternalEnergies",
                ("code_specific_energy", ["specific_thermal_energy"], None),
            ),
            ("Densities", ("code_mass / code_length**3", ["density"], None)),
            ("SmoothingLengths", ("code_length", ["smoothing_length"], None)),
        )
        super().__init__(ds, field_list, slice_info)

    def setup_particle_fields(self, ptype, *args, **kwargs):
        super().setup_particle_fields(ptype, *args, **kwargs)

        if ptype in ("PartType0", "Gas"):
            self.setup_gas_particle_fields(ptype)

    def setup_gas_particle_fields(self, ptype):
        self.alias((ptype, "temperature"), (ptype, "Temperatures"))
        self.alias(("gas", "temperature"), (ptype, "Temperatures"))

        for ax in ("x", "y", "z"):
            self.alias((ptype, ax), (ptype, "particle_position_" + ax))
