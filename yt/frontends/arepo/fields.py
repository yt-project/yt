from yt.frontends.gadget.api import GadgetFieldInfo
from yt.fields.magnetic_field import \
    setup_magnetic_field_aliases


class ArepoFieldInfo(GadgetFieldInfo):
    known_particle_fields = GadgetFieldInfo.known_particle_fields + \
                            (("smoothing_length", ("code_length", [], None)),
                             ("MagneticField",
                              ("code_magnetic", ["particle_magnetic_field"], None)),
                             )

    def setup_gas_particle_fields(self, ptype):
        super(ArepoFieldInfo, self).setup_gas_particle_fields(ptype)

        magnetic_field = "MagneticField"
        if (ptype, magnetic_field) in self.field_list:
            setup_magnetic_field_aliases(
                self, ptype, magnetic_field
            )
