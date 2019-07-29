from yt.frontends.gadget.api import GadgetFieldInfo
from yt.fields.magnetic_field import \
    setup_magnetic_field_aliases

metal_elements = ["He", "C", "N", "O", "Ne",
                  "Mg", "Si", "Fe"]

class ArepoFieldInfo(GadgetFieldInfo):
    known_particle_fields = GadgetFieldInfo.known_particle_fields + \
                            (("smoothing_length", ("code_length", [], None)),
                             ("MagneticField",
                              ("code_magnetic", ["particle_magnetic_field"], None)),
                             ("GFM_Metallicity", ("", ["metallicity"], None)),
                             ("GFM_Metals_00", ("", ["H_fraction"], None)),
                             ("GFM_Metals_01", ("", ["He_fraction"], None)),
                             ("GFM_Metals_02", ("", ["C_fraction"], None)),
                             ("GFM_Metals_03", ("", ["N_fraction"], None)),
                             ("GFM_Metals_04", ("", ["O_fraction"], None)),
                             ("GFM_Metals_05", ("", ["Ne_fraction"], None)),
                             ("GFM_Metals_06", ("", ["Mg_fraction"], None)),
                             ("GFM_Metals_07", ("", ["Si_fraction"], None)),
                             ("GFM_Metals_08", ("", ["Fe_fraction"], None)),
                             )

    def setup_gas_particle_fields(self, ptype):
        super(ArepoFieldInfo, self).setup_gas_particle_fields(ptype)
        if ("PartType0", "GFM_Metals_00") in self.field_list:
            self.nuclei_names = metal_elements
            self.species_names = ["H"] + metal_elements

        magnetic_field = "MagneticField"
        if (ptype, magnetic_field) in self.field_list:
            setup_magnetic_field_aliases(
                self, ptype, magnetic_field
            )

