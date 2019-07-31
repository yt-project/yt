from yt.frontends.gadget.api import GadgetFieldInfo
from yt.fields.magnetic_field import \
    setup_magnetic_field_aliases
from yt.fields.species_fields import \
    add_species_field_by_fraction

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
            self.species_names = ["H"]
            if (ptype, "NeutralHydrogenAbundance") in self.field_list:
                self.species_names += ["H_p0", "H_p1"]
            self.species_names += metal_elements

        if (ptype, "MagneticField") in self.field_list:
            setup_magnetic_field_aliases(
                self, ptype, "MagneticField"
            )

        if (ptype, "NeutralHydrogenAbundance") in self.field_list:
            def _h_p0_fraction(field, data):
                return data[ptype, "GFM_Metals_00"] * \
                       data[ptype, "NeutralHydrogenAbundance"]

            self.add_field((ptype, "H_p0_fraction"),
                           sampling_type="particle",
                           function=_h_p0_fraction,
                           units="")

            def _h_p1_fraction(field, data):
                return data[ptype, "GFM_Metals_00"] * \
                       (1.0-data[ptype, "NeutralHydrogenAbundance"])

            self.add_field((ptype, "H_p1_fraction"),
                           sampling_type="particle",
                           function=_h_p1_fraction,
                           units="")

            add_species_field_by_fraction(self, ptype, "H_p0")
            add_species_field_by_fraction(self, ptype, "H_p1")

        if (ptype, "ElectronAbundance") in self.field_list:
            def _el_number_density(field, data):
                return data[ptype, "ElectronAbundance"] * \
                       data[ptype, "H_number_density"]
            self.add_field((ptype, "El_number_density"),
                           sampling_type="particle",
                           function=_el_number_density,
                           units=self.ds.unit_system["number_density"])
