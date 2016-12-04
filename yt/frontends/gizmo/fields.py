"""
Gizmo-specific fields




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, yt Development Team
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.fields.field_info_container import \
    FieldInfoContainer
from yt.fields.particle_fields import \
    add_volume_weighted_smoothed_field
from yt.fields.species_fields import \
    add_species_field_by_density, \
    setup_species_fields
from yt.frontends.gadget.fields import \
    GadgetFieldInfo
from yt.frontends.sph.fields import \
    SPHFieldInfo

metal_elements = ["He", "C", "N", "O", "Ne",
                  "Mg", "Si", "S", "Ca", "Fe"]

class GizmoFieldInfo(GadgetFieldInfo):
    known_particle_fields = (
        ("Mass", ("code_mass", ["particle_mass"], None)),
        ("Masses", ("code_mass", ["particle_mass"], None)),
        ("Coordinates", ("code_length", ["particle_position"], None)),
        ("Velocity", ("code_velocity", ["particle_velocity"], None)),
        ("Velocities", ("code_velocity", ["particle_velocity"], None)),
        ("ParticleIDs", ("", ["particle_index"], None)),
        ("InternalEnergy", ("code_velocity ** 2", ["thermal_energy"], None)),
        ("SmoothingLength", ("code_length", ["smoothing_length"], None)),
        ("Density", ("code_mass / code_length**3", ["density"], None)),
        ("MaximumTemperature", ("K", [], None)),
        ("Temperature", ("K", ["temperature"], None)),
        ("Epsilon", ("code_length", [], None)),
        ("Metals", ("code_metallicity", ["metallicity"], None)),
        ("Metallicity", ("code_metallicity", ["metallicity"], None)),
        ("Phi", ("code_length", [], None)),
        ("StarFormationRate", ("Msun / yr", [], None)),
        ("FormationTime", ("code_time", ["creation_time"], None)),
        ("Metallicity_00", ("", ["metallicity"], None)),
        ("Metallicity_01", ("", ["He_metallicity"], None)),
        ("Metallicity_02", ("", ["C_metallicity"], None)),
        ("Metallicity_03", ("", ["N_metallicity"], None)),
        ("Metallicity_04", ("", ["O_metallicity"], None)),
        ("Metallicity_05", ("", ["Ne_metallicity"], None)),
        ("Metallicity_06", ("", ["Mg_metallicity"], None)),
        ("Metallicity_07", ("", ["Si_metallicity"], None)),
        ("Metallicity_08", ("", ["S_metallicity"], None)),
        ("Metallicity_09", ("", ["Ca_metallicity"], None)),
        ("Metallicity_10", ("", ["Fe_metallicity"], None)),
    )

    def __init__(self, *args, **kwargs):
        super(SPHFieldInfo, self).__init__(*args, **kwargs)
        if ("PartType0", "Metallicity_00") in self.field_list:
            self.nuclei_names = metal_elements
            self.species_names = ["H", "H_p1"] + metal_elements

    def setup_particle_fields(self, ptype):
        FieldInfoContainer.setup_particle_fields(self, ptype)
        if ptype in ("PartType0",):
            self.setup_gas_particle_fields(ptype)
            setup_species_fields(self, ptype)

    def setup_gas_particle_fields(self, ptype):
        super(GizmoFieldInfo, self).setup_gas_particle_fields(ptype)

        def _h_density(field, data):
            x_H = 1.0 - data[(ptype, "He_metallicity")] - \
              data[(ptype, "metallicity")]
            return x_H * data[(ptype, "density")] * \
              data[(ptype, "NeutralHydrogenAbundance")]

        self.add_field(
            (ptype, "H_density"),
            function=_h_density,
            particle_type=True,
            units=self.ds.unit_system["density"])
        add_species_field_by_density(self, ptype, "H", particle_type=True)
        for suffix in ["density", "fraction", "mass", "number_density"]:
            self.alias((ptype, "H_p0_%s" % suffix), (ptype, "H_%s" % suffix))

        def _h_p1_density(field, data):
            x_H = 1.0 - data[(ptype, "He_metallicity")] - \
              data[(ptype, "metallicity")]
            return x_H * data[(ptype, "density")] * \
              (1.0 - data[(ptype, "NeutralHydrogenAbundance")])

        self.add_field(
            (ptype, "H_p1_density"),
            function=_h_p1_density,
            particle_type=True,
            units=self.ds.unit_system["density"])
        add_species_field_by_density(self, ptype, "H_p1", particle_type=True)

        def _nuclei_mass_density_field(field, data):
            species = field.name[1][:field.name[1].find("_")]
            return data[ptype, "density"] * \
              data[ptype, "%s_metallicity" % species]

        num_neighbors = 64
        for species in ['H', 'H_p0', 'H_p1']:
            for suf in ["_density", "_number_density"]:
                field = "%s%s" % (species, suf)
                fn = add_volume_weighted_smoothed_field(
                    ptype, "particle_position", "particle_mass",
                    "smoothing_length", "density", field,
                    self, num_neighbors)
                self.alias(("gas", field), fn[0])

        for species in self.nuclei_names:
            self.add_field(
                (ptype, "%s_nuclei_mass_density" % species),
                function=_nuclei_mass_density_field,
                particle_type=True,
                units=self.ds.unit_system["density"])

            for suf in ["_nuclei_mass_density", "_metallicity"]:
                field = "%s%s" % (species, suf)
                fn = add_volume_weighted_smoothed_field(
                    ptype, "particle_position", "particle_mass",
                    "smoothing_length", "density", field,
                    self, num_neighbors)

                self.alias(("gas", field), fn[0])
