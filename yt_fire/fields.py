"""
FIRE-specific fields




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, Britton Smith.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.fields.species_fields import \
    add_species_field_by_density
from yt.frontends.gadget.fields import \
    GadgetFieldInfo
from yt.frontends.sph.fields import \
    SPHFieldInfo

class FIREFieldInfo(GadgetFieldInfo):
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
        ("Metallicity_01", ("", ["He_nuclei_fraction"], None)),
        ("Metallicity_02", ("", ["C_nuclei_fraction"], None)),
        ("Metallicity_03", ("", ["N_nuclei_fraction"], None)),
        ("Metallicity_04", ("", ["O_nuclei_fraction"], None)),
        ("Metallicity_05", ("", ["Ne_nuclei_fraction"], None)),
        ("Metallicity_06", ("", ["Mg_nuclei_fraction"], None)),
        ("Metallicity_07", ("", ["Si_nuclei_fraction"], None)),
        ("Metallicity_08", ("", ["S_nuclei_fraction"], None)),
        ("Metallicity_09", ("", ["Ca_nuclei_fraction"], None)),
        ("Metallicity_10", ("", ["Fe_nuclei_fraction"], None)),
    )

    def __init__(self, *args, **kwargs):
        super(SPHFieldInfo, self).__init__(*args, **kwargs)
        if ("PartType0", "Metallicity_00") in self.field_list:
            self.nuclei_names = ["He", "C", "N", "O", "Ne", "Mg", "Si", "S",
                                 "Ca", "Fe"]

    def setup_gas_particle_fields(self, ptype):
        super(FIREFieldInfo, self).setup_gas_particle_fields(ptype)
        self.alias((ptype, "temperature"), (ptype, "Temperature"))

        def _h_density(field, data):
            x_H = 1.0 - data[(ptype, "He_nuclei_fraction")] - \
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
            self.alias((ptype, "H_%s" % suffix), (ptype, "H_p0_%s" % suffix))

        def _h_p1_density(field, data):
            x_H = 1.0 - data[(ptype, "He_nuclei_fraction")] - \
              data[(ptype, "metallicity")]
            return x_H * data[(ptype, "density")] * \
              (1.0 - data[(ptype, "NeutralHydrogenAbundance")])

        self.add_field(
            (ptype, "H_p1_density"),
            function=_h_p1_density,
            particle_type=True,
            units=self.ds.unit_system["density"])
        add_species_field_by_density(self, ptype, "H_p1", particle_type=True)

        for species in self.nuclei_names:
            def _nuclei_density_field(field, data):
                return data[ptype, "density"] * \
                  data[ptype, "%s_nuclei_fraction" % species]

            self.add_field(
                (ptype, "%s_nuclei_density" % species),
                function=_h_p1_density,
                particle_type=True,
                units=self.ds.unit_system["density"])
