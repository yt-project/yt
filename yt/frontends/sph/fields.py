"""
SPH fields




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.fields.field_info_container import \
    FieldInfoContainer
from yt.fields.species_fields import \
    setup_species_fields

class SPHFieldInfo(FieldInfoContainer):
    known_other_fields = ()

    known_particle_fields = (
        ("Mass", ("code_mass", ["particle_mass"], None)),
        ("Masses", ("code_mass", ["particle_mass"], None)),
        ("Coordinates", ("code_length", ["particle_position"], None)),
        ("Velocity", ("code_velocity", ["particle_velocity"], None)),
        ("Velocities", ("code_velocity", ["particle_velocity"], None)),
        ("ParticleIDs", ("", ["particle_index"], None)),
        ("InternalEnergy", ("code_velocity**2", ["thermal_energy"], None)),
        ("SmoothingLength", ("code_length", ["smoothing_length"], None)),
        ("Density", ("code_mass / code_length**3", ["density"], None)),
        ("MaximumTemperature", ("K", [], None)),
        ("Temperature", ("K", ["temperature"], None)),
        ("Epsilon", ("code_length", [], None)),
        ("Metals", ("code_metallicity", ["metallicity"], None)),
        ("Metallicity", ("code_metallicity", ["metallicity"], None)),
        ("Phi", ("code_length", [], None)),
        ("StarFormationRate", ("code_mass / code_time", [], None)),
        ("FormationTime", ("code_time", ["creation_time"], None)),
        # These are metallicity fields that get discovered for FIRE simulations
        ("Metallicity_00", ("", ["metallicity"], None)),
        ("Metallicity_01", ("", ["He_fraction"], None)),
        ("Metallicity_02", ("", ["C_fraction"], None)),
        ("Metallicity_03", ("", ["N_fraction"], None)),
        ("Metallicity_04", ("", ["O_fraction"], None)),
        ("Metallicity_05", ("", ["Ne_fraction"], None)),
        ("Metallicity_06", ("", ["Mg_fraction"], None)),
        ("Metallicity_07", ("", ["Si_fraction"], None)),
        ("Metallicity_08", ("", ["S_fraction"], None)),
        ("Metallicity_09", ("", ["Ca_fraction"], None)),
        ("Metallicity_10", ("", ["Fe_fraction"], None)),
    )

    def __init__(self, *args, **kwargs):
        super(SPHFieldInfo, self).__init__(*args, **kwargs)
        # Special case for FIRE
        if ("PartType0", "Metallicity_00") in self.field_list:
            self.species_names += ["He", "C", "N", "O", "Ne", "Mg", "Si", "S",
                "Ca", "Fe"]

    def setup_particle_fields(self, ptype, *args, **kwargs):
        super(SPHFieldInfo, self).setup_particle_fields(ptype, *args, **kwargs)
        setup_species_fields(self, ptype)
