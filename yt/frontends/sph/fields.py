"""
SPH fields




"""
from __future__ import absolute_import

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
    )

    def setup_particle_fields(self, ptype, *args, **kwargs):
        super(SPHFieldInfo, self).setup_particle_fields(ptype, *args, **kwargs)
        setup_species_fields(self, ptype)
