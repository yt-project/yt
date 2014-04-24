"""
SDF-specific fields




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import os
import numpy as np

from yt.funcs import *

from yt.fields.field_info_container import \
    FieldInfoContainer

from yt.config import ytcfg
from yt.utilities.physical_constants import mh
from yt.fields.species_fields import \
    add_species_field_by_fraction, \
    add_species_field_by_density, \
    setup_species_fields

from yt.fields.particle_fields import \
    add_volume_weighted_smoothed_field

class SDFFieldInfo(FieldInfoContainer):
    known_other_fields = ()

    known_particle_fields = (
        ("Mass", ("code_mass", ["particle_mass"], None)),
        ("Masses", ("code_mass", ["particle_mass"], None)),
        ("Coordinates", ("code_length", ["particle_position"], None)),
        ("Velocity", ("code_velocity", ["particle_velocity"], None)),
        ("Velocities", ("code_velocity", ["particle_velocity"], None)),
        ("ParticleIDs", ("", ["particle_index"], None)),
        ("InternalEnergy", ("", ["thermal_energy"], None)),
        ("SmoothingLength", ("code_length", ["smoothing_length"], None)),
        ("Density", ("code_mass / code_length**3", ["density"], None)),
        ("MaximumTemperature", ("K", [], None)),
        ("Temperature", ("K", ["temperature"], None)),
        ("Epsilon", ("code_length", [], None)),
        ("Metals", ("code_metallicity", ["metallicity"], None)),
        ("Phi", ("code_length", [], None)),
        ("FormationTime", ("code_time", ["creation_time"], None)),
    )
