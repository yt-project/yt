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
        ("mass", ("code_mass", ["particle_mass"], None)),
        ("x", ("code_length", ["particle_position_x"], None)),
        ("y", ("code_length", ["particle_position_y"], None)),
        ("z", ("code_length", ["particle_position_z"], None)),
        ("vx", ("code_velocity", ["particle_velocity_x"], None)),
        ("vy", ("code_velocity", ["particle_velocity_y"], None)),
        ("vz", ("code_velocity", ["particle_velocity_z"], None)),
        ("ident", ("", ["particle_index"], None)),
    )
