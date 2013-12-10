"""
OWLS-specific fields




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.funcs import *
from yt.fields.field_info_container import \
    FieldInfoContainer

class GadgetHDF5FieldInfo(FieldInfoContainer):
    known_other_fields = ()

    known_particle_fields = (
        ("Coordinates", ("code_length", [], None)),
        ("Velocities", ("code_length / code_time", [], None)),
        ("Velocity", ("code_length / code_time", [], None)),
        ("Mass", ("code_mass", [], None)),
        ("ParticleIDs", ("", ["particle_ids", "particle_identifiers"], None)),
    )

class GadgetFieldInfo(FieldInfoContainer):
    known_other_fields = ()

    known_particle_fields = ()

class TipsyFieldInfo(FieldInfoContainer):
    known_other_fields = ()

    known_particle_fields = ()

class OWLSFieldInfo(FieldInfoContainer):
    known_other_fields = ()
    known_particle_fields = ()
