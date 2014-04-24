"""
Data structures for a generic SDF frontend




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import h5py
import numpy as np
import stat
import weakref
import struct
import glob
import time
import os
import types

from yt.utilities.logger import ytLogger as mylog
from yt.geometry.particle_geometry_handler import \
    ParticleIndex
from yt.data_objects.static_output import \
    Dataset, ParticleFile
from yt.utilities.physical_constants import \
    G, \
    cm_per_kpc, \
    mass_sun_cgs
from yt.utilities.cosmology import Cosmology
from .fields import \
    SDFFieldInfo

class SDFFile(ParticleFile):
    def __init__(self, pf, io, filename, file_id):
        pass

class SDFParticleDataset(Dataset):
    pass
