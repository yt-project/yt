"""
SDF data-file handling function




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import glob
import h5py
import numpy as np
from yt.funcs import *
from yt.utilities.exceptions import *

from yt.utilities.io_handler import \
    BaseIOHandler

from yt.utilities.fortran_utils import read_record
from yt.utilities.lib.geometry_utils import compute_morton

from yt.geometry.oct_container import _ORDER_MAX

class IOHandlerSDF(BaseIOHandler):
    _dataset_type = "SDF"

    def _read_fluid_selection(self, chunks, selector, fields, size):
        raise NotImplementedError

    def _read_particle_coords(self, chunks, ptf):
        pass

    def _read_particle_fields(self, chunks, ptf, selector):
        pass

    def _initialize_index(self, data_file, regions):
        pass

    def _count_particles(self, data_file):
        pass

    def _identify_fields(self, data_file):
        return fields, {}

