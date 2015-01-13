"""
Octree geometry handler




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
import string, re, gc, time, cPickle
import weakref

from itertools import chain, izip

from yt.funcs import *
from yt.utilities.logger import ytLogger as mylog
from yt.arraytypes import blankRecordArray
from yt.config import ytcfg
from yt.fields.field_info_container import NullFunc
from yt.geometry.geometry_handler import Index, YTDataChunk
from yt.utilities.definitions import MAXLEVEL
from yt.utilities.io_handler import io_registry
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelAnalysisInterface

from yt.data_objects.data_containers import data_object_registry

class OctreeIndex(Index):
    """The Index subclass for oct AMR datasets"""
    def _setup_geometry(self):
        mylog.debug("Initializing Octree Geometry Handler.")
        self._initialize_oct_handler()

    def get_smallest_dx(self):
        """
        Returns (in code units) the smallest cell size in the simulation.
        """
        return (self.dataset.domain_width /
                (self.dataset.domain_dimensions * 2**(self.max_level))).min()

    def convert(self, unit):
        return self.dataset.conversion_factors[unit]
