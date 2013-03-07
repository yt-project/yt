"""
Octree geometry handler

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2007-2011 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import h5py
import numpy as np
import string, re, gc, time, cPickle
import weakref

from itertools import chain, izip

from yt.funcs import *
from yt.utilities.logger import ytLogger as mylog
from yt.arraytypes import blankRecordArray
from yt.config import ytcfg
from yt.data_objects.field_info_container import NullFunc
from yt.geometry.geometry_handler import GeometryHandler, YTDataChunk
from yt.utilities.definitions import MAXLEVEL
from yt.utilities.io_handler import io_registry
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelAnalysisInterface, parallel_splitter

from yt.data_objects.data_containers import data_object_registry

class OctreeGeometryHandler(GeometryHandler):

    def _setup_geometry(self):
        mylog.debug("Initializing Octree Geometry Handler.")
        self._initialize_oct_handler()

    def get_smallest_dx(self):
        """
        Returns (in code units) the smallest cell size in the simulation.
        """
        return (self.parameter_file.domain_width /
                (2**self.max_level)).min()

    def convert(self, unit):
        return self.parameter_file.conversion_factors[unit]

    def find_max(self, field, finest_levels = 3):
        """
        Returns (value, center) of location of maximum for a given field.
        """
        if (field, finest_levels) in self._max_locations:
            return self._max_locations[(field, finest_levels)]
        mv, pos = self.find_max_cell_location(field, finest_levels)
        self._max_locations[(field, finest_levels)] = (mv, pos)
        return mv, pos

    def find_max_cell_location(self, field, finest_levels = 3):
        source = self.all_data()
        if finest_levels is not False:
            source.min_level = self.max_level - finest_levels
        mylog.debug("Searching for maximum value of %s", field)
        max_val, maxi, mx, my, mz = \
            source.quantities["MaxLocation"](field)
        mylog.info("Max Value is %0.5e at %0.16f %0.16f %0.16f", 
              max_val, mx, my, mz)
        self.pf.parameters["Max%sValue" % (field)] = max_val
        self.pf.parameters["Max%sPos" % (field)] = "%s" % ((mx,my,mz),)
        return max_val, np.array((mx,my,mz), dtype='float64')

