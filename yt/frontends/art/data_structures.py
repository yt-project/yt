"""
ART-specific data structures

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: UCSD
Author: Christopher Moody <cemoody@ucsc.edu>
Affiliation: UCSC
Homepage: http://yt-project.org/
License:
  Copyright (C) 2010-2011 Matthew Turk.  All Rights Reserved.

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
import numpy as np
import stat
import weakref
import cStringIO

from yt.funcs import *
from yt.data_objects.grid_patch import \
      AMRGridPatch
from yt.geometry.oct_geometry_handler import \
    OctreeGeometryHandler
from yt.geometry.geometry_handler import \
    GeometryHandler, YTDataChunk
from yt.data_objects.static_output import \
    StaticOutput

from .definitions import *
from yt.utilities.definitions import \
    mpc_conversion, sec_conversion
from yt.utilities.lib import \
    get_box_grids_level
from yt.utilities.io_handler import \
    io_registry
from yt.data_objects.field_info_container import \
    FieldInfoContainer, NullFunc
import yt.utilities.fortran_utils as fpu

class ARTDomainFile(object):
    _last_mask = None
    _last_selector_id = None
    _hydro_offset = None
    _level_count = None
    nvar = 6
    def __init__(self, pf, domain_id):
        self.pf = pf
        self.domain_id = domain_id
        self._read_amr_header()
    
	@property
    def level_count(self):
        if self._level_count is not None: return self._level_count
        self.hydro_offset
        return self._level_count

    @property
    def hydro_offset(self):
        if self._hydro_offset is not None: return self._hydro_offset
		#open the file and figure the level_count and offsets


