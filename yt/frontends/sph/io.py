"""
Gadget-specific data-file handling function

Author: Christopher E Moody <juxtaposicion@gmail.com>
Affiliation: UC Santa Cruz
Homepage: http://yt-project.org/
License:
  Copyright (C) 2010-2011 Christopher E Moody, Matthew Turk.  All Rights Reserved.

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
from yt.funcs import *

from yt.utilities.io_handler import \
    BaseIOHandler

class IOHandlerOWLS(BaseIOHandler):
    _data_style = "OWLS"

    def _read_fluid_selection(self, chunks, selector, fields, size):
        raise NotImplementedError

    def _read_particle_selection_by_type(self, chunks, selector, fields):
        rv = {}
        # We first need a set of masks for each particle type
        ptf = defaultdict(lambda: list)
        psize = defaultdict(lambda: 0)
        for ftype, fname in fields:
            ptf[ftype].append(fname)
        for chunk in chunks: # Will be OWLS domains
            for subset in chunk.objs:
                for ptype, field_list in sorted(ptf.items()):
                    f = h5py.File(subset.domain.filename, "r")
                    coords = f["/PartType%s/Coordinates" % ptype][:]
                    psize[ptype] += subset.count_particles(selector,
                                coords[:,0], coords[:,1], coords[:,2])
                    del coords
        # Now we have all the sizes, and we can allocate
        ind = {}
        for field in fields:
            rv[field] = np.empty(size, dtype="float64")
            ind[field] = 0
        for chunk in chunks: # Will be OWLS domains
            for subset in chunk.objs:
                for ptype, field_list in sorted(ptf.items()):
                    f = h5py.File(subset.domain.filename, "r")
                    g = f["/PartType%s" % ptype]
                    coords = g["Coordinates"][:]
                    mask = subset.select_particles(selector,
                                coords[:,0], coords[:,1], coords[:,2])
                    del coords
                    if mask is None: continue
                    for field in field_list:
                        data = g[field][mask,...]
                        my_ind = ind[ptype, field]
                        rv[ptype, field][my_ind:my_ind + data.size] = data
                        ind[ptype, field] += data.shape[0]
        return rv
