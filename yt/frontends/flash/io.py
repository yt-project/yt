"""
FLASH-specific IO functions

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: UCSD
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

import numpy as na
import h5py

from yt.utilities.io_handler import \
    BaseIOHandler
from yt.utilities.logger import ytLogger as mylog

class IOHandlerFLASH(BaseIOHandler):
    _particle_reader = False
    _data_style = "flash_hdf5"

    def __init__(self, pf, *args, **kwargs):
        self._num_per_stride = kwargs.pop("num_per_stride", 1000000)
        BaseIOHandler.__init__(self, *args, **kwargs)
        # Now we cache the particle fields
        self.pf = pf
        self._handle = pf._handle
        try :
            particle_fields = [s[0].strip() for s in
                               self._handle["/particle names"][:]]
            self._particle_fields = dict([("particle_" + s, i) for i, s in
                                          enumerate(particle_fields)])
        except KeyError:
            self._particle_fields = {}

    def _read_particles(self, fields_to_read, type, args, grid_list,
            count_list, conv_factors):
        pass

    def _select_particles(self, grid, field):
        f = self._handle
        npart = f["/tracer particles"].shape[0]
        total_selected = 0
        start = 0
        stride = 1e6
        blki = self._particle_fields["particle_blk"]
        bi = grid.id - grid._id_offset
        fi = self._particle_fields[field]
        tr = []
        while start < npart:
            end = min(start + stride - 1, npart)
            gi = f["/tracer particles"][start:end,blki] == bi
            tr.append(f["/tracer particles"][gi,fi])
            start = end
        return na.concatenate(tr)

    def _read_data_set(self, grid, field):
        f = self._handle
        if field in self._particle_fields:
            if grid.NumberOfParticles == 0: return na.array([], dtype='float64')
            start = self.pf.h._particle_indices[grid.id - grid._id_offset]
            end = self.pf.h._particle_indices[grid.id - grid._id_offset + 1]
            fi = self._particle_fields[field]
            tr = f["/tracer particles"][start:end, fi]
        else:
            tr = f["/%s" % field][grid.id - grid._id_offset,:,:,:].transpose()
        return tr.astype("float64")

    def _read_data_slice(self, grid, field, axis, coord):
        sl = [slice(None), slice(None), slice(None)]
        sl[axis] = slice(coord, coord + 1)
        f = self._handle
        tr = f["/%s" % field][grid.id - grid._id_offset].transpose()[sl]
        return tr.astype("float64")

    def _read_fluid_selection(self, chunks, selector, fields, size):
        if any((ftype != "gas" for ftype, fname in fields)):
            raise NotImplementedError
        f = self._handle
        rv = {}
        for field in fields:
            ftype, fname = field
            rv[field] = na.empty(size, dtype=f["/%s" % fname].dtype)
        ng = sum(len(c.objs) for c in chunks)
        mylog.debug("Reading %s cells of %s fields in %s blocks",
                    size, [f2 for f1, f2 in fields], ng)
        for field in fields:
            ftype, fname = field
            ds = f["/%s" % fname]
            ind = 0
            for chunk in chunks:
                for g in chunk.objs:
                    mask = g.select(selector) # caches
                    if mask is None: continue
                    data = ds[g.id - g._id_offset,:,:,:].transpose()[mask]
                    rv[field][ind:ind+data.size] = data
                    ind += data.size
        return rv
