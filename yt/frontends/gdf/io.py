"""
The data-file handling functions

Author: Samuel W. Skillman <samskillman@gmail.com>
Affiliation: University of Colorado at Boulder
Author: Matthew Turk <matthewturk@gmail.com>
Author: J. S. Oishi <jsoishi@gmail.com>
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
import numpy as np
from yt.funcs import \
    mylog
from yt.utilities.io_handler import \
    BaseIOHandler


def field_dname(grid_id, field_name):
    return "/data/grid_%010i/%s" % (grid_id, field_name)


# TODO all particle bits were removed
class IOHandlerGDFHDF5(BaseIOHandler):
    _data_style = "grid_data_format"
    _offset_string = 'data:offsets=0'
    _data_string = 'data:datatype=0'

    def __init__(self, pf, *args, **kwargs):
        # TODO check if _num_per_stride is needed
        self._num_per_stride = kwargs.pop("num_per_stride", 1000000)
        BaseIOHandler.__init__(self, *args, **kwargs)
        self.pf = pf
        self._handle = pf._handle


    def _read_data_set(self, grid, field):
        data = self._handle[field_dname(grid.id, field)][:, :, :]
        # TODO transpose data if needed (grid.pf.field_ordering)
        return data.astype("float64")

    def _read_data_slice(self, grid, field, axis, coord):
        slc = [slice(None), slice(None), slice(None)]
        slc[axis] = slice(coord, coord + 1)
        # TODO transpose data if needed
        data = self._handle[field_dname(grid.id, field)][slc]
        return data.astype("float64")

    def _read_fluid_selection(self, chunks, selector, fields, size):
        chunks = list(chunks)
        # TODO ????
        #if any((ftype != "gas" for ftype, fname in fields)):
        #    raise NotImplementedError
        fhandle = self._handle
        rv = {}
        for field in fields:
            ftype, fname = field
            rv[field] = np.empty(
                size, dtype=fhandle[field_dname(0, fname)].dtype)
        ngrids = sum(len(chunk.objs) for chunk in chunks)
        mylog.debug("Reading %s cells of %s fields in %s blocks",
                    size, [fname for ftype, fname in fields], ngrids)
        for field in fields:
            ftype, fname = field
            ind = 0
            for chunk in chunks:
                for grid in chunk.objs:
                    mask = grid.select(selector)  # caches
                    if mask is None:
                        continue
                    # TODO transpose if needed
                    data = fhandle[field_dname(grid.id, fname)][mask]
                    rv[field][ind:ind + data.size] = data
                    ind += data.size
        return rv
