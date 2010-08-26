"""
FLASH-specific IO functions

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: UCSD
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2010 Matthew Turk.  All Rights Reserved.

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

from yt.utilities.io_handler import \
    BaseIOHandler

class IOHandlerFLASH(BaseIOHandler):
    _data_style = "flash_hdf5"

    def __init__(self, *args, **kwargs):
        BaseIOHandler.__init__(self, *args, **kwargs)

    def _read_data_set(self, grid, field):
        f = h5py.File(grid.pf.parameter_filename, "r")
        tr = f["/%s" % field][grid.id - grid._id_offset,:,:,:].transpose()
        return tr.astype("float64")

    def _read_data_slice(self, grid, field, axis, coord):
        sl = [slice(None), slice(None), slice(None)]
        sl[axis] = slice(coord, coord + 1)
        f = h5py.File(grid.pf.parameter_filename, "r")
        tr = f["/%s" % field][grid.id - grid._id_offset].transpose()[sl]
        return tr.astype("float64")

