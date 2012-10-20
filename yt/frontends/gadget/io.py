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

from yt.utilities.io_handler import \
    BaseIOHandler

class IOHandlerGadget(BaseIOHandler):
    _data_style = 'gadget_infrastructure'
    def _read_data_set(self, grid, field):
        data = []
        fh = h5py.File(grid.filename,mode='r')
        for ptype in grid.particle_types:
            address = '/data/grid_%010i/particles/%s/%s' % (grid.id, ptype, field)
            data.append(fh[address][:])
        if len(data) > 0:
            data = np.concatenate(data)
        fh.close()
        return np.array(data)
    def _read_field_names(self,grid): 
        adr = grid.Address
        fh = h5py.File(grid.filename,mode='r')
        rets = cPickle.loads(fh['/root'].attrs['fieldnames'])
        fh.close()
        return rets

    def _read_data_slice(self,grid, field, axis, coord):
        #how would we implement axis here?
        dat = self._read_data_set(grid,field)
        return dat[coord]

