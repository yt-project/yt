"""
Gadget-specific data-file handling function

Author: Christopher E Moody <juxtaposicion@gmail.com>
Affiliation: UC Santa Cruz
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2010 Christopher E Moody, Matthew Turk.  All Rights Reserved.

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

class IOHandlerGadget(BaseIOHandler):
    _data_style = 'gadget_hdf5'
    def _read_data_set(self, grid, field):
        adr = grid.Address
        fh = h5py.File(grid.filename,mode='r')
        if 'particles' in fh[adr].keys():
            adr2 = adr+'/particles'
            return fh[adr2][field]
        return None
    def _read_field_names(self,grid): 
        adr = grid.Address
        fh = h5py.File(grid.filename,mode='r')
        rets = cPickle.loads(fh['/root'].attrs['fieldnames'])
        return rets

    def _read_data_slice(self,grid, field, axis, coord):
        adr = grid.Address
        fh = h5py.File(grid.filename,mode='r')
        if 'particles' in fh[adr].keys():
            adr2 = adr+'/particles'
            return fh[adr2][field][coord,axis]
        return None

