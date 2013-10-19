"""
Gadget-specific data-file handling function



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

from yt.utilities.io_handler import \
    BaseIOHandler

class IOHandlerGadget(BaseIOHandler):
    _data_style = 'gadget_infrastructure'
    def _read_data(self, grid, field):
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
        dat = self._read_data(grid,field)
        return dat[coord]

