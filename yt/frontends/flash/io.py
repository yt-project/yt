"""
FLASH-specific IO functions



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
import h5py

from yt.utilities.io_handler import \
    BaseIOHandler

class IOHandlerFLASH(BaseIOHandler):
    _particle_reader = False
    _data_style = "flash_hdf5"

    def __init__(self, pf, *args, **kwargs):
        self._num_per_stride = kwargs.pop("num_per_stride", 1000000)
        BaseIOHandler.__init__(self, *args, **kwargs)
        # Now we cache the particle fields
        self.pf = pf
        self._handle = pf._handle
        self._particle_handle = pf._particle_handle
        
        try :
            particle_fields = [s[0].strip() for s in
                               self._particle_handle["/particle names"][:]]
            self._particle_fields = dict([("particle_" + s, i) for i, s in
                                          enumerate(particle_fields)])
        except KeyError:
            self._particle_fields = {}

    def _read_particles(self, fields_to_read, type, args, grid_list,
            count_list, conv_factors):
        pass

    def _read_data(self, grid, field):
        f = self._handle
        f_part = self._particle_handle
        if field in self._particle_fields:
            if grid.NumberOfParticles == 0: return np.array([], dtype='float64')
            start = self.pf.h._particle_indices[grid.id - grid._id_offset]
            end = self.pf.h._particle_indices[grid.id - grid._id_offset + 1]
            fi = self._particle_fields[field]
            tr = f_part["/tracer particles"][start:end, fi]
        else:
            tr = f["/%s" % field][grid.id - grid._id_offset,:,:,:].transpose()
        return tr.astype("float64")
