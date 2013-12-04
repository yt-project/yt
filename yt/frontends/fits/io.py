"""
FITS-specific IO functions
"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
try:
    import astropy.io.fits as pyfits
except ImportError:
    pass

from yt.utilities.math_utils import prec_accum

from yt.utilities.io_handler import \
    BaseIOHandler
from yt.utilities.logger import ytLogger as mylog

class IOHandlerFITS(BaseIOHandler):
    _particle_reader = False
    _data_style = "fits"

    def __init__(self, pf):
        super(IOHandlerFITS, self).__init__(pf)
        self.pf = pf
        self._handle = pf._handle
        
    def _read_particles(self, fields_to_read, type, args, grid_list,
            count_list, conv_factors):
        pass

    def _read_data_set(self, grid, field):
        f = self._handle
        if self.pf.dimensionality == 2:
            nx,ny = f[field].data.tranpose().shape
            tr = f[field].data.transpose().reshape(nx,ny,1)
        elif self.pf.dimensionality == 3:
            tr = f[field].data.transpose()
        return tr.astype("float64")

    def _read_data_slice(self, grid, field, axis, coord):
        sl = [slice(None), slice(None), slice(None)]
        sl[axis] = slice(coord, coord + 1)
        f = self._handle
        if self.pf.dimensionality == 2:
            nx,ny = f[field].data.transpose().shape
            tr = f[field].data.transpose().reshape(nx,ny,1)[sl]
        elif self.pf.dimensionality == 3:
            tr = f[field].data.transpose()[sl]
        return tr.astype("float64")

    def _read_fluid_selection(self, chunks, selector, fields, size):
        chunks = list(chunks)
        if any((ftype != "gas" for ftype, fname in fields)):
            raise NotImplementedError
        f = self._handle
        rv = {}
        dt = "float64"
        for field in fields:
            rv[field] = np.empty(size, dtype=dt)
        ng = sum(len(c.objs) for c in chunks)
        mylog.debug("Reading %s cells of %s fields in %s blocks",
                    size, [f2 for f1, f2 in fields], ng)
        for field in fields:
            ftype, fname = field
            ds = f[fname].data.astype("float64")
            ind = 0
            for chunk in chunks:
                for g in chunk.objs:
                    if self.pf.dimensionality == 2:
                        nx,ny = ds.transpose().shape
                        data = ds.transpose().reshape(nx,ny,1)
                    elif self.pf.dimensionality == 3:
                        data = ds.transpose()
                    ind += g.select(selector, data, rv[field], ind) # caches
        return rv
