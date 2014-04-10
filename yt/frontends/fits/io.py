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

from yt.utilities.math_utils import prec_accum

from yt.utilities.io_handler import \
    BaseIOHandler
from yt.utilities.logger import ytLogger as mylog

class IOHandlerFITS(BaseIOHandler):
    _particle_reader = False
    _dataset_type = "fits"

    def __init__(self, pf):
        super(IOHandlerFITS, self).__init__(pf)
        self.pf = pf
        self._handle = pf._handle
        self.folded = False
        if self.pf.folded_axis is not None:
            self.folded = True

    def _read_particles(self, fields_to_read, type, args, grid_list,
            count_list, conv_factors):
        pass

    def _read_fluid_selection(self, chunks, selector, fields, size):
        chunks = list(chunks)
        if any((ftype != "fits" for ftype, fname in fields)):
            raise NotImplementedError
        f = self._handle
        rv = {}
        dt = "float64"
        for field in fields:
            rv[field] = np.empty(size, dtype=dt)
        ng = sum(len(c.objs) for c in chunks)
        mylog.debug("Reading %s cells of %s fields in %s grids",
                    size, [f2 for f1, f2 in fields], ng)
        for field in fields:
            ftype, fname = field
            ds = f[fname]
            ind = 0
            for chunk in chunks:
                for g in chunk.objs:
                    centering = np.array([0.5]*3)
                    if self.folded:
                        centering[-1] = 0.0
                    start = (g.LeftEdge.ndarray_view()-centering).astype("int")
                    end = (g.RightEdge.ndarray_view()-centering).astype("int")
                    if self.folded:
                        my_off = self.pf.line_database.get(fname, 0)\
                            + self.pf.folded_width/2

                        start[-1] += my_off
                        end[-1] += my_off
                        mylog.debug("Reading from " + str(start) + str(end))
                    if self.pf.dimensionality == 2:
                        nx, ny = g.ActiveDimensions[:2]
                        nz = 1
                        data = np.zeros((nx,ny,nz))
                        data[:,:,0] = ds.data[start[1]:end[1],start[0]:end[0]].transpose()
                    elif self.pf.dimensionality == 3:
                        data = ds.data[start[2]:end[2],start[1]:end[1],start[0]:end[0]].transpose()
                    if self.pf.mask_nans: data[np.isnan(data)] = 0.0
                    ind += g.select(selector, data.astype("float64"), rv[field], ind)
        return rv

class IOHandlerFITSXYV(IOHandlerFITS):
    _particle_reader = False
    _dataset_type = "xyv_fits"

    def __init__(self, pf):
        super(IOHandlerFITSXYV,self).__init__(pf)

    def _read_fluid_selection(self, chunks, selector, fields, size):
        chunks = list(chunks)
        if any((ftype != "xyv_fits" for ftype, fname in fields)):
            raise NotImplementedError
        f = self._handle
        rv = {}
        dt = "float64"
        for field in fields:
            rv[field] = np.empty(size, dtype=dt)
        ng = sum(len(c.objs) for c in chunks)
        mylog.debug("Reading %s cells of %s fields in %s grids",
                    size, [f2 for f1, f2 in fields], ng)
        for field in fields:
            ftype, fname = field
            if self.pf.four_dims:
                #ds = f[fname.split("_")[0]]
                ds = f['primary']
            else:
                ds = f['primary']
            ind = 0
            for chunk in chunks:
                for g in chunk.objs:
                    centering = np.array([0.5]*3)
                    if self.folded:
                        centering[-1] = 0.0
                    start = (g.LeftEdge.ndarray_view()-centering).astype("int")
                    end = (g.RightEdge.ndarray_view()-centering).astype("int")
                    if self.folded:
                        my_off = self.pf.line_database.get(fname, 0)\
                            + self.pf.folded_width/2
                        start[-1] += my_off
                        end[-1] += my_off
                        mylog.debug("Reading from " + str(start) + str(end))
                    if self.pf.four_dims:
                        idx = self.pf.index._field_map[fname]
                        data = ds.data[idx,start[2]:end[2],start[1]:end[1],start[0]:end[0]].transpose()
                    else:
                        data = ds.data[start[2]:end[2],start[1]:end[1],start[0]:end[0]].transpose()
                    if self.pf.mask_nans: data[np.isnan(data)] = 0.0
                    ind += g.select(selector, data.astype("float64"), rv[field], ind)
        return rv
