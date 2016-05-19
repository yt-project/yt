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

from yt.utilities.io_handler import \
    BaseIOHandler
from yt.utilities.logger import ytLogger as mylog

class IOHandlerFITS(BaseIOHandler):
    _particle_reader = False
    _dataset_type = "fits"

    def __init__(self, ds):
        super(IOHandlerFITS, self).__init__(ds)
        self.ds = ds
        self._handle = ds._handle

    def _read_particles(self, fields_to_read, type, args, grid_list,
            count_list, conv_factors):
        pass

    def _read_particle_coords(self, chunks, ptf):
        pdata = self.ds._handle[self.ds.first_image].data
        assert(len(ptf) == 1)
        ptype = list(ptf.keys())[0]
        x = np.asarray(pdata.field("X"), dtype="=f8")
        y = np.asarray(pdata.field("Y"), dtype="=f8")
        z = np.ones(x.shape)
        x = (x-0.5)/self.ds.reblock+0.5
        y = (y-0.5)/self.ds.reblock+0.5
        yield ptype, (x,y,z)

    def _read_particle_fields(self, chunks, ptf, selector):
        pdata = self.ds._handle[self.ds.first_image].data
        assert(len(ptf) == 1)
        ptype = list(ptf.keys())[0]
        field_list = ptf[ptype]
        x = np.asarray(pdata.field("X"), dtype="=f8")
        y = np.asarray(pdata.field("Y"), dtype="=f8")
        z = np.ones(x.shape)
        x = (x-0.5)/self.ds.reblock+0.5
        y = (y-0.5)/self.ds.reblock+0.5
        mask = selector.select_points(x, y, z, 0.0)
        if mask is None: return
        for field in field_list:
            fd = field.split("_")[-1]
            data = pdata.field(fd.upper())
            if fd in ["x","y"]:
                data = (data.copy()-0.5)/self.ds.reblock+0.5
            yield (ptype, field), data[mask]

    def _read_fluid_selection(self, chunks, selector, fields, size):
        chunks = list(chunks)
        if any((ftype != "fits" for ftype, fname in fields)):
            raise NotImplementedError
        rv = {}
        dt = "float64"
        for field in fields:
            rv[field] = np.empty(size, dtype=dt)
        ng = sum(len(c.objs) for c in chunks)
        mylog.debug("Reading %s cells of %s fields in %s grids",
                    size, [f2 for f1, f2 in fields], ng)
        dx = self.ds.domain_width/self.ds.domain_dimensions
        for field in fields:
            ftype, fname = field
            f = self.ds.index._file_map[fname]
            ds = f[self.ds.index._ext_map[fname]]
            bzero, bscale = self.ds.index._scale_map[fname]
            ind = 0
            for chunk in chunks:
                for g in chunk.objs:
                    start = ((g.LeftEdge-self.ds.domain_left_edge)/dx).to_ndarray().astype("int")
                    end = start + g.ActiveDimensions
                    slices = [slice(start[i],end[i]) for i in range(3)]
                    if self.ds.dimensionality == 2:
                        nx, ny = g.ActiveDimensions[:2]
                        nz = 1
                        data = np.zeros((nx,ny,nz))
                        data[:,:,0] = ds.data[slices[1],slices[0]].transpose()
                    elif self.ds.naxis == 4:
                        idx = self.ds.index._axis_map[fname]
                        data = ds.data[idx,slices[2],slices[1],slices[0]].transpose()
                    else:
                        data = ds.data[slices[2],slices[1],slices[0]].transpose()
                    if fname in self.ds.nan_mask:
                        data[np.isnan(data)] = self.ds.nan_mask[fname]
                    elif "all" in self.ds.nan_mask:
                        data[np.isnan(data)] = self.ds.nan_mask["all"]
                    data = bzero + bscale*data
                    ind += g.select(selector, data.astype("float64"), rv[field], ind)
        return rv
