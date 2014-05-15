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

    def __init__(self, pf):
        super(IOHandlerFITS, self).__init__(pf)
        self.pf = pf
        self._handle = pf._handle
        if self.pf.line_width is not None:
            self.line_db = self.pf.line_database
            self.dz = self.pf.line_width/self.domain_dimensions[self.pf.spec_axis]
        else:
            self.line_db = None
            self.dz = 1.

    def _read_particles(self, fields_to_read, type, args, grid_list,
            count_list, conv_factors):
        pass

    def _read_particle_coords(self, chunks, ptf):
        pdata = self.pf._handle[self.pf.first_image].data
        assert(len(ptf) == 1)
        ptype = ptf.keys()[0]
        x = np.asarray(pdata.field("X"), dtype="=f8")
        y = np.asarray(pdata.field("Y"), dtype="=f8")
        z = np.ones(x.shape)
        x = (x-0.5)/self.pf.reblock+0.5
        y = (y-0.5)/self.pf.reblock+0.5
        yield ptype, (x,y,z)

    def _read_particle_fields(self, chunks, ptf, selector):
        pdata = self.pf._handle[self.pf.first_image].data
        assert(len(ptf) == 1)
        ptype = ptf.keys()[0]
        field_list = ptf[ptype]
        x = np.asarray(pdata.field("X"), dtype="=f8")
        y = np.asarray(pdata.field("Y"), dtype="=f8")
        z = np.ones(x.shape)
        x = (x-0.5)/self.pf.reblock+0.5
        y = (y-0.5)/self.pf.reblock+0.5
        mask = selector.select_points(x, y, z, 0.0)
        if mask is None: return
        for field in field_list:
            fd = field.split("_")[-1]
            data = pdata.field(fd.upper())
            if fd in ["x","y"]:
                data = (data.copy()-0.5)/self.pf.reblock+0.5
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
        dx = self.pf.domain_width/self.pf.domain_dimensions
        for field in fields:
            ftype, fname = field
            tmp_fname = fname
            if fname in self.pf.line_database:
                fname = self.pf.field_list[0][1]
            f = self.pf.index._file_map[fname]
            ds = f[self.pf.index._ext_map[fname]]
            bzero, bscale = self.pf.index._scale_map[fname]
            fname = tmp_fname
            ind = 0
            for chunk in chunks:
                for g in chunk.objs:
                    start = ((g.LeftEdge-self.pf.domain_left_edge)/dx).to_ndarray().astype("int")
                    end = start + g.ActiveDimensions
                    if self.line_db is not None and fname in self.line_db:
                        my_off = self.line_db.get(fname).in_units(self.pf.spec_unit).value
                        my_off = my_off - 0.5*self.pf.line_width
                        my_off = int((my_off-self.pf.freq_begin)/self.dz)
                        my_off = max(my_off, 0)
                        my_off = min(my_off, self.pf.dims[self.pf.spec_axis]-1)
                        start[self.pf.spec_axis] += my_off
                        end[self.pf.spec_axis] += my_off
                        mylog.debug("Reading from " + str(start) + str(end))
                    slices = [slice(start[i],end[i]) for i in xrange(3)]
                    if self.pf.reversed:
                        new_start = self.pf.dims[self.pf.spec_axis]-1-start[self.pf.spec_axis]
                        new_end = max(self.pf.dims[self.pf.spec_axis]-1-end[self.pf.spec_axis],0)
                        slices[self.pf.spec_axis] = slice(new_start,new_end,-1)
                    if self.pf.dimensionality == 2:
                        nx, ny = g.ActiveDimensions[:2]
                        nz = 1
                        data = np.zeros((nx,ny,nz))
                        data[:,:,0] = ds.data[slices[1],slices[0]].transpose()
                    elif self.pf.naxis == 4:
                        idx = self.pf.index._axis_map[fname]
                        data = ds.data[idx,slices[2],slices[1],slices[0]].transpose()
                    else:
                        data = ds.data[slices[2],slices[1],slices[0]].transpose()
                    if self.line_db is not None:
                        nz1 = data.shape[self.pf.spec_axis]
                        nz2 = g.ActiveDimensions[self.pf.spec_axis]
                        if nz1 != nz2:
                            old_data = data.copy()
                            data = np.zeros(g.ActiveDimensions)
                            data[:,:,nz2-nz1:] = old_data
                    if fname in self.pf.nan_mask:
                        data[np.isnan(data)] = self.pf.nan_mask[fname]
                    elif "all" in self.pf.nan_mask:
                        data[np.isnan(data)] = self.pf.nan_mask["all"]
                    data = bzero + bscale*data
                    ind += g.select(selector, data.astype("float64"), rv[field], ind)
        return rv
