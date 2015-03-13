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
        if self.ds.line_width is not None:
            self.line_db = self.ds.line_database
            self.dz = self.ds.line_width/self.domain_dimensions[self.ds.spec_axis]
        else:
            self.line_db = None
            self.dz = 1.

    def _read_particles(self, fields_to_read, type, args, grid_list,
            count_list, conv_factors):
        pass

    def _read_particle_coords(self, chunks, ptf):
        pdata = self.ds._handle[self.ds.first_image].data
        assert(len(ptf) == 1)
        ptype = ptf.keys()[0]
        x = np.asarray(pdata.field("X"), dtype="=f8")
        y = np.asarray(pdata.field("Y"), dtype="=f8")
        z = np.ones(x.shape)
        x = (x-0.5)/self.ds.reblock+0.5
        y = (y-0.5)/self.ds.reblock+0.5
        yield ptype, (x,y,z)

    def _read_particle_fields(self, chunks, ptf, selector):
        pdata = self.ds._handle[self.ds.first_image].data
        assert(len(ptf) == 1)
        ptype = ptf.keys()[0]
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
            tmp_fname = fname
            if fname in self.ds.line_database:
                fname = self.ds.field_list[0][1]
            f = self.ds.index._file_map[fname]
            ds = f[self.ds.index._ext_map[fname]]
            bzero, bscale = self.ds.index._scale_map[fname]
            fname = tmp_fname
            ind = 0
            for chunk in chunks:
                for g in chunk.objs:
                    start = ((g.LeftEdge-self.ds.domain_left_edge)/dx).to_ndarray().astype("int")
                    end = start + g.ActiveDimensions
                    if self.line_db is not None and fname in self.line_db:
                        my_off = self.line_db.get(fname).in_units(self.ds.spec_unit).value
                        my_off = my_off - 0.5*self.ds.line_width
                        my_off = int((my_off-self.ds.freq_begin)/self.dz)
                        my_off = max(my_off, 0)
                        my_off = min(my_off, self.ds.dims[self.ds.spec_axis]-1)
                        start[self.ds.spec_axis] += my_off
                        end[self.ds.spec_axis] += my_off
                        mylog.debug("Reading from " + str(start) + str(end))
                    slices = [slice(start[i],end[i]) for i in range(3)]
                    if self.ds.reversed:
                        new_start = self.ds.dims[self.ds.spec_axis]-1-start[self.ds.spec_axis]
                        new_end = max(self.ds.dims[self.ds.spec_axis]-1-end[self.ds.spec_axis],0)
                        slices[self.ds.spec_axis] = slice(new_start,new_end,-1)
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
                    if self.line_db is not None:
                        nz1 = data.shape[self.ds.spec_axis]
                        nz2 = g.ActiveDimensions[self.ds.spec_axis]
                        if nz1 != nz2:
                            old_data = data.copy()
                            data = np.zeros(g.ActiveDimensions)
                            data[:,:,nz2-nz1:] = old_data
                    if fname in self.ds.nan_mask:
                        data[np.isnan(data)] = self.ds.nan_mask[fname]
                    elif "all" in self.ds.nan_mask:
                        data[np.isnan(data)] = self.ds.nan_mask["all"]
                    data = bzero + bscale*data
                    ind += g.select(selector, data.astype("float64"), rv[field], ind)
        return rv
