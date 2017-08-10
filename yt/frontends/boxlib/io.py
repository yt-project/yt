"""
AMReX/Boxlib data-file handling functions



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import os
import numpy as np
from collections import defaultdict

from yt.utilities.io_handler import \
    BaseIOHandler
from yt.funcs import mylog
from yt.frontends.chombo.io import parse_orion_sinks
from yt.geometry.selection_routines import GridSelector

def _remove_raw(all_fields, raw_fields):
    centered_fields = set(all_fields)
    for raw in raw_fields:
        centered_fields.discard(raw)
    return list(centered_fields)

class IOHandlerBoxlib(BaseIOHandler):

    _dataset_type = "boxlib_native"

    def __init__(self, ds, *args, **kwargs):
        self.ds = ds

    def _read_fluid_selection(self, chunks, selector, fields, size):
        chunks = list(chunks)
        if any(( not (ftype == "boxlib" or ftype == 'raw') for ftype, fname in fields)):
            raise NotImplementedError
        rv = {}
        raw_fields = []
        for field in fields:
            if field[0] == "raw":
                nodal_flag = self.ds.nodal_flags[field[1]]
                num_nodes = 2**sum(nodal_flag)
                rv[field] = np.empty((size, num_nodes), dtype="float64")
                raw_fields.append(field)
            else:
                rv[field] = np.empty(size, dtype="float64")
        centered_fields = _remove_raw(fields, raw_fields)
        ng = sum(len(c.objs) for c in chunks)
        mylog.debug("Reading %s cells of %s fields in %s grids",
                    size, [f2 for f1, f2 in fields], ng)
        ind = 0
        for chunk in chunks:
            data = self._read_chunk_data(chunk, centered_fields)
            for g in chunk.objs:
                for field in fields:
                    if field in centered_fields:
                        ds = data[g.id].pop(field)
                    else:
                        ds = self._read_raw_field(g, field)
                    nd = g.select(selector, ds, rv[field], ind) # caches
                ind += nd
                data.pop(g.id)
        return rv

    def _read_raw_field(self, grid, field):
        field_name = field[1]
        base_dir = self.ds.index.raw_file

        box_list = self.ds.index.raw_field_map[field_name][0]
        fn_list = self.ds.index.raw_field_map[field_name][1]
        offset_list = self.ds.index.raw_field_map[field_name][2]

        lev = grid.Level        
        filename = base_dir + "Level_%d/" % lev + fn_list[grid.id]
        offset = offset_list[grid.id]
        box = box_list[grid.id]
        
        lo = box[0]
        hi = box[1]
        shape = hi - lo + 1
        with open(filename, "rb") as f:
            f.seek(offset)
            f.readline()  # always skip the first line
            arr = np.fromfile(f, 'float64', np.product(shape))
            arr = arr.reshape(shape, order='F')
        return arr

    def _read_chunk_data(self, chunk, fields):
        data = {}
        grids_by_file = defaultdict(list)
        if len(chunk.objs) == 0: return data
        for g in chunk.objs:
            if g.filename is None:
                continue
            grids_by_file[g.filename].append(g)
        dtype = self.ds.index._dtype
        bpr = dtype.itemsize
        for filename in grids_by_file:
            grids = grids_by_file[filename]
            grids.sort(key = lambda a: a._offset)
            f = open(filename, "rb")
            for grid in grids:
                data[grid.id] = {}
                local_offset = grid._get_offset(f) - f.tell()
                count = grid.ActiveDimensions.prod()
                size = count * bpr
                for field in self.ds.index.field_order:
                    if field in fields:
                        # We read it ...
                        f.seek(local_offset, os.SEEK_CUR)
                        v = np.fromfile(f, dtype=dtype, count=count)
                        v = v.reshape(grid.ActiveDimensions, order='F')
                        data[grid.id][field] = v
                        local_offset = 0
                    else:
                        local_offset += size
        return data

    def _read_particle_selection(self, chunks, selector, fields):
        rv = {}
        chunks = list(chunks)
        unions = self.ds.particle_unions

        rv = {f: np.array([]) for f in fields}
        for chunk in chunks:
            for grid in chunk.objs:
                for ftype, fname in fields:
                    if ftype in unions:
                        for subtype in unions[ftype]:
                            data = self._read_particles(grid, selector,
                                                        subtype, fname)
                            rv[ftype, fname] = np.concatenate((data,
                                                               rv[ftype, fname]))
                    else:
                        data = self._read_particles(grid, selector,
                                                    ftype, fname)
                        rv[ftype, fname] = np.concatenate((data,
                                                           rv[ftype, fname]))
        return rv

    def _read_particles(self, grid, selector, ftype, name):

        npart = grid._pdata[ftype]["NumberOfParticles"]
        if npart == 0:
            return np.array([])

        fn = grid._pdata[ftype]["particle_filename"]
        offset = grid._pdata[ftype]["offset"]
        pheader = self.ds.index.particle_headers[ftype]
        
        # handle the case that this is an integer field
        int_fnames = [fname for _, fname in pheader.known_int_fields]
        if name in int_fnames:
            ind = int_fnames.index(name)
            fn = grid._pdata[ftype]["particle_filename"]
            with open(fn, "rb") as f:

                # read in the position fields for selection
                f.seek(offset + 
                       pheader.particle_int_dtype.itemsize * npart)
                rdata = np.fromfile(f, pheader.real_type, pheader.num_real * npart)
                x = np.asarray(rdata[0::pheader.num_real], dtype=np.float64)
                y = np.asarray(rdata[1::pheader.num_real], dtype=np.float64)
                if (grid.ds.dimensionality == 2):
                    z = np.ones_like(y)
                    z *= 0.5*(grid.LeftEdge[2] + grid.RightEdge[2])
                else:
                    z = np.asarray(rdata[2::pheader.num_real], dtype=np.float64)
                mask = selector.select_points(x, y, z, 0.0)

                if mask is None:
                    return np.array([])
                
                # read in the data we want
                f.seek(offset)
                idata = np.fromfile(f, pheader.int_type, pheader.num_int * npart)
                data = np.asarray(idata[ind::pheader.num_int], dtype=np.float64)
                return data[mask].flatten()

        # handle case that this is a real field
        real_fnames = [fname for _, fname in pheader.known_real_fields]
        if name in real_fnames:
            ind = real_fnames.index(name)
            with open(fn, "rb") as f:

                # read in the position fields for selection
                f.seek(offset + 
                       pheader.particle_int_dtype.itemsize * npart)
                rdata = np.fromfile(f, pheader.real_type, pheader.num_real * npart)
                x = np.asarray(rdata[0::pheader.num_real], dtype=np.float64)
                y = np.asarray(rdata[1::pheader.num_real], dtype=np.float64)
                if (grid.ds.dimensionality == 2):
                    z = np.ones_like(y)
                    z *= 0.5*(grid.LeftEdge[2] + grid.RightEdge[2])
                else:
                    z = np.asarray(rdata[2::pheader.num_real], dtype=np.float64)
                mask = selector.select_points(x, y, z, 0.0)

                if mask is None:
                    return np.array([])

                data = np.asarray(rdata[ind::pheader.num_real], dtype=np.float64)
                return data[mask].flatten()


class IOHandlerOrion(IOHandlerBoxlib):
    _dataset_type = "orion_native"

    _particle_filename = None
    @property
    def particle_filename(self):
        fn = self.ds.output_dir + "/StarParticles"
        if not os.path.exists(fn):
            fn = self.ds.output_dir + "/SinkParticles"
        self._particle_filename = fn
        return self._particle_filename

    _particle_field_index = None
    @property
    def particle_field_index(self):

        index = parse_orion_sinks(self.particle_filename)

        self._particle_field_index = index
        return self._particle_field_index

    def _read_particle_selection(self, chunks, selector, fields):
        rv = {}
        chunks = list(chunks)

        if isinstance(selector, GridSelector):

            if not (len(chunks) == len(chunks[0].objs) == 1):
                raise RuntimeError

            grid = chunks[0].objs[0]

            for ftype, fname in fields:
                rv[ftype, fname] = self._read_particles(grid, fname)

            return rv

        rv = {f: np.array([]) for f in fields}
        for chunk in chunks:
            for grid in chunk.objs:
                for ftype, fname in fields:
                    data = self._read_particles(grid, fname)
                    rv[ftype, fname] = np.concatenate((data, rv[ftype, fname]))
        return rv

    def _read_particles(self, grid, field):
        """
        parses the Orion Star Particle text files

        """

        particles = []

        if grid.NumberOfParticles == 0:
            return np.array(particles)

        def read(line, field):
            entry = line.strip().split(' ')[self.particle_field_index[field]]
            return np.float(entry)

        try:
            lines = self._cached_lines
            for num in grid._particle_line_numbers:
                line = lines[num]
                particles.append(read(line, field))
            return np.array(particles)
        except AttributeError:
            fn = self.particle_filename
            with open(fn, 'r') as f:
                lines = f.readlines()
                self._cached_lines = lines
                for num in grid._particle_line_numbers:
                    line = lines[num]
                    particles.append(read(line, field))
            return np.array(particles)
