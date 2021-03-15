from itertools import groupby

import numpy as np

from yt.utilities.io_handler import BaseIOHandler
from yt.utilities.logger import ytLogger as mylog


# http://stackoverflow.com/questions/2361945/detecting-consecutive-integers-in-a-list
def grid_sequences(grids):
    g_iter = sorted(grids, key=lambda g: g.id)
    for _, g in groupby(enumerate(g_iter), lambda i_x1: i_x1[0] - i_x1[1].id):
        seq = list(v[1] for v in g)
        yield seq


ii = [0, 1, 0, 1, 0, 1, 0, 1]
jj = [0, 0, 1, 1, 0, 0, 1, 1]
kk = [0, 0, 0, 0, 1, 1, 1, 1]


class IOHandlerAthenaPP(BaseIOHandler):
    _particle_reader = False
    _dataset_type = "athena_pp"

    def __init__(self, ds):
        super().__init__(ds)
        self._handle = ds._handle

    def _read_particles(
        self, fields_to_read, type, args, grid_list, count_list, conv_factors
    ):
        pass

    def _read_fluid_selection(self, chunks, selector, fields, size):
        chunks = list(chunks)
        if any((ftype != "athena_pp" for ftype, fname in fields)):
            raise NotImplementedError
        f = self._handle
        rv = {}
        for field in fields:
            # Always use *native* 64-bit float.
            rv[field] = np.empty(size, dtype="=f8")
        ng = sum(len(c.objs) for c in chunks)
        mylog.debug(
            "Reading %s cells of %s fields in %s blocks",
            size,
            [f2 for f1, f2 in fields],
            ng,
        )
        last_dname = None
        for field in fields:
            ftype, fname = field
            dname, fdi = self.ds._field_map[fname]
            if dname != last_dname:
                ds = f[f"/{dname}"]
            ind = 0
            for chunk in chunks:
                if self.ds.logarithmic:
                    for mesh in chunk.objs:
                        nx, ny, nz = mesh.mesh_dims // self.ds.index.mesh_factors
                        data = np.empty(mesh.mesh_dims, dtype="=f8")
                        for n, id in enumerate(mesh.mesh_blocks):
                            data[
                                ii[n] * nx : (ii[n] + 1) * nx,
                                jj[n] * ny : (jj[n] + 1) * ny,
                                kk[n] * nz : (kk[n] + 1) * nz,
                            ] = ds[fdi, id, :, :, :].transpose()
                        ind += mesh.select(selector, data, rv[field], ind)  # caches
                else:
                    for gs in grid_sequences(chunk.objs):
                        start = gs[0].id - gs[0]._id_offset
                        end = gs[-1].id - gs[-1]._id_offset + 1
                        data = ds[fdi, start:end, :, :, :].transpose()
                        for i, g in enumerate(gs):
                            ind += g.select(selector, data[..., i], rv[field], ind)
            last_dname = dname
        return rv

    def _read_chunk_data(self, chunk, fields):
        if self.ds.logarithmic:
            pass
        f = self._handle
        rv = {}
        for g in chunk.objs:
            rv[g.id] = {}
        if len(fields) == 0:
            return rv
        for field in fields:
            ftype, fname = field
            dname, fdi = self.ds._field_map[fname]
            ds = f[f"/{dname}"]
            for gs in grid_sequences(chunk.objs):
                start = gs[0].id - gs[0]._id_offset
                end = gs[-1].id - gs[-1]._id_offset + 1
                data = ds[fdi, start:end, :, :, :].transpose()
                for i, g in enumerate(gs):
                    rv[g.id][field] = np.asarray(data[..., i], "=f8")
        return rv
