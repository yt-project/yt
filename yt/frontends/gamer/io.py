from itertools import groupby

import numpy as np

from yt.geometry.selection_routines import AlwaysSelector
from yt.utilities.io_handler import BaseIOHandler
from yt.utilities.logger import ytLogger as mylog

# -----------------------------------------------------------------------------
# GAMER shares a similar HDF5 format, and thus io.py as well, with FLASH
# -----------------------------------------------------------------------------


# group grids with consecutive indices together to improve the I/O performance
# --> grids are assumed to be sorted into ascending numerical order already
def grid_sequences(grids):
    for _k, g in groupby(enumerate(grids), lambda i_x: i_x[0] - i_x[1].id):
        seq = list(v[1] for v in g)
        yield seq


def particle_sequences(grids):
    for _k, g in groupby(enumerate(grids), lambda i_x: i_x[0] - i_x[1].id):
        seq = list(v[1] for v in g)
        yield seq[0], seq[-1]


class IOHandlerGAMER(BaseIOHandler):
    _particle_reader = False
    _dataset_type = "gamer"

    def __init__(self, ds):
        super().__init__(ds)
        self._handle = ds._handle
        self._group_grid = ds._group_grid
        self._group_particle = ds._group_particle
        self._field_dtype = "float64"  # fixed even when FLOAT8 is off
        self._particle_handle = ds._particle_handle
        self.patch_size = ds.parameters["PatchSize"] * ds.refine_by
        self.pgroup = ds.refine_by**3  # number of patches in a patch group

    def _read_particle_coords(self, chunks, ptf):
        chunks = list(chunks)  # generator --> list
        p_idx = self.ds.index._particle_indices

        # shortcuts
        par_posx = self._group_particle["ParPosX"]
        par_posy = self._group_particle["ParPosY"]
        par_posz = self._group_particle["ParPosZ"]

        # currently GAMER does not support multiple particle types
        assert len(ptf) == 1
        ptype = list(ptf.keys())[0]

        for chunk in chunks:
            for g1, g2 in particle_sequences(chunk.objs):
                start = p_idx[g1.id]
                end = p_idx[g2.id + 1]
                x = np.asarray(par_posx[start:end], dtype=self._field_dtype)
                y = np.asarray(par_posy[start:end], dtype=self._field_dtype)
                z = np.asarray(par_posz[start:end], dtype=self._field_dtype)
                yield ptype, (x, y, z), 0.0

    def _read_particle_fields(self, chunks, ptf, selector):
        chunks = list(chunks)  # generator --> list
        p_idx = self.ds.index._particle_indices

        # shortcuts
        par_posx = self._group_particle["ParPosX"]
        par_posy = self._group_particle["ParPosY"]
        par_posz = self._group_particle["ParPosZ"]

        # currently GAMER does not support multiple particle types
        assert len(ptf) == 1
        ptype = list(ptf.keys())[0]
        pfields = ptf[ptype]

        for chunk in chunks:
            for g1, g2 in particle_sequences(chunk.objs):
                start = p_idx[g1.id]
                end = p_idx[g2.id + 1]
                x = np.asarray(par_posx[start:end], dtype=self._field_dtype)
                y = np.asarray(par_posy[start:end], dtype=self._field_dtype)
                z = np.asarray(par_posz[start:end], dtype=self._field_dtype)

                mask = selector.select_points(x, y, z, 0.0)
                if mask is None:
                    continue

                for field in pfields:
                    data = self._group_particle[field][start:end]
                    yield (ptype, field), data[mask]

    def _read_fluid_selection(self, chunks, selector, fields, size):
        chunks = list(chunks)  # generator --> list

        if any((ftype != "gamer" for ftype, fname in fields)):
            raise NotImplementedError

        rv = {}
        for field in fields:
            rv[field] = np.empty(size, dtype=self._field_dtype)

        ng = sum(len(c.objs) for c in chunks)  # c.objs is a list of grids
        mylog.debug(
            "Reading %s cells of %s fields in %s grids",
            size,
            [f2 for f1, f2 in fields],
            ng,
        )

        # shortcuts
        ps2 = self.patch_size
        ps1 = ps2 // 2

        for field in fields:
            ds = self._group_grid[field[1]]
            offset = 0
            for chunk in chunks:
                for gs in grid_sequences(chunk.objs):
                    start = (gs[0].id) * self.pgroup
                    end = (gs[-1].id + 1) * self.pgroup
                    buf = ds[start:end, :, :, :]
                    ngrid = len(gs)
                    data = np.empty((ngrid, ps2, ps2, ps2), dtype=self._field_dtype)

                    for g in range(ngrid):
                        pid0 = g * self.pgroup
                        data[g, 0:ps1, 0:ps1, 0:ps1] = buf[pid0 + 0, :, :, :]
                        data[g, 0:ps1, 0:ps1, ps1:ps2] = buf[pid0 + 1, :, :, :]
                        data[g, 0:ps1, ps1:ps2, 0:ps1] = buf[pid0 + 2, :, :, :]
                        data[g, ps1:ps2, 0:ps1, 0:ps1] = buf[pid0 + 3, :, :, :]
                        data[g, 0:ps1, ps1:ps2, ps1:ps2] = buf[pid0 + 4, :, :, :]
                        data[g, ps1:ps2, ps1:ps2, 0:ps1] = buf[pid0 + 5, :, :, :]
                        data[g, ps1:ps2, 0:ps1, ps1:ps2] = buf[pid0 + 6, :, :, :]
                        data[g, ps1:ps2, ps1:ps2, ps1:ps2] = buf[pid0 + 7, :, :, :]

                    data = data.transpose()

                    for i, g in enumerate(gs):
                        offset += g.select(selector, data[..., i], rv[field], offset)
        return rv

    def _read_chunk_data(self, chunk, fields):
        rv = {}
        if len(chunk.objs) == 0:
            return rv

        for g in chunk.objs:
            rv[g.id] = {}

        # Split into particles and non-particles
        fluid_fields, particle_fields = [], []
        for ftype, fname in fields:
            if ftype in self.ds.particle_types:
                particle_fields.append((ftype, fname))
            else:
                fluid_fields.append((ftype, fname))

        # particles
        if len(particle_fields) > 0:
            selector = AlwaysSelector(self.ds)
            rv.update(self._read_particle_selection([chunk], selector, particle_fields))

        # fluid
        if len(fluid_fields) == 0:
            return rv

        ps2 = self.patch_size
        ps1 = ps2 // 2

        for field in fluid_fields:
            ds = self._group_grid[field[1]]

            for gs in grid_sequences(chunk.objs):
                start = (gs[0].id) * self.pgroup
                end = (gs[-1].id + 1) * self.pgroup
                buf = ds[start:end, :, :, :]
                ngrid = len(gs)
                data = np.empty((ngrid, ps2, ps2, ps2), dtype=self._field_dtype)

                for g in range(ngrid):
                    pid0 = g * self.pgroup
                    data[g, 0:ps1, 0:ps1, 0:ps1] = buf[pid0 + 0, :, :, :]
                    data[g, 0:ps1, 0:ps1, ps1:ps2] = buf[pid0 + 1, :, :, :]
                    data[g, 0:ps1, ps1:ps2, 0:ps1] = buf[pid0 + 2, :, :, :]
                    data[g, ps1:ps2, 0:ps1, 0:ps1] = buf[pid0 + 3, :, :, :]
                    data[g, 0:ps1, ps1:ps2, ps1:ps2] = buf[pid0 + 4, :, :, :]
                    data[g, ps1:ps2, ps1:ps2, 0:ps1] = buf[pid0 + 5, :, :, :]
                    data[g, ps1:ps2, 0:ps1, ps1:ps2] = buf[pid0 + 6, :, :, :]
                    data[g, ps1:ps2, ps1:ps2, ps1:ps2] = buf[pid0 + 7, :, :, :]

                data = data.transpose()

                for i, g in enumerate(gs):
                    rv[g.id][field] = data[..., i]
        return rv
