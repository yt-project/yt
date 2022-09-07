import os
from collections import defaultdict

import numpy as np

from yt.frontends.chombo.io import parse_orion_sinks
from yt.funcs import mylog
from yt.geometry.selection_routines import GridSelector
from yt.utilities.io_handler import BaseIOHandler


def _remove_raw(all_fields, raw_fields):
    centered_fields = set(all_fields)
    for raw in raw_fields:
        centered_fields.discard(raw)
    return list(centered_fields)


class IOHandlerBoxlib(BaseIOHandler):

    _dataset_type = "boxlib_native"

    def __init__(self, ds, *args, **kwargs):
        super().__init__(ds)

    def _read_fluid_selection(self, chunks, selector, fields, size):
        chunks = list(chunks)
        if any((not (ftype == "boxlib" or ftype == "raw") for ftype, fname in fields)):
            raise NotImplementedError
        rv = {}
        raw_fields = []
        for field in fields:
            if field[0] == "raw":
                nodal_flag = self.ds.nodal_flags[field[1]]
                num_nodes = 2 ** sum(nodal_flag)
                rv[field] = np.empty((size, num_nodes), dtype="float64")
                raw_fields.append(field)
            else:
                rv[field] = np.empty(size, dtype="float64")
        centered_fields = _remove_raw(fields, raw_fields)
        ng = sum(len(c.objs) for c in chunks)
        mylog.debug(
            "Reading %s cells of %s fields in %s grids",
            size,
            [f2 for f1, f2 in fields],
            ng,
        )
        ind = 0
        for chunk in chunks:
            data = self._read_chunk_data(chunk, centered_fields)
            for g in chunk.objs:
                for field in fields:
                    if field in centered_fields:
                        ds = data[g.id].pop(field)
                    else:
                        ds = self._read_raw_field(g, field)
                    nd = g.select(selector, ds, rv[field], ind)
                ind += nd
                data.pop(g.id)
        return rv

    def _read_raw_field(self, grid, field):
        field_name = field[1]
        base_dir = self.ds.index.raw_file

        nghost = self.ds.index.raw_field_nghost[field_name]
        box_list = self.ds.index.raw_field_map[field_name][0]
        fn_list = self.ds.index.raw_field_map[field_name][1]
        offset_list = self.ds.index.raw_field_map[field_name][2]

        lev = grid.Level
        filename = os.path.join(base_dir, f"Level_{lev}", fn_list[grid.id])
        offset = offset_list[grid.id]
        box = box_list[grid.id]

        lo = box[0] - nghost
        hi = box[1] + nghost
        shape = hi - lo + 1
        with open(filename, "rb") as f:
            f.seek(offset)
            f.readline()  # always skip the first line
            arr = np.fromfile(f, "float64", np.product(shape))
            arr = arr.reshape(shape, order="F")
        return arr[
            tuple(
                slice(None) if (nghost[dim] == 0) else slice(nghost[dim], -nghost[dim])
                for dim in range(self.ds.dimensionality)
            )
        ]

    def _read_chunk_data(self, chunk, fields):
        data = {}
        grids_by_file = defaultdict(list)
        if len(chunk.objs) == 0:
            return data
        for g in chunk.objs:
            if g.filename is None:
                continue
            grids_by_file[g.filename].append(g)
        dtype = self.ds.index._dtype
        bpr = dtype.itemsize
        for filename in grids_by_file:
            grids = grids_by_file[filename]
            grids.sort(key=lambda a: a._offset)
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
                        v = v.reshape(grid.ActiveDimensions, order="F")
                        data[grid.id][field] = v
                        local_offset = 0
                    else:
                        local_offset += size
            f.close()
        return data

    def _read_particle_coords(self, chunks, ptf):
        yield from (
            (ptype, xyz, 0.0)
            for ptype, xyz in self._read_particle_fields(chunks, ptf, None)
        )

    def _read_particle_fields(self, chunks, ptf, selector):
        for chunk in chunks:  # These should be organized by grid filename
            for g in chunk.objs:
                for ptype, field_list in sorted(ptf.items()):
                    npart = g._pdata[ptype]["NumberOfParticles"]
                    if npart == 0:
                        continue

                    fn = g._pdata[ptype]["particle_filename"]
                    offset = g._pdata[ptype]["offset"]
                    pheader = self.ds.index.particle_headers[ptype]

                    with open(fn, "rb") as f:
                        # read in the position fields for selection
                        f.seek(offset + pheader.particle_int_dtype.itemsize * npart)
                        rdata = np.fromfile(
                            f, pheader.real_type, pheader.num_real * npart
                        )

                        # Allow reading particles in 1, 2, and 3 dimensions,
                        # setting the appropriate default for unused dimensions.
                        pos = []
                        for idim in [1, 2, 3]:
                            if g.ds.dimensionality >= idim:
                                pos.append(
                                    np.asarray(
                                        rdata[idim - 1 :: pheader.num_real],
                                        dtype=np.float64,
                                    )
                                )
                            else:
                                center = 0.5 * (
                                    g.LeftEdge[idim - 1] + g.RightEdge[idim - 1]
                                )
                                pos.append(np.full(npart, center, dtype=np.float64))
                        x, y, z = pos

                        if selector is None:
                            # This only ever happens if the call is made from
                            # _read_particle_coords.
                            yield ptype, (x, y, z)
                            continue
                        mask = selector.select_points(x, y, z, 0.0)
                        if mask is None:
                            continue
                        for field in field_list:
                            # handle the case that this is an integer field
                            int_fnames = [
                                fname for _, fname in pheader.known_int_fields
                            ]
                            if field in int_fnames:
                                ind = int_fnames.index(field)
                                f.seek(offset)
                                idata = np.fromfile(
                                    f, pheader.int_type, pheader.num_int * npart
                                )
                                data = np.asarray(
                                    idata[ind :: pheader.num_int], dtype=np.float64
                                )
                                yield (ptype, field), data[mask].flatten()

                            # handle case that this is a real field
                            real_fnames = [
                                fname for _, fname in pheader.known_real_fields
                            ]
                            if field in real_fnames:
                                ind = real_fnames.index(field)
                                data = np.asarray(
                                    rdata[ind :: pheader.num_real], dtype=np.float64
                                )
                                yield (ptype, field), data[mask].flatten()


class IOHandlerOrion(IOHandlerBoxlib):
    _dataset_type = "orion_native"

    _particle_filename = None

    @property
    def particle_filename(self):
        fn = os.path.join(self.ds.output_dir, "StarParticles")
        if not os.path.exists(fn):
            fn = os.path.join(self.ds.output_dir, "SinkParticles")
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
            entry = line.strip().split(" ")[self.particle_field_index[field]]
            return float(entry)

        try:
            lines = self._cached_lines
            for num in grid._particle_line_numbers:
                line = lines[num]
                particles.append(read(line, field))
            return np.array(particles)
        except AttributeError:
            fn = self.particle_filename
            with open(fn) as f:
                lines = f.readlines()
                self._cached_lines = lines
                for num in grid._particle_line_numbers:
                    line = lines[num]
                    particles.append(read(line, field))
            return np.array(particles)
