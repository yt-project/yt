import numpy as np

from yt.utilities.exceptions import YTDomainOverflow
from yt.utilities.io_handler import BaseIOHandler, BaseParticleIOHandler
from yt.utilities.logger import ytLogger as mylog


class IOHandlerStream(BaseIOHandler):

    _dataset_type = "stream"
    _vector_fields = {"particle_velocity": 3, "particle_position": 3}

    def __init__(self, ds):
        self.fields = ds.stream_handler.fields
        self.field_units = ds.stream_handler.field_units
        super().__init__(ds)

    def _read_data_set(self, grid, field):
        # This is where we implement processor-locking
        tr = self.fields[grid.id][field]
        if callable(tr):
            tr = tr(grid, field)
        # If it's particles, we copy.
        if len(tr.shape) == 1:
            return tr.copy()
        # New in-place unit conversion breaks if we don't copy first
        return tr

    def _read_fluid_selection(self, chunks, selector, fields, size):
        chunks = list(chunks)
        if any((ftype not in self.ds.fluid_types for ftype, fname in fields)):
            raise NotImplementedError
        rv = {}
        for field in fields:
            rv[field] = self.ds.arr(np.empty(size, dtype="float64"))

        ng = sum(len(c.objs) for c in chunks)
        mylog.debug(
            "Reading %s cells of %s fields in %s blocks",
            size,
            [f2 for f1, f2 in fields],
            ng,
        )
        for field in fields:
            ftype, fname = field
            ind = 0
            for chunk in chunks:
                for g in chunk.objs:
                    ds = self.fields[g.id][ftype, fname]
                    if callable(ds):
                        ds = ds(g, field)
                    ind += g.select(selector, ds, rv[field], ind)  # caches
        return rv

    def _read_particle_coords(self, chunks, ptf):
        chunks = list(chunks)
        for chunk in chunks:
            for g in chunk.objs:
                if g.NumberOfParticles == 0:
                    continue
                gf = self.fields[g.id]
                for ptype in sorted(ptf):
                    if (ptype, "particle_position") in gf:
                        x, y, z = gf[ptype, "particle_position"].T
                    else:
                        x, y, z = (
                            gf[ptype, f"particle_position_{ax}"]
                            for ax in self.ds.coordinates.axis_order
                        )
                    yield ptype, (x, y, z), 0.0

    def _read_particle_fields(self, chunks, ptf, selector):
        chunks = list(chunks)
        for chunk in chunks:
            for g in chunk.objs:
                if g.NumberOfParticles == 0:
                    continue
                gf = self.fields[g.id]
                for ptype, field_list in sorted(ptf.items()):
                    if (ptype, "particle_position") in gf:
                        x, y, z = gf[ptype, "particle_position"].T
                    else:
                        x, y, z = (
                            gf[ptype, f"particle_position_{ax}"]
                            for ax in self.ds.coordinates.axis_order
                        )
                    mask = selector.select_points(x, y, z, 0.0)
                    if mask is None:
                        continue
                    for field in field_list:
                        data = np.asarray(gf[ptype, field])
                        yield (ptype, field), data[mask]

    @property
    def _read_exception(self):
        return KeyError


class StreamParticleIOHandler(BaseParticleIOHandler):
    _dataset_type = "stream_particles"
    _vector_fields = {"particle_velocity": 3, "particle_position": 3}

    def __init__(self, ds):
        self.fields = ds.stream_handler.fields
        super().__init__(ds)

    def _read_particle_coords(self, chunks, ptf):
        for data_file in sorted(
            self._get_data_files(chunks), key=lambda x: (x.filename, x.start)
        ):
            f = self.fields[data_file.filename]
            # This double-reads
            for ptype in sorted(ptf):
                yield ptype, (
                    f[ptype, "particle_position_x"],
                    f[ptype, "particle_position_y"],
                    f[ptype, "particle_position_z"],
                ), 0.0

    def _read_smoothing_length(self, chunks, ptf, ptype):
        for data_file in sorted(
            self._get_data_files(chunks), key=lambda x: (x.filename, x.start)
        ):
            f = self.fields[data_file.filename]
            return f[ptype, "smoothing_length"]

    def _get_data_files(self, chunks):
        data_files = set()
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        return data_files

    def _read_particle_data_file(self, data_file, ptf, selector=None):

        return_data = {}
        f = self.fields[data_file.filename]
        for ptype, field_list in sorted(ptf.items()):
            if (ptype, "particle_position") in f:
                ppos = f[ptype, "particle_position"]
                x = ppos[:, 0]
                y = ppos[:, 1]
                z = ppos[:, 2]
            else:
                x, y, z = (f[ptype, f"particle_position_{ax}"] for ax in "xyz")
            if (ptype, "smoothing_length") in self.ds.field_list:
                hsml = f[ptype, "smoothing_length"]
            else:
                hsml = 0.0

            if selector:
                mask = selector.select_points(x, y, z, hsml)
            if mask is None:
                continue
            for field in field_list:
                data = f[ptype, field]
                if selector:
                    data = data[mask]

                return_data[(ptype, field)] = data

        return return_data

    def _yield_coordinates(self, data_file, needed_ptype=None):
        # self.fields[g.id][fname] is the pattern here
        for ptype in self.ds.particle_types_raw:
            if needed_ptype is not None and needed_ptype is not ptype:
                continue
            try:
                pos = np.column_stack(
                    [
                        self.fields[data_file.filename][
                            (ptype, f"particle_position_{ax}")
                        ]
                        for ax in "xyz"
                    ]
                )
            except KeyError:
                pos = self.fields[data_file.filename][ptype, "particle_position"]
            if np.any(pos.min(axis=0) < data_file.ds.domain_left_edge) or np.any(
                pos.max(axis=0) > data_file.ds.domain_right_edge
            ):
                raise YTDomainOverflow(
                    pos.min(axis=0),
                    pos.max(axis=0),
                    data_file.ds.domain_left_edge,
                    data_file.ds.domain_right_edge,
                )
            yield ptype, pos

    def _get_smoothing_length(self, data_file, dtype, shape):
        ptype = self.ds._sph_ptypes[0]
        return self.fields[data_file.filename][ptype, "smoothing_length"]

    def _count_particles(self, data_file):
        pcount = {}
        for ptype in self.ds.particle_types_raw:
            pcount[ptype] = 0
        # stream datasets only have one "file"
        if data_file.file_id > 0:
            return pcount
        for ptype in self.ds.particle_types_raw:
            d = self.fields[data_file.filename]
            try:
                pcount[ptype] = d[ptype, "particle_position_x"].size
            except KeyError:
                pcount[ptype] = d[ptype, "particle_position"].shape[0]
        return pcount

    def _identify_fields(self, data_file):
        return self.fields[data_file.filename].keys(), {}


class IOHandlerStreamHexahedral(BaseIOHandler):
    _dataset_type = "stream_hexahedral"
    _vector_fields = {"particle_velocity": 3, "particle_position": 3}

    def __init__(self, ds):
        self.fields = ds.stream_handler.fields
        super().__init__(ds)

    def _read_fluid_selection(self, chunks, selector, fields, size):
        chunks = list(chunks)
        assert len(chunks) == 1
        chunk = chunks[0]
        rv = {}
        for field in fields:
            ftype, fname = field
            rv[field] = np.empty(size, dtype="float64")
        ngrids = sum(len(chunk.objs) for chunk in chunks)
        mylog.debug(
            "Reading %s cells of %s fields in %s blocks",
            size,
            [fn for ft, fn in fields],
            ngrids,
        )
        for field in fields:
            ind = 0
            ftype, fname = field
            for chunk in chunks:
                for g in chunk.objs:
                    ds = self.fields[g.mesh_id].get(field, None)
                    if ds is None:
                        ds = self.fields[g.mesh_id][fname]
                    ind += g.select(selector, ds, rv[field], ind)  # caches
        return rv


class IOHandlerStreamOctree(BaseIOHandler):
    _dataset_type = "stream_octree"
    _vector_fields = {"particle_velocity": 3, "particle_position": 3}

    def __init__(self, ds):
        self.fields = ds.stream_handler.fields
        super().__init__(ds)

    def _read_fluid_selection(self, chunks, selector, fields, size):
        rv = {}
        ind = 0
        chunks = list(chunks)
        assert len(chunks) == 1
        for chunk in chunks:
            assert len(chunk.objs) == 1
            for subset in chunk.objs:
                field_vals = {}
                for field in fields:
                    field_vals[field] = self.fields[
                        subset.domain_id - subset._domain_offset
                    ][field]
                subset.fill(field_vals, rv, selector, ind)
        return rv


class IOHandlerStreamUnstructured(BaseIOHandler):
    _dataset_type = "stream_unstructured"

    def __init__(self, ds):
        self.fields = ds.stream_handler.fields
        super().__init__(ds)

    def _read_fluid_selection(self, chunks, selector, fields, size):
        chunks = list(chunks)
        rv = {}
        for field in fields:
            ftype, fname = field
            if ftype == "all":
                ci = np.concatenate(
                    [mesh.connectivity_indices for mesh in self.ds.index.mesh_union]
                )
            else:
                mesh_id = int(ftype[-1]) - 1
                m = self.ds.index.meshes[mesh_id]
                ci = m.connectivity_indices
            num_elem = ci.shape[0]
            if fname in self.ds._node_fields:
                nodes_per_element = ci.shape[1]
                rv[field] = np.empty((num_elem, nodes_per_element), dtype="float64")
            else:
                rv[field] = np.empty(num_elem, dtype="float64")
        for field in fields:
            ind = 0
            ftype, fname = field
            if ftype == "all":
                objs = [mesh for mesh in self.ds.index.mesh_union]
            else:
                mesh_ids = [int(ftype[-1])]
                chunk = chunks[mesh_ids[0] - 1]
                objs = chunk.objs
            for g in objs:
                ds = self.fields[g.mesh_id].get(field, None)
                if ds is None:
                    f = ("connect%d" % (g.mesh_id + 1), fname)
                    ds = self.fields[g.mesh_id][f]
                ind += g.select(selector, ds, rv[field], ind)  # caches
            rv[field] = rv[field][:ind]
        return rv
