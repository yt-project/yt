from typing import Dict

import numpy as np

from yt.geometry.selection_routines import GridSelector
from yt.utilities.io_handler import BaseIOHandler
from yt.utilities.logger import ytLogger as mylog
from yt.utilities.on_demand_imports import _h5py as h5py

_convert_mass = ("particle_mass", "mass")

_particle_position_names: Dict[str, str] = {}


class IOHandlerPackedHDF5(BaseIOHandler):

    _dataset_type = "enzo_packed_3d"
    _base = slice(None)
    _field_dtype = "float64"

    def _read_field_names(self, grid):
        if grid.filename is None:
            return []
        f = h5py.File(grid.filename, mode="r")
        try:
            group = f["/Grid%08i" % grid.id]
        except KeyError:
            group = f
        fields = []
        dtypes = set()
        add_io = "io" in grid.ds.particle_types
        add_dm = "DarkMatter" in grid.ds.particle_types
        for name, v in group.items():
            # NOTE: This won't work with 1D datasets or references.
            # For all versions of Enzo I know about, we can assume all floats
            # are of the same size.  So, let's grab one.
            if not hasattr(v, "shape") or v.dtype == "O":
                continue
            elif len(v.dims) == 1:
                if grid.ds.dimensionality == 1:
                    fields.append(("enzo", str(name)))
                elif add_io:
                    fields.append(("io", str(name)))
                elif add_dm:
                    fields.append(("DarkMatter", str(name)))
            else:
                fields.append(("enzo", str(name)))
                dtypes.add(v.dtype)

        if len(dtypes) == 1:
            # Now, if everything we saw was the same dtype, we can go ahead and
            # set it here.  We do this because it is a HUGE savings for 32 bit
            # floats, since our numpy copying/casting is way faster than
            # h5py's, for some reason I don't understand.  This does *not* need
            # to be correct -- it will get fixed later -- it just needs to be
            # okay for now.
            self._field_dtype = list(dtypes)[0]
        f.close()
        return fields

    @property
    def _read_exception(self):
        return (KeyError,)

    def _read_particle_coords(self, chunks, ptf):
        yield from (
            (ptype, xyz, 0.0)
            for ptype, xyz in self._read_particle_fields(chunks, ptf, None)
        )

    def _read_particle_fields(self, chunks, ptf, selector):
        chunks = list(chunks)
        for chunk in chunks:  # These should be organized by grid filename
            f = None
            for g in chunk.objs:
                if g.filename is None:
                    continue
                if f is None:
                    # print("Opening (read) %s" % g.filename)
                    f = h5py.File(g.filename, mode="r")
                nap = sum(g.NumberOfActiveParticles.values())
                if g.NumberOfParticles == 0 and nap == 0:
                    continue
                ds = f.get("/Grid%08i" % g.id)
                for ptype, field_list in sorted(ptf.items()):
                    if ptype == "io":
                        if g.NumberOfParticles == 0:
                            continue
                        pds = ds
                    elif ptype == "DarkMatter":
                        if g.NumberOfActiveParticles[ptype] == 0:
                            continue
                        pds = ds
                    elif not g.NumberOfActiveParticles[ptype]:
                        continue
                    else:
                        for pname in ["Active Particles", "Particles"]:
                            pds = ds.get(f"{pname}/{ptype}")
                            if pds is not None:
                                break
                        else:
                            raise RuntimeError(
                                "Could not find active particle group in data."
                            )
                    pn = _particle_position_names.get(ptype, r"particle_position_%s")
                    x, y, z = (
                        np.asarray(pds.get(pn % ax)[()], dtype="=f8") for ax in "xyz"
                    )
                    if selector is None:
                        # This only ever happens if the call is made from
                        # _read_particle_coords.
                        yield ptype, (x, y, z)
                        continue
                    mask = selector.select_points(x, y, z, 0.0)
                    if mask is None:
                        continue
                    for field in field_list:
                        data = np.asarray(pds.get(field)[()], "=f8")
                        if field in _convert_mass:
                            data *= g.dds.prod(dtype="f8")
                        yield (ptype, field), data[mask]
            if f:
                f.close()

    def io_iter(self, chunks, fields):
        h5_dtype = self._field_dtype
        for chunk in chunks:
            fid = None
            filename = -1
            for obj in chunk.objs:
                if obj.filename is None:
                    continue
                if obj.filename != filename:
                    # Note one really important thing here: even if we do
                    # implement LRU caching in the _read_obj_field function,
                    # we'll still be doing file opening and whatnot.  This is a
                    # problem, but one we can return to.
                    if fid is not None:
                        fid.close()
                    fid = h5py.h5f.open(
                        obj.filename.encode("latin-1"), h5py.h5f.ACC_RDONLY
                    )
                    filename = obj.filename
                for field in fields:
                    nodal_flag = self.ds.field_info[field].nodal_flag
                    dims = obj.ActiveDimensions[::-1] + nodal_flag[::-1]
                    data = np.empty(dims, dtype=h5_dtype)
                    yield field, obj, self._read_obj_field(obj, field, (fid, data))
        if fid is not None:
            fid.close()

    def _read_obj_field(self, obj, field, fid_data):
        if fid_data is None:
            fid_data = (None, None)
        fid, data = fid_data
        if fid is None:
            close = True
            fid = h5py.h5f.open(obj.filename.encode("latin-1"), h5py.h5f.ACC_RDONLY)
        else:
            close = False
        if data is None:
            data = np.empty(obj.ActiveDimensions[::-1], dtype=self._field_dtype)
        ftype, fname = field
        try:
            node = "/Grid%08i/%s" % (obj.id, fname)
            dg = h5py.h5d.open(fid, node.encode("latin-1"))
        except KeyError:
            if fname == "Dark_Matter_Density":
                data[:] = 0
                return data.T
            raise
        dg.read(h5py.h5s.ALL, h5py.h5s.ALL, data)
        # I don't know why, but on some installations of h5py this works, but
        # on others, nope.  Doesn't seem to be a version thing.
        # dg.close()
        if close:
            fid.close()
        return data.T


class IOHandlerPackedHDF5GhostZones(IOHandlerPackedHDF5):
    _dataset_type = "enzo_packed_3d_gz"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        NGZ = self.ds.parameters.get("NumberOfGhostZones", 3)
        self._base = (slice(NGZ, -NGZ), slice(NGZ, -NGZ), slice(NGZ, -NGZ))

    def _read_obj_field(self, *args, **kwargs):
        return super()._read_obj_field(*args, **kwargs)[self._base]


class IOHandlerInMemory(BaseIOHandler):

    _dataset_type = "enzo_inline"

    def __init__(self, ds, ghost_zones=3):
        self.ds = ds
        import enzo

        self.enzo = enzo
        self.grids_in_memory = enzo.grid_data
        self.old_grids_in_memory = enzo.old_grid_data
        self.my_slice = (
            slice(ghost_zones, -ghost_zones),
            slice(ghost_zones, -ghost_zones),
            slice(ghost_zones, -ghost_zones),
        )
        BaseIOHandler.__init__(self, ds)

    def _read_field_names(self, grid):
        fields = []
        add_io = "io" in grid.ds.particle_types
        for name, v in self.grids_in_memory[grid.id].items():
            # NOTE: This won't work with 1D datasets or references.
            if not hasattr(v, "shape") or v.dtype == "O":
                continue
            elif v.ndim == 1:
                if grid.ds.dimensionality == 1:
                    fields.append(("enzo", str(name)))
                elif add_io:
                    fields.append(("io", str(name)))
            else:
                fields.append(("enzo", str(name)))
        return fields

    def _read_fluid_selection(self, chunks, selector, fields, size):
        rv = {}
        # Now we have to do something unpleasant
        chunks = list(chunks)
        if isinstance(selector, GridSelector):
            if not (len(chunks) == len(chunks[0].objs) == 1):
                raise RuntimeError
            g = chunks[0].objs[0]
            for ftype, fname in fields:
                rv[(ftype, fname)] = self.grids_in_memory[g.id][fname].swapaxes(0, 2)
            return rv
        if size is None:
            size = sum(g.count(selector) for chunk in chunks for g in chunk.objs)
        for field in fields:
            ftype, fname = field
            fsize = size
            rv[field] = np.empty(fsize, dtype="float64")
        ng = sum(len(c.objs) for c in chunks)
        mylog.debug(
            "Reading %s cells of %s fields in %s grids",
            size,
            [f2 for f1, f2 in fields],
            ng,
        )

        ind = 0
        for chunk in chunks:
            for g in chunk.objs:
                # We want a *hard error* here.
                # if g.id not in self.grids_in_memory: continue
                for field in fields:
                    ftype, fname = field
                    data_view = self.grids_in_memory[g.id][fname][
                        self.my_slice
                    ].swapaxes(0, 2)
                    nd = g.select(selector, data_view, rv[field], ind)
                ind += nd
        assert ind == fsize
        return rv

    def _read_particle_coords(self, chunks, ptf):
        chunks = list(chunks)
        for chunk in chunks:  # These should be organized by grid filename
            for g in chunk.objs:
                if g.id not in self.grids_in_memory:
                    continue
                nap = sum(g.NumberOfActiveParticles.values())
                if g.NumberOfParticles == 0 and nap == 0:
                    continue
                for ptype in sorted(ptf):
                    x, y, z = (
                        self.grids_in_memory[g.id]["particle_position_x"],
                        self.grids_in_memory[g.id]["particle_position_y"],
                        self.grids_in_memory[g.id]["particle_position_z"],
                    )
                    yield ptype, (x, y, z), 0.0

    def _read_particle_fields(self, chunks, ptf, selector):
        chunks = list(chunks)
        for chunk in chunks:  # These should be organized by grid filename
            for g in chunk.objs:
                if g.id not in self.grids_in_memory:
                    continue
                nap = sum(g.NumberOfActiveParticles.values())
                if g.NumberOfParticles == 0 and nap == 0:
                    continue
                for ptype, field_list in sorted(ptf.items()):
                    x, y, z = (
                        self.grids_in_memory[g.id]["particle_position_x"],
                        self.grids_in_memory[g.id]["particle_position_y"],
                        self.grids_in_memory[g.id]["particle_position_z"],
                    )
                    mask = selector.select_points(x, y, z, 0.0)
                    if mask is None:
                        continue
                    for field in field_list:
                        data = self.grids_in_memory[g.id][field]
                        if field in _convert_mass:
                            data = data * g.dds.prod(dtype="f8")
                        yield (ptype, field), data[mask]


class IOHandlerPacked2D(IOHandlerPackedHDF5):

    _dataset_type = "enzo_packed_2d"
    _particle_reader = False

    def _read_data_set(self, grid, field):
        f = h5py.File(grid.filename, mode="r")
        ds = f["/Grid%08i/%s" % (grid.id, field)][:]
        f.close()
        return ds.transpose()[:, :, None]

    def _read_fluid_selection(self, chunks, selector, fields, size):
        rv = {}
        # Now we have to do something unpleasant
        chunks = list(chunks)
        if isinstance(selector, GridSelector):
            if not (len(chunks) == len(chunks[0].objs) == 1):
                raise RuntimeError
            g = chunks[0].objs[0]
            f = h5py.File(g.filename, mode="r")
            gds = f.get("/Grid%08i" % g.id)
            for ftype, fname in fields:
                rv[(ftype, fname)] = np.atleast_3d(gds.get(fname)[()].transpose())
            f.close()
            return rv
        if size is None:
            size = sum(g.count(selector) for chunk in chunks for g in chunk.objs)
        for field in fields:
            ftype, fname = field
            fsize = size
            rv[field] = np.empty(fsize, dtype="float64")
        ng = sum(len(c.objs) for c in chunks)
        mylog.debug(
            "Reading %s cells of %s fields in %s grids",
            size,
            [f2 for f1, f2 in fields],
            ng,
        )
        ind = 0
        for chunk in chunks:
            f = None
            for g in chunk.objs:
                if f is None:
                    # print("Opening (count) %s" % g.filename)
                    f = h5py.File(g.filename, mode="r")
                gds = f.get("/Grid%08i" % g.id)
                if gds is None:
                    gds = f
                for field in fields:
                    ftype, fname = field
                    ds = np.atleast_3d(gds.get(fname)[()].transpose())
                    nd = g.select(selector, ds, rv[field], ind)  # caches
                ind += nd
            f.close()
        return rv


class IOHandlerPacked1D(IOHandlerPackedHDF5):

    _dataset_type = "enzo_packed_1d"
    _particle_reader = False

    def _read_data_set(self, grid, field):
        f = h5py.File(grid.filename, mode="r")
        ds = f["/Grid%08i/%s" % (grid.id, field)][:]
        f.close()
        return ds.transpose()[:, None, None]
