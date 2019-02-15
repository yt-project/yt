"""
Enzo-P-specific IO functions



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2017, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.utilities.exceptions import \
    YTException
from yt.utilities.io_handler import \
    BaseIOHandler
from yt.extern.six import b, iteritems
from yt.utilities.on_demand_imports import _h5py as h5py
import numpy as np

_particle_position_names = {}

class EnzoPIOHandler(BaseIOHandler):

    _dataset_type = "enzo_p"
    _base = slice(None)
    _field_dtype = "float64"
    _sep = "_"

    def __init__(self, *args, **kwargs):
        super(EnzoPIOHandler, self).__init__(*args, **kwargs)
        self._base = self.ds.dimensionality * \
          (slice(self.ds.ghost_zones,
                 -self.ds.ghost_zones),)

    def _read_field_names(self, grid):
        if grid.filename is None: return []
        f = h5py.File(grid.filename, "r")
        try:
            group = f[grid.block_name]
        except KeyError:
            raise YTException(
                message="Grid %s is missing from data file %s." %
                (grid.block_name, grid.filename), ds=self.ds)
        fields = []
        ptypes = set()
        dtypes = set()
        # keep one field for each particle type so we can count later
        sample_pfields = {}
        for name, v in iteritems(group):
            if not hasattr(v, "shape") or v.dtype == "O":
                continue
            # mesh fields are "field <name>"
            if name.startswith("field"):
                _, fname = name.split(self._sep, 1)
                fields.append(("enzop", fname))
                dtypes.add(v.dtype)
            # particle fields are "particle <type> <name>"
            else:
                _, ftype, fname = name.split(self._sep, 2)
                fields.append((ftype, fname))
                ptypes.add(ftype)
                dtypes.add(v.dtype)
                if ftype not in sample_pfields:
                    sample_pfields[ftype] = fname
        self.sample_pfields = sample_pfields

        if len(dtypes) == 1:
            # Now, if everything we saw was the same dtype, we can go ahead and
            # set it here.  We do this because it is a HUGE savings for 32 bit
            # floats, since our numpy copying/casting is way faster than
            # h5py's, for some reason I don't understand.  This does *not* need
            # to be correct -- it will get fixed later -- it just needs to be
            # okay for now.
            self._field_dtype = list(dtypes)[0]
        f.close()
        return fields, ptypes

    def _read_particle_coords(self, chunks, ptf):
        for rv in self._read_particle_fields(chunks, ptf, None):
            yield rv

    def _read_particle_fields(self, chunks, ptf, selector):
        chunks = list(chunks)
        dc = self.ds.domain_center.in_units("code_length").d
        for chunk in chunks: # These should be organized by grid filename
            f = None
            for g in chunk.objs:
                if g.filename is None:
                    continue
                if f is None:
                    f = h5py.File(g.filename, "r")
                if g.particle_count is None:
                    fnstr = "%s/%s" % \
                      (g.block_name, self._sep.join(["particle", "%s", "%s"]))
                    g.particle_count = \
                      dict((ptype, f.get(fnstr %
                            (ptype, self.sample_pfields[ptype])).size)
                            for ptype in self.sample_pfields)
                    g.total_particles = sum(g.particle_count.values())
                if g.total_particles == 0:
                    continue
                group = f.get(g.block_name)
                for ptype, field_list in sorted(ptf.items()):
                    pn = self._sep.join(
                        ["particle", ptype, "%s"])
                    if g.particle_count[ptype] == 0:
                        continue
                    coords = \
                      tuple(np.asarray(group.get(pn % ax)[()], dtype="=f8")
                            for ax in 'xyz'[:self.ds.dimensionality])
                    for i in range(self.ds.dimensionality, 3):
                        coords += \
                          (dc[i] * np.ones(g.particle_count[ptype], dtype="f8"),)
                    if selector is None:
                        # This only ever happens if the call is made from
                        # _read_particle_coords.
                        yield ptype, coords
                        continue
                    coords += (0.0,)
                    mask = selector.select_points(*coords)
                    if mask is None:
                        continue
                    for field in field_list:
                        data = np.asarray(group.get(pn % field)[()], "=f8")
                        yield (ptype, field), data[mask]
            if f: f.close()

    def io_iter(self, chunks, fields):
        for chunk in chunks:
            fid = None
            filename = -1
            for obj in chunk.objs:
                if obj.filename is None: continue
                if obj.filename != filename:
                    # Note one really important thing here: even if we do
                    # implement LRU caching in the _read_obj_field function,
                    # we'll still be doing file opening and whatnot.  This is a
                    # problem, but one we can return to.
                    if fid is not None:
                        fid.close()
                    fid = h5py.h5f.open(b(obj.filename), h5py.h5f.ACC_RDONLY)
                    filename = obj.filename
                for field in fields:
                    data = None
                    yield field, obj, self._read_obj_field(
                        obj, field, (fid, data))
        if fid is not None:
            fid.close()

    def _read_obj_field(self, obj, field, fid_data):
        if fid_data is None: fid_data = (None, None)
        fid, data = fid_data
        if fid is None:
            close = True
            fid = h5py.h5f.open(b(obj.filename), h5py.h5f.ACC_RDONLY)
        else:
            close = False
        ftype, fname = field
        node = "/%s/field%s%s" % (obj.block_name, self._sep, fname)
        dg = h5py.h5d.open(fid, b(node))
        rdata = np.empty(self.ds.grid_dimensions[:self.ds.dimensionality][::-1],
                         dtype=self._field_dtype)
        dg.read(h5py.h5s.ALL, h5py.h5s.ALL, rdata)
        if close:
            fid.close()
        data = rdata[self._base].T
        if self.ds.dimensionality < 3:
            nshape = data.shape + (1,)*(3 - self.ds.dimensionality)
            data  = np.reshape(data, nshape)
        return data
