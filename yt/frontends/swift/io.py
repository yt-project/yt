"""
OWLS data-file handling function




"""
from __future__ import print_function

#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.utilities.on_demand_imports import _h5py as h5py
import numpy as np

from yt.frontends.sph.io import \
    IOHandlerSPH

class IOHandlerSwift(IOHandlerSPH):
    _dataset_type = "swift"

    def __init__(self, ds, *args, **kwargs):
        super(IOHandlerSwift, self).__init__(ds, *args, **kwargs)

    def _read_fluid_selection(self, chunks, selector, fields, size):
        raise NotImplementedError

    def _read_particle_coords(self, chunks, ptf):
        # This will read chunks and yield the results.
        chunks = list(chunks)
        data_files = set([])
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        for data_file in sorted(data_files, key=lambda x: x.filename):
            si, ei = data_file.start, data_file.end
            f = h5py.File(data_file.filename, "r")
            # This double-reads
            for ptype, field_list in sorted(ptf.items()):
                if data_file.total_particles[ptype] == 0:
                    continue
                x = f["/%s/Coordinates" % ptype][si:ei, 0].astype("float64")
                y = f["/%s/Coordinates" % ptype][si:ei, 1].astype("float64")
                z = f["/%s/Coordinates" % ptype][si:ei, 2].astype("float64")
                if ptype == self.ds._sph_ptype:
                    hsml = f[
                        "/%s/SmoothingLength" % ptype][si:ei].astype("float64")
                else:
                    hsml = 0.0
                yield ptype, (x, y, z), hsml
            f.close()

    def _yield_coordinates(self, data_file, needed_ptype=None):
        si, ei = data_file.start, data_file.end
        f = h5py.File(data_file.filename)
        pcount = f["/Header"].attrs["NumPart_ThisFile"][:].astype("int")
        np.clip(pcount - si, 0, ei - si, out=pcount)
        pcount = pcount.sum()
        for key in f.keys():
            if not key.startswith("PartType"):
                continue
            if "Coordinates" not in f[key]:
                continue
            if needed_ptype and key != needed_ptype:
                continue
            ds = f[key]["Coordinates"][si:ei,...]
            dt = ds.dtype.newbyteorder("N") # Native
            pos = np.empty(ds.shape, dtype=dt)
            pos[:] = ds
            yield key, pos
        f.close()

    def _get_smoothing_length(self, data_file, position_dtype, position_shape):
        ptype = self.ds._sph_ptype
        ind = int(ptype[-1])
        si, ei = data_file.start, data_file.end
        with h5py.File(data_file.filename) as f:
            pcount = f["/Header"].attrs["NumPart_ThisFile"][ind].astype("int")
            pcount = np.clip(pcount - si, 0, ei - si)
            ds = f[ptype]["SmoothingLength"][si:ei,...]
            dt = ds.dtype.newbyteorder("N") # Native
            if position_dtype is not None and dt < position_dtype:
                # Sometimes positions are stored in double precision
                # but smoothing lengths are stored in single precision.
                # In these cases upcast smoothing length to double precision
                # to avoid ValueErrors when we pass these arrays to Cython.
                dt = position_dtype
            hsml = np.empty(ds.shape, dtype=dt)
            hsml[:] = ds
            return hsml

    def _read_particle_fields(self, chunks, ptf, selector):
        # Now we have all the sizes, and we can allocate
        data_files = set([])
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)

        for data_file in sorted(data_files, key=lambda x: x.filename):
            si, ei = data_file.start, data_file.end
            f = h5py.File(data_file.filename, "r")
            for ptype, field_list in sorted(ptf.items()):
                if data_file.total_particles[ptype] == 0:
                    continue
                g = f["/%s" % ptype]
                coords = g["Coordinates"][si:ei].astype("float64")
                if ptype == 'PartType0':
                    hsmls = g["SmoothingLength"][si:ei].astype("float64")
                else:
                    hsmls = 0.0
                mask = selector.select_points(
                    coords[:,0], coords[:,1], coords[:,2], hsmls)
                del coords
                if mask is None:
                    continue
                for field in field_list:
                    if field in ("Mass", "Masses"):
                        data = g[self.ds._particle_mass_name][si:ei][mask, ...]
                    else:
                        data = g[field][si:ei][mask, ...]

                    yield (ptype, field), data
            f.close()

    def _count_particles(self, data_file):
        si, ei = data_file.start, data_file.end
        f = h5py.File(data_file.filename, "r")
        pcount = f["/Header"].attrs["NumPart_ThisFile"][:].astype("int")
        f.close()
        if None not in (si, ei):
            np.clip(pcount - si, 0, ei - si, out=pcount)
        npart = dict(("PartType%s" % (i), v) for i, v in enumerate(pcount))
        return npart

    def _identify_fields(self, data_file):
        f = h5py.File(data_file.filename, "r")
        fields = []
        cname = self.ds._particle_coordinates_name  # Coordinates
        mname = self.ds._particle_mass_name  # Coordinates

        for key in f.keys():
            if not key.startswith("PartType"):
                continue

            g = f[key]
            if cname not in g:
                continue

            ptype = str(key)
            for k in g.keys():
                    kk = k
                    # a little hack to deal with swift masses
                    if str(kk) == mname:
                        fields.append((ptype, "Mass"))
                        continue
                    if not hasattr(g[kk], "shape"):
                        continue
                    if len(g[kk].shape) > 1:
                        self._vector_fields[kk] = g[kk].shape[1]
                    fields.append((ptype, str(kk)))

        f.close()
        return fields, {}
