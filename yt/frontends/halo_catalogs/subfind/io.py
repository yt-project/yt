"""
Subfind data-file handling function




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import h5py
import numpy as np

from yt.utilities.exceptions import *
from yt.funcs import mylog

from yt.utilities.io_handler import \
    BaseIOHandler

from yt.utilities.lib.geometry_utils import compute_morton

from yt.geometry.oct_container import _ORDER_MAX

class IOHandlerSubfindHDF5(BaseIOHandler):
    _dataset_type = "subfind_hdf5"

    def _read_fluid_selection(self, chunks, selector, fields, size):
        raise NotImplementedError

    def _read_particle_coords(self, chunks, ptf):
        # This will read chunks and yield the results.
        chunks = list(chunks)
        data_files = set([])
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        for data_file in data_files:
            with h5py.File(data_file.filename, "r") as f:
                for ptype, field_list in sorted(ptf.items()):
                    pcount = f[ptype].attrs["Number_of_groups"]
                    coords = f[ptype]["CenterOfMass"].value.astype("float64")
                    coords = np.resize(coords, (pcount, 3))
                    x = coords[:, 0]
                    y = coords[:, 1]
                    z = coords[:, 2]
                    yield ptype, (x, y, z)

    def _read_particle_fields(self, chunks, ptf, selector):
        # Now we have all the sizes, and we can allocate
        chunks = list(chunks)
        data_files = set([])
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        for data_file in data_files:
            with h5py.File(data_file.filename, "r") as f:
                for ptype, field_list in sorted(ptf.items()):
                    pcount = f[ptype].attrs["Number_of_groups"]
                    coords = f[ptype]["CenterOfMass"].value.astype("float64")
                    coords = np.resize(coords, (pcount, 3))
                    x = coords[:, 0]
                    y = coords[:, 1]
                    z = coords[:, 2]
                    mask = selector.select_points(x, y, z)
                    del x, y, z
                    if mask is None: continue
                    for field in field_list:
                        if field in f[ptype].keys():
                            field_data = f[ptype][field].value.astype("float64")
                        else:
                            fname = field[:field.rfind("_")]
                            field_data = f[ptype][fname].value.astype("float64")
                            my_div = field_data.size / pcount
                            if my_div > 1:
                                field_data = np.resize(field_data, (pcount, my_div))
                                findex = int(field[field.rfind("_") + 1:])
                                field_data = field_data[:, findex]
                        data = field_data[mask]
                        yield (ptype, field), data

    def _initialize_index(self, data_file, regions):
        pcount = sum(self._count_particles(data_file).values())
        morton = np.empty(pcount, dtype='uint64')
        mylog.debug("Initializing index % 5i (% 7i particles)",
                    data_file.file_id, pcount)
        ind = 0
        with h5py.File(data_file.filename, "r") as f:
            if not f.keys(): return None
            dx = np.finfo(f["FOF"]['CenterOfMass'].dtype).eps
            dx = 2.0*self.pf.quan(dx, "code_length")
            for ptype in ["FOF", "SUBFIND"]:
                my_pcount = f[ptype].attrs["Number_of_groups"]
                pos = f[ptype]["CenterOfMass"].value.astype("float64")
                pos = np.resize(pos, (my_pcount, 3))
                pos = data_file.pf.arr(pos, "code_length")
                
                # These are 32 bit numbers, so we give a little lee-way.
                # Otherwise, for big sets of particles, we often will bump into the
                # domain edges.  This helps alleviate that.
                np.clip(pos, self.pf.domain_left_edge + dx,
                             self.pf.domain_right_edge - dx, pos)
                if np.any(pos.min(axis=0) < self.pf.domain_left_edge) or \
                   np.any(pos.max(axis=0) > self.pf.domain_right_edge):
                    raise YTDomainOverflow(pos.min(axis=0),
                                           pos.max(axis=0),
                                           self.pf.domain_left_edge,
                                           self.pf.domain_right_edge)
                regions.add_data_file(pos, data_file.file_id)
                morton[ind:ind+pos.shape[0]] = compute_morton(
                    pos[:,0], pos[:,1], pos[:,2],
                    data_file.pf.domain_left_edge,
                    data_file.pf.domain_right_edge)
                ind += pos.shape[0]
        return morton

    def _count_particles(self, data_file):
        with h5py.File(data_file.filename, "r") as f:
            data_file.particle_count = \
              dict([(ptype, f[ptype].attrs["Number_of_groups"])
                    for ptype in self.pf.particle_types_raw])
            return data_file.particle_count

    def _identify_fields(self, data_file):
        pcount = self._count_particles(data_file)
        fields = []
        with h5py.File(data_file.filename, "r") as f:
            for ptype in self.pf.particle_types_raw:
                fields.extend(h5_field_list(f[ptype], ptype, pcount[ptype]))
        return fields, {}

def h5_field_list(fh, ptype, pcount):
    fields = []
    for field in fh.keys():
        if isinstance(fh[field], h5py.Group):
            fields.extend(h5_field_list(fh[field], ptype, pcount))
        else:
            if fh[field].size % pcount:
                pass
                # mylog.warn("Cannot add field (%s, %s) with size %d for %d particles." % \
                #            (ptype, field, fh[field].size, pcount))
            else:
                my_div = fh[field].size / pcount
                fname = fh[field].name[fh[field].name.find(ptype) + len(ptype) + 1:]
                if my_div > 1:
                    for i in range(my_div):
                        fields.append((ptype, "%s_%d" % (fname, i)))
                else:
                    fields.append((ptype, fname))
    return fields
