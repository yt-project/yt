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
                        data = f[ptype][field][mask].astype("float64")
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
            return dict([(ptype, f[ptype].attrs["Number_of_groups"]) \
                         for ptype in "FOF", "SUBFIND"])

    def _identify_fields(self, data_file):
        with h5py.File(data_file.filename, "r") as f:
            fields = []
            for ptype in ["FOF", "SUBFIND"]:
                for ax in "xyz":
                    fields.append((ptype, "particle_position_%s" % ax))
        return fields, {}
