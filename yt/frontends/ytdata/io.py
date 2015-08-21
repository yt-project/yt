"""
YTData data-file handling function




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import h5py
import numpy as np

from yt.utilities.exceptions import *
from yt.funcs import \
    mylog

from yt.utilities.io_handler import \
    BaseIOHandler

from yt.utilities.lib.geometry_utils import \
    compute_morton

from yt.geometry.oct_container import \
    _ORDER_MAX

class IOHandlerYTDataHDF5(BaseIOHandler):
    _dataset_type = "ytdata_hdf5"

    def _read_fluid_selection(self, chunks, selector, fields, size):
        raise NotImplementedError

    def _read_particle_coords(self, chunks, ptf):
        # This will read chunks and yield the results.
        chunks = list(chunks)
        data_files = set([])
        # Only support halo reading for now.
        assert(len(ptf) == 1)
        assert(list(ptf.keys())[0] == "grid")
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        for data_file in sorted(data_files):
            all_count = self._count_particles(data_file)
            pcount = all_count["grid"]
            with h5py.File(data_file.filename, "r") as f:
                x = f["grid"]['x'].value.astype("float64")
                y = f["grid"]['y'].value.astype("float64")
                z = f["grid"]['z'].value.astype("float64")
                yield "grid", (x, y, z)

    def _read_particle_fields(self, chunks, ptf, selector):
        # Now we have all the sizes, and we can allocate
        chunks = list(chunks)
        data_files = set([])
        # Only support halo reading for now.
        assert(len(ptf) == 1)
        assert(list(ptf.keys())[0] == "grid")
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        for data_file in sorted(data_files):
            all_count = self._count_particles(data_file)
            pcount = all_count["grid"]
            with h5py.File(data_file.filename, "r") as f:
                for ptype, field_list in sorted(ptf.items()):
                    x = f["grid"]['x'].value.astype("float64")
                    y = f["grid"]['y'].value.astype("float64")
                    z = f["grid"]['z'].value.astype("float64")
                    mask = selector.select_points(x, y, z, 0.0)
                    del x, y, z
                    if mask is None: continue
                    for field in field_list:
                        data = f["grid"][field][mask].astype("float64")
                        yield (ptype, field), data

    def _initialize_index(self, data_file, regions):
        all_count = self._count_particles(data_file)
        pcount = all_count["grid"]
        morton = np.empty(pcount, dtype='uint64')
        mylog.debug("Initializing index % 5i (% 7i particles)",
                    data_file.file_id, pcount)
        ind = 0
        with h5py.File(data_file.filename, "r") as f:
            if not f["grid"].keys(): return None
            pos = np.empty((pcount, 3), dtype="float64")
            pos = data_file.ds.arr(pos, "code_length")
            dx = np.finfo(f["grid"]['x'].dtype).eps
            dx = 2.0*self.ds.quan(dx, "code_length")
            pos[:,0] = f["grid"]["x"].value
            pos[:,1] = f["grid"]["y"].value
            pos[:,2] = f["grid"]["z"].value
            # These are 32 bit numbers, so we give a little lee-way.
            # Otherwise, for big sets of particles, we often will bump into the
            # domain edges.  This helps alleviate that.
            np.clip(pos, self.ds.domain_left_edge + dx,
                         self.ds.domain_right_edge - dx, pos)
            if np.any(pos.min(axis=0) < self.ds.domain_left_edge) or \
               np.any(pos.max(axis=0) > self.ds.domain_right_edge):
                raise YTDomainOverflow(pos.min(axis=0),
                                       pos.max(axis=0),
                                       self.ds.domain_left_edge,
                                       self.ds.domain_right_edge)
            regions.add_data_file(pos, data_file.file_id)
            morton[ind:ind+pos.shape[0]] = compute_morton(
                pos[:,0], pos[:,1], pos[:,2],
                data_file.ds.domain_left_edge,
                data_file.ds.domain_right_edge)
        return morton

    def _count_particles(self, data_file):
        with h5py.File(data_file.filename, "r") as f:
            return {"grid": f["grid"].attrs["num_elements"]}

    def _identify_fields(self, data_file):
        with h5py.File(data_file.filename, "r") as f:
            fields = [("grid", str(field)) for field in f["grid"]]
            units = dict([(("grid", str(field)), 
                           f["grid"][field].attrs["units"]) for field in f["grid"]])
        return fields, units
