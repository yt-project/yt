"""
ytree data-file handling function




"""

#-----------------------------------------------------------------------------
# Copyright (c) yt Development Team. All rights reserved.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.utilities.on_demand_imports import _h5py as h5py
import numpy as np

from yt.utilities.exceptions import YTDomainOverflow
from yt.funcs import \
    mylog

from yt.utilities.io_handler import \
    BaseIOHandler

from yt.utilities.lib.geometry_utils import \
    compute_morton


class IOHandlerYTreeHDF5(BaseIOHandler):
    _dataset_type = "ytree_arbor"

    def _read_fluid_selection(self, chunks, selector, fields, size):
        raise NotImplementedError

    def _read_particle_coords(self, chunks, ptf):
        # This will read chunks and yield the results.
        chunks = list(chunks)
        data_files = set([])
        # Only support halo reading for now.
        assert(len(ptf) == 1)
        assert(list(ptf.keys())[0] == "halos")
        ptype = 'halos'
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        for data_file in sorted(data_files):
            with h5py.File(data_file.filename, "r") as f:
                units = self.ds._field_dict['position_x']['units']
                pos = data_file._get_particle_positions(ptype, f=f)
                x, y, z = (self.ds.arr(pos[:, i], units)
                           for i in range(3))
                yield "halos", (x, y, z)

    def _read_particle_fields(self, chunks, ptf, selector):
        # Now we have all the sizes, and we can allocate
        chunks = list(chunks)
        data_files = set([])
        # Only support halo reading for now.
        assert(len(ptf) == 1)
        assert(list(ptf.keys())[0] == "halos")
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        for data_file in sorted(data_files):
            with h5py.File(data_file.filename, "r") as f:
                for ptype, field_list in sorted(ptf.items()):
                    units = self.ds._field_dict['position_x']['units']
                    pos = data_file._get_particle_positions(ptype, f=f)
                    x, y, z = (self.ds.arr(pos[:, i], units)
                               for i in range(3))
                    mask = selector.select_points(x, y, z, 0.0)
                    del x, y, z
                    if mask is None: continue
                    for field in field_list:
                        data = f['data'][field][mask].astype("float64")
                        yield (ptype, field), data

    def _initialize_index(self, data_file, regions):
        pcount = data_file.particle_count
        morton = np.empty(pcount, dtype='uint64')
        mylog.debug("Initializing index % 5i (% 7i particles)",
                    data_file.file_id, pcount)
        ind = 0
        if pcount == 0: return None
        ptype = 'halos'
        with h5py.File(data_file.filename, "r") as f:
            units = self.ds._field_dict['position_x']['units']
            pos = data_file._get_particle_positions(ptype, f=f)
            pos = data_file.ds.arr(pos, units). to("code_length")
            dle = self.ds.domain_left_edge.to("code_length")
            dre = self.ds.domain_right_edge.to("code_length")
            if np.any(pos.min(axis=0) < dle) or \
               np.any(pos.max(axis=0) > dre):
                raise YTDomainOverflow(pos.min(axis=0),
                                       pos.max(axis=0),
                                       dle, dre)
            regions.add_data_file(pos, data_file.file_id)
            morton[ind:ind+pos.shape[0]] = compute_morton(
                pos[:,0], pos[:,1], pos[:,2], dle, dre)
        return morton

    def _count_particles(self, data_file):
        return {'halos': data_file.particle_count}

    def _identify_fields(self, data_file):
        fields = [('halos', field) for field in self.ds._field_dict]
        units = dict((('halos', field), self.ds._field_dict[field]['units'])
                     for field in self.ds._field_dict)
        return fields, units
