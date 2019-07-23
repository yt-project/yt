"""
Rockstar data-file handling function




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
import os

from yt.funcs import mylog

from yt.utilities.exceptions import \
    YTDomainOverflow

from yt.utilities.io_handler import \
    BaseIOHandler

from .definitions import halo_dts
from yt.utilities.lib.geometry_utils import compute_morton

from operator import attrgetter

class IOHandlerRockstarBinary(BaseIOHandler):
    _dataset_type = "rockstar_binary"

    def __init__(self, *args, **kwargs):
        super(IOHandlerRockstarBinary, self).__init__(*args, **kwargs)
        self._halo_dt = halo_dts[self.ds.parameters['format_revision']]

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
        for data_file in sorted(data_files,key=attrgetter("filename")):
            pcount = data_file.header['num_halos']
            if pcount == 0:
                continue
            with open(data_file.filename, "rb") as f:
                pos = data_file._get_particle_positions(ptype, f=f)
                yield "halos", (pos[:, i] for i in range(3))

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
        for data_file in sorted(data_files,key=attrgetter("filename")):
            pcount = data_file.header['num_halos']
            if pcount == 0:
                continue
            with open(data_file.filename, "rb") as f:
                for ptype, field_list in sorted(ptf.items()):
                    pos = data_file._get_particle_positions(ptype, f=f)
                    x, y, z = (pos[:, i] for i in range(3))
                    mask = selector.select_points(x, y, z, 0.0)
                    del x, y, z
                    f.seek(data_file._position_offset, os.SEEK_SET)
                    halos = np.fromfile(f, dtype=self._halo_dt, count = pcount)
                    if mask is None: continue
                    for field in field_list:
                        data = halos[field][mask].astype("float64")
                        yield (ptype, field), data

    def _initialize_index(self, data_file, regions):
        pcount = data_file.header["num_halos"]
        morton = np.empty(pcount, dtype='uint64')
        mylog.debug("Initializing index % 5i (% 7i particles)",
                    data_file.file_id, pcount)
        if pcount == 0:
            return morton
        ind = 0
        ptype = 'halos'
        with open(data_file.filename, "rb") as f:
            pos = data_file._get_particle_positions(ptype, f=f)
            pos = data_file.ds.arr(pos, "code_length")
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
        return {'halos': data_file.header['num_halos']}

    def _identify_fields(self, data_file):
        fields = [("halos", f) for f in self._halo_dt.fields if
                  "padding" not in f]
        return fields, {}
