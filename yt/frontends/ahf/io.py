"""
AHF-specific IO functions



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2017, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from operator import attrgetter

import numpy as np

from yt.funcs import \
    mylog
from yt.utilities.exceptions import \
    YTDomainOverflow
from yt.utilities.io_handler import \
    BaseIOHandler
from yt.utilities.lib.geometry_utils import \
    compute_morton


class IOHandlerAHFHalos(BaseIOHandler):
    _particle_reader = False
    _dataset_type = 'ahf'

    def _read_fluid_selection(self, chunks, selector, fields, size):
        raise NotImplementedError

    def _read_particle_coords(self, chunks, ptf):
        # This needs to *yield* a series of tuples of (ptype, (x, y, z)).
        # chunks is a list of chunks, and ptf is a dict where the keys are
        # ptypes and the values are lists of fields.
        for data_file in self._get_data_files(chunks, ptf):
            pos = data_file._get_particle_positions('halos')
            x, y, z = (pos[:, i] for i in range(3))
            yield 'halos', (x, y, z)

    def _read_particle_fields(self, chunks, ptf, selector):
        # This gets called after the arrays have been allocated.  It needs to
        # yield ((ptype, field), data) where data is the masked results of
        # reading ptype, field and applying the selector to the data read in.
        # Selector objects have a .select_points(x,y,z) that returns a mask, so
        # you need to do your masking here.
        for data_file in self._get_data_files(chunks, ptf):
            cols = []
            for field_list in ptf.values():
                cols.extend(field_list)
            cols = list(set(cols))
            halos = data_file.read_data(usecols=cols)
            pos = data_file._get_particle_positions('halos')
            x, y, z = (pos[:, i] for i in range(3))
            yield 'halos', (x, y, z)
            mask = selector.select_points(x, y, z, 0.0)
            del x, y, z
            if mask is None: continue
            for ptype, field_list in sorted(ptf.items()):
                for field in field_list:
                    data = halos[field][mask].astype('float64')
                    yield (ptype, field), data

    def _initialize_index(self, data_file, regions):
        halos = data_file.read_data(usecols=['ID'])
        pcount = len(halos['ID'])
        morton = np.empty(pcount, dtype='uint64')
        mylog.debug('Initializing index % 5i (% 7i particles)',
                    data_file.file_id, pcount)
        if pcount == 0:
            return morton
        ind = 0
        pos = data_file._get_particle_positions('halos')
        pos = data_file.ds.arr(pos, 'code_length')
        dle = self.ds.domain_left_edge
        dre = self.ds.domain_right_edge
        if np.any(pos.min(axis=0) < dle) or np.any(pos.max(axis=0) > dre):
            raise YTDomainOverflow(pos.min(axis=0),
                                   pos.max(axis=0),
                                   dle, dre)
        regions.add_data_file(pos, data_file.file_id)
        morton[ind:ind+pos.shape[0]] = compute_morton(
            pos[:, 0], pos[:, 1], pos[:, 2], dle, dre)
        return morton

    def _count_particles(self, data_file):
        halos = data_file.read_data(usecols=['ID'])
        return {'halos': len(halos['ID'])}

    def _identify_fields(self, data_file):
        fields = [('halos', f) for f in data_file.col_names]
        return fields, {}

    # Helper methods

    def _get_data_files(self, chunks, ptf):
        # Only support halo reading for now.
        assert len(ptf) == 1
        assert list(ptf.keys())[0] == 'halos'
        # Get data_files
        chunks = list(chunks)
        data_files = set([])
        for chunk in chunks:
            for obj in chunk.objs:
                data_files.update(obj.data_files)
        data_files = sorted(data_files, key=attrgetter('filename'))
        for data_file in data_files:
            yield data_file
