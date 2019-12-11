"""
AdaptaHOP data-file handling function




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

from yt.utilities.lib.geometry_utils import compute_morton
from yt.utilities.cython_fortran_utils import FortranFile

from operator import attrgetter

from .definitions import \
    HEADER_ATTRIBUTES, HALO_ATTRIBUTES

class IOHandlerAdaptaHOPBinary(BaseIOHandler):
    _dataset_type = "adaptahop_binary"

    def _read_fluid_selection(self, chunks, selector, fields, size):
        raise NotImplementedError

    def _read_chunk_data(self, chunk, fields):
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
        pcount = data_file.ds.parameters['nhalos'] + data_file.ds.parameters['nsubs']
        morton = np.empty(pcount, dtype='uint64')
        mylog.debug("Initializing index % 5i (% 7i particles)",
                    data_file.file_id, pcount)
        if pcount == 0:
            return morton
        ind = 0
        ptype = 'halos'

        pos = self._get_particle_positions()
        print(pos.min(axis=0), pos.max(axis=0))
        pos = data_file.ds.arr(pos, "cm") + self.ds.domain_width.to('Mpc') / 2  # FIXME
        pos = pos.to('code_length')
        print(pos.min(axis=0), pos.max(axis=0), self.ds.domain_width.to('Mpc'))
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
        nhalos = data_file.ds.parameters['nhalos'] + data_file.ds.parameters['nsubs']
        return {'halos': nhalos}

    def _identify_fields(self, data_file):
        fields = []
        for attr, _1, _2 in HALO_ATTRIBUTES:
            if isinstance(attr, str):
                fields.append(('halos', attr))
            else:
                for a in attr:
                    fields.append(('halos', a))
        return fields, {}


    # Specific to AdaptaHOP
    def _get_particle_positions(self):
        with FortranFile(self.ds.parameter_filename) as fpu:
            params = fpu.read_attrs(HEADER_ATTRIBUTES)

            todo = _todo_from_attributes(('x', 'y', 'z'))

            nhalos = params['nhalos'] + params['nsubs']
            data = np.zeros((nhalos, 3))
            for ihalo in range(nhalos):
                for it in todo:
                    if isinstance(it, int):
                        fpu.skip(it)
                    else:
                        dt = fpu.read_attrs(it)
                        data[ihalo, 0] = dt['x']
                        data[ihalo, 1] = dt['y']
                        data[ihalo, 2] = dt['z']

        return data * 3.08e24

            
def _todo_from_attributes(attributes):
    iskip = 0
    todo = []

    attributes = set(attributes)

    for i, (attrs, l, k) in enumerate(HALO_ATTRIBUTES):
        if not isinstance(attrs, tuple):
            attrs = (attrs, )
        ok = False
        for attr in attrs:
            if attr in attributes:
                ok = True
                break

        if i == 0:
            state = 'read' if ok else 'skip'

        if ok:
            if state == 'skip':
                # Switched from skip to read, store skip information and start
                # new read list
                todo.append(iskip)
                todo.append([])
                iskip = 0
            todo[-1].append((attrs, l, k))
            state = 'read'
        else:
            iskip += 1
            state = 'skip'

    if state == 'skip' and iskip > 0:
        todo.append(iskip)

    return todo
