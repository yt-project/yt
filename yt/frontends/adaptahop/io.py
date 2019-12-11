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
            pcount = data_file.ds.parameters['nhalos'] + data_file.ds.parameters['nsubs']
            if pcount == 0:
                continue
            pos = self._get_particle_positions()
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
        for data_file in sorted(data_files, key=attrgetter("filename")):
            pcount = data_file.ds.parameters['nhalos'] + data_file.ds.parameters['nsubs']
            if pcount == 0:
                continue
            ptype = 'halos'
            field_list0 = sorted(ptf[ptype], key=_find_attr_position)
            field_list_pos = ['particle_position_%s' % k for k in 'xyz']
            field_list = sorted(
                set(field_list0 + field_list_pos),
                key=_find_attr_position
            )

            with FortranFile(self.ds.parameter_filename) as fpu:
                params = fpu.read_attrs(HEADER_ATTRIBUTES)

                todo = _todo_from_attributes(
                    field_list
                )

                attr2pos = {f: i for i, f in enumerate(field_list)}

                nhalos = params['nhalos'] + params['nsubs']
                data = np.zeros((nhalos, len(field_list)))
                for ihalo in range(nhalos):
                    for it in todo:
                        if isinstance(it, int):
                            fpu.skip(it)
                        else:
                            for k, v in fpu.read_attrs(it).items():
                                data[ihalo, attr2pos[k]] = v
            ipos = [field_list.index(k) for k in field_list_pos]
            x, y, z = (data[:, i] for i in ipos)
            mask = selector.select_points(x, y, z, 0.0)
            del x, y, z

            if mask is None: continue
            for field in field_list0:
                i = field_list.index(field)
                yield (ptype, field), data[mask, i].astype('float64')

    def _initialize_index(self, data_file, regions):
        pcount = data_file.ds.parameters['nhalos'] + data_file.ds.parameters['nsubs']
        morton = np.empty(pcount, dtype='uint64')
        mylog.debug("Initializing index % 5i (% 7i particles)",
                    data_file.file_id, pcount)
        if pcount == 0:
            return morton
        ind = 0

        pos = self._get_particle_positions()
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

    #-----------------------------------------------------
    # Specific to AdaptaHOP
    def _get_particle_positions(self):
        """Read the particles and return them in code_units"""
        with FortranFile(self.ds.parameter_filename) as fpu:
            params = fpu.read_attrs(HEADER_ATTRIBUTES)

            todo = _todo_from_attributes(
                ('particle_position_x', 'particle_position_y', 'particle_position_z'))

            nhalos = params['nhalos'] + params['nsubs']
            data = np.zeros((nhalos, 3))
            for ihalo in range(nhalos):
                for it in todo:
                    if isinstance(it, int):
                        fpu.skip(it)
                    else:
                        # Small optimisation here: we can read as vector
                        # dt = fpu.read_attrs(it)
                        # data[ihalo, 0] = dt['particle_position_x']
                        # data[ihalo, 1] = dt['particle_position_y']
                        # data[ihalo, 2] = dt['particle_position_z']
                        data[ihalo, :] = fpu.read_vector(it[0][-1])

        data = self.ds.arr(data, "code_length") + self.ds.domain_width / 2  # FIXME
        return data


def _todo_from_attributes(attributes):
    # Helper function to generate a list of read-skip instructions given a list of
    # attributes. This is used to skip fields most of the fields when reading
    # the tree_brick files.
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

def _find_attr_position(key):
    j = 0
    for attrs, l, k in HALO_ATTRIBUTES:
        if not isinstance(attrs, tuple):
            attrs = (attrs, )
        for a in attrs:
            if key == a:
                return j
            j += 1
    raise KeyError