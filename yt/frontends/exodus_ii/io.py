"""
ExodusII-specific IO functions



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
from yt.utilities.io_handler import \
    BaseIOHandler
from yt.utilities.file_handler import \
    NetCDF4FileHandler


class IOHandlerExodusII(BaseIOHandler):
    _particle_reader = False
    _dataset_type = "exodus_ii"
    _INDEX_OFFSET = 1

    def __init__(self, ds):
        self.filename = ds.index_filename
        exodus_ii_handler = NetCDF4FileHandler(self.filename)
        self.handler = exodus_ii_handler.dataset
        super(IOHandlerExodusII, self).__init__(ds)
        self.node_fields = ds._get_nod_names()
        self.elem_fields = ds._get_elem_names()

    def _read_particle_coords(self, chunks, ptf):
        pass

    def _read_particle_fields(self, chunks, ptf, selector):
        pass

    def _read_fluid_selection(self, chunks, selector, fields, size):
        # This needs to allocate a set of arrays inside a dictionary, where the
        # keys are the (ftype, fname) tuples and the values are arrays that
        # have been masked using whatever selector method is appropriate.  The
        # dict gets returned at the end and it should be flat, with selected
        # data.  Note that if you're reading grid data, you might need to
        # special-case a grid selector object.
        chunks = list(chunks)
        rv = {}
        for field in fields:
            ftype, fname = field
            ci = self.handler.variables[ftype][:] - self._INDEX_OFFSET
            num_elem = ci.shape[0]
            if fname in self.node_fields:
                nodes_per_element = ci.shape[1]
                rv[field] = np.empty((num_elem, nodes_per_element), dtype="float64")
            elif fname in self.elem_fields:
                rv[field] = np.empty(num_elem, dtype="float64")
        for field in fields:
            ind = 0
            ftype, fname = field
            mesh_id = int(ftype[-1])
            chunk = chunks[mesh_id - 1]
            ci = self.handler.variables[ftype][:] - self._INDEX_OFFSET
            if fname in self.node_fields:
                field_ind = self.node_fields.index(fname)
                fdata = self.handler.variables['vals_nod_var%d' % (field_ind + 1)]
                data = fdata[self.ds.step][ci]
                for g in chunk.objs:
                    ind += g.select(selector, data, rv[field], ind)  # caches
            if fname in self.elem_fields:
                field_ind = self.elem_fields.index(fname)
                fdata = self.handler.variables['vals_elem_var%deb%s' %
                                               (field_ind + 1, mesh_id)][:]
                data = fdata[self.ds.step, :]
                for g in chunk.objs:
                    ind += g.select(selector, data, rv[field], ind)  # caches
        return rv

    def _read_chunk_data(self, chunk, fields):
        # This reads the data from a single chunk, and is only used for
        # caching.
        pass
