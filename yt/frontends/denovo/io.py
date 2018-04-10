"""
Denovo-specific IO functions



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.utilities.io_handler import \
    BaseIOHandler
from yt.utilities.on_demand_imports import _h5py as h5py
import numpy as np
from yt.utilities.logger import ytLogger as mylog

def field_dname(field_name):
    return "/denovo/{0}".format(field_name)

class IOHandlerDenovoHDF5(BaseIOHandler):
    _particle_reader = False
    _dataset_type = 'denovo'

    def __init__(self,ds):
        super(IOHandlerDenovoHDF5, self).__init__(ds)
        self._handle = ds._handle
        try:
            self.groups = ds.parameters['mesh_g']
        except KeyError:
            self.groups = False

    def _read_field_names(self, grid):
        # This function is used to pull in all of the fields from the Denovo
        # file relevant to what we can plot here.
        pass



    def _read_particle_coords(self, chunks, ptf):
        # At this time Denovo has no particles, so this function will not
        # return particle coords.
        pass

    def _read_particle_fields(self, chunks, ptf, selector):
        # At this time Denovo has no particle fields, so this function
        # will not return any particle fields.
        pass

    def _read_fluid_selection(self, chunks, selector, fields, size):
        # This needs to allocate a set of arrays inside a dictionary, where the
        # keys are the (ftype, fname) tuples and the values are arrays that
        # have been masked using whatever selector method is appropriate.  The
        # dict gets returned at the end and it should be flat, with selected
        # data.  Note that if you're reading grid data, you might need to
        # special-case a grid selector object.
        # Also note that "chunks" is a generator for multiple chunks, each of
        # which contains a list of grids. The returned numpy arrays should be
        # in 64-bit float and contiguous along the z direction. Therefore, for
        # a C-like input array with the dimension [x][y][z] or a
        # Fortran-like input array with the dimension (z,y,x), a matrix
        # transpose is required (e.g., using np_array.transpose() or
        # np_array.swapaxes(0,2)).

        chunks = list(chunks)
        assert(len(chunks) == 1)

        fhandle = self._handle
        rv = {}
        for field in fields:
            ftype, fname = field
            rv[field] = np.empty(size, dtype=fhandle[field_dname(fname)].dtype)
            ngrids = sum(len(chunk.objs) for chunk in chunks)
            mylog.debug("Reading %s cells of %s fields in %s blocks",
                    size, [fname for ft, fn in fields], ngrids)
        for field in fields:
            ftype, fname = field
            # right now only pull in 0th group of the flux to see if things
            # work.
            ds = np.array(fhandle[field_dname(fname)][0,...],
                    dtype="float64").transpose().copy()
            ind = 0
            for chunk in chunks:
                for g in chunk.objs:
                    ind += g.select(selector, ds, rv[field], ind) # caches
        return rv

    def _read_chunk_data(self, chunk, fields):
        # This reads the data from a single chunk without doing any selection,
        # and is only used for caching data that might be used by multiple
        # different selectors later. For instance, this can speed up ghost zone
        # computation.
        pass
