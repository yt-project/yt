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
        # will not return any particle fields. This may be updated later if
        # support for output data from shift is added.
        pass

    def _read_fluid_selection(self, chunks, selector, fields, size):
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
            if 'egroup' in ftype:
                # trying to get the energy group number from the egroup string
                # isn't that pretty, but it works for now.
                group_no = int(ftype.split('_')[-1])
                mylog.debug("reading in {} data for energy group {}".format(fname,
                        group_no))
                ds = np.array(fhandle[field_dname(fname)][group_no,...],
                        dtype="float64").transpose().copy()
                ind = 0
                for chunk in chunks:
                    for g in chunk.objs:
                        ind += g.select(selector, ds, rv[field], ind) # caches
            else:
                ds = np.array(fhandle[field_dname(fname)][0,...],
                        dtype="float64").transpose().copy()
                ind = 0
                for chunk in chunks:
                    for g in chunk.objs:
                        ind += g.select(selector, ds, rv[field], ind) # caches

        return rv

    def _read_chunk_data(self, chunk, fields):
        pass
