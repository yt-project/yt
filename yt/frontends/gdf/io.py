"""
The data-file handling functions



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
import h5py
from yt.funcs import \
    mylog
from yt.utilities.io_handler import \
    BaseIOHandler


def field_dname(grid_id, field_name):
    return "/data/grid_%010i/%s" % (grid_id, field_name)


# TODO all particle bits were removed
class IOHandlerGDFHDF5(BaseIOHandler):
    _data_style = "grid_data_format"
    _offset_string = 'data:offsets=0'
    _data_string = 'data:datatype=0'

    def __init__(self, pf, *args, **kwargs):
        super(IOHandlerGDFHDF5, self).__init__(pf)
        # Now we cache the particle fields

    def _read_fluid_selection(self, chunks, selector, fields, size):
        chunks = list(chunks)
        if any((ftype != "gas" for ftype, fname in fields)):
            raise NotImplementedError
        fhandle = h5py.File(self.pf.parameter_filename, "r")
        rv = {}
        for field in fields:
            ftype, fname = field
            rv[field] = np.empty(
                size, dtype=fhandle[field_dname(0, fname)].dtype)
        ngrids = sum(len(chunk.objs) for chunk in chunks)
        mylog.debug("Reading %s cells of %s fields in %s blocks",
                    size, [fname for ftype, fname in fields], ngrids)
        for field in fields:
            ftype, fname = field
            ind = 0
            for chunk in chunks:
                for grid in chunk.objs:
                    data = fhandle[field_dname(grid.id, fname)][:]
                    if self.pf.field_ordering == 1:
                        data = data.swapaxes(0, 2)
                    ind += grid.select(selector, data, rv[field], ind) # caches
        fhandle.close()
        return rv
