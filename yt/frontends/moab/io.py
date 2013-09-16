"""MOAB-specific fields


"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
from yt.funcs import mylog
from yt.utilities.io_handler import BaseIOHandler


def field_dname(field_name):
    return "/tstt/elements/Hex8/tags/{0}".format(field_name)


# TODO all particle bits were removed
class IOHandlerMoabH5MHex8(BaseIOHandler):
    _data_style = "moab_hex8"

    def __init__(self, pf, *args, **kwargs):
        # TODO check if _num_per_stride is needed
        BaseIOHandler.__init__(self, *args, **kwargs)
        self.pf = pf
        self._handle = pf._handle

    def _read_fluid_selection(self, chunks, selector, fields, size):
        chunks = list(chunks)
        fhandle = self._handle
        rv = {}
        for field in fields:
            ftype, fname = field
            rv[field] = np.empty(size, dtype=fhandle[field_dname(fname)].dtype)
        ngrids = sum(len(chunk.objs) for chunk in chunks)
        mylog.debug("Reading %s cells of %s fields in %s blocks",
                    size, [fname for ftype, fname in fields], ngrids)
        x, y, z = fhandle['/tstt/nodes/coordinates'][:].T
        mask = selector.select_points(x, y, z)
        for field in fields:
            ftype, fname = field
            data = fhandle[field_dname(fname)][mask]
            rv[field][:] = data
        return rv
