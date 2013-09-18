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
        assert size
        chunks = list(chunks)
        assert(len(chunks) == 1)
        fhandle = self._handle
        rv = {}
        for field in fields:
            ftype, fname = field
            rv[field] = np.empty(size, dtype=fhandle[field_dname(fname)].dtype)
        ngrids = sum(len(chunk.objs) for chunk in chunks)
        mylog.debug("Reading %s cells of %s fields in %s blocks",
                    size, [fname for ftype, fname in fields], ngrids)
        for field in fields:
            ftype, fname = field
            ds = np.array(fhandle[field_dname(fname)][:], dtype="float64")
            ind = 0
            for chunk in chunks:
                for g in chunk.objs:
                    ind += g.select(selector, ds, rv[field], ind) # caches
        return rv

class IOHandlerMoabPyneHex8(BaseIOHandler):
    _data_style = "moab_hex8_pyne"

    def __init__(self, pf, *args, **kwargs):
        # TODO check if _num_per_stride is needed
        BaseIOHandler.__init__(self, *args, **kwargs)
        self.pf = pf

    def _read_fluid_selection(self, chunks, selector, fields, size):
        assert size
        chunks = list(chunks)
        assert(len(chunks) == 1)
        fhandle = self._handle
        rv = {}
        for field in fields:
            ftype, fname = field
            rv[field] = np.empty(size, dtype=fhandle[field_dname(fname)].dtype)
        ngrids = sum(len(chunk.objs) for chunk in chunks)
        mylog.debug("Reading %s cells of %s fields in %s blocks",
                    size, [fname for ftype, fname in fields], ngrids)
        for field in fields:
            ftype, fname = field
            ds = np.array(fhandle[field_dname(fname)][:], dtype="float64")
            ind = 0
            for chunk in chunks:
                for g in chunk.objs:
                    ind += g.select(selector, ds, rv[field], ind) # caches
        return rv
