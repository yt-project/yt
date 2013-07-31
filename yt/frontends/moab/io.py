"""The data-file handling functions

Author: Anthony Scopatz <scopatz@gmail.com>
Affiliation: The University of Wisconsin-Madison

"""
import numpy as np
from yt.funcs import mylog
from yt.utilities.io_handler import BaseIOHandler


def field_dname(field_name):
    return "/tstt/elements/Hex8/tags/{0}".format(field_name)


# TODO all particle bits were removed
class IOHandlerMoabH5MHex8(BaseIOHandler):
    _data_style = "h5m"
    _offset_string = 'data:offsets=0'
    _data_string = 'data:datatype=0'

    def __init__(self, pf, *args, **kwargs):
        # TODO check if _num_per_stride is needed
        self._num_per_stride = kwargs.pop("num_per_stride", 1000000)
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
