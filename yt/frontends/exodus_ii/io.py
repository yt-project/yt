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
from .util import ExodusIIData
import numpy as np

from yt.utilities.io_handler import \
    BaseIOHandler

class IOHandlerExodusII(BaseIOHandler):
    _particle_reader = False
    _dataset_type = "exodus_ii"
    _node_types = ("diffused", "convected")
    _INDEX_OFFSET = 1

    def __init__(self):
        self.filename = filename
        self.ds       = ExodusIIData(filename)
        self.ds.read()

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
        chunks = list(chunks) # chunks in this case correspond to mesh_id or slices in the ExodusII data
        rv = {}
        for field in fields:
            rv[field] = self.ds.arr(np.empty(size, dtype="float64"))

        for field in fields:
            ftype, fname = field
            ind = _INDEX_OFFSET
            for chunk in chunks

        return rv

    def _read_chunk_data(self, chunk, fields):
        # This reads the data from a single chunk, and is only used for
        # caching.
        pass
