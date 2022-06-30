"""
CF-Radial-specific IO functions



"""

import numpy as np

from yt.utilities.io_handler import BaseIOHandler


class CFRadialIOHandler(BaseIOHandler):
    _particle_reader = False
    _dataset_type = "cf_radial"

    def _read_fluid_selection(self, chunks, selector, fields, size):
        # This needs to allocate a set of arrays inside a dictionary, where the
        # keys are the (ftype, fname) tuples and the values are arrays that
        # have been masked using whatever selector method is appropriate.  The
        # dict gets returned at the end and it should be flat, with selected
        # data.

        rv = {field: np.empty(size, dtype="float64") for field in fields}

        offset = 0

        with self.ds._handle() as xr_ds_handle:
            for field in fields:
                for chunk in chunks:
                    for grid in chunk.objs:
                        variable = xr_ds_handle.variables[field[1]]
                        data = variable.values[0, ...].T
                        offset += grid.select(selector, data, rv[field], offset)
        return rv
