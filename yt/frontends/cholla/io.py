import numpy as np

from yt.utilities.io_handler import BaseIOHandler
from yt.utilities.on_demand_imports import _h5py as h5py


class ChollaIOHandler(BaseIOHandler):
    _particle_reader = False
    _dataset_type = "cholla"

    def _read_particle_coords(self, chunks, ptf):
        raise NotImplementedError

    def _read_particle_fields(self, chunks, ptf, selector):
        raise NotImplementedError

    def _read_fluid_selection(self, chunks, selector, fields, size):
        data = {}
        for field in fields:
            data[field] = np.empty(size, dtype="float64")

        with h5py.File(self.ds.parameter_filename, "r") as fh:
            ind = 0
            for chunk in chunks:
                for grid in chunk.objs:
                    nd = 0
                    for field in fields:
                        ftype, fname = field
                        values = fh[fname][:].astype("=f8")
                        nd = grid.select(selector, values, data[field], ind)
                    ind += nd
        return data

    def _read_chunk_data(self, chunk, fields):
        raise NotImplementedError
