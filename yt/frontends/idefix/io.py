import numpy as np

from yt.frontends.idefix.dmpfile_io import read_idefix_dmpfile
from yt.utilities.io_handler import BaseIOHandler


class IdefixIOHandler(BaseIOHandler):
    _particle_reader = False
    _dataset_type = "idefix"

    def __init__(self, ds):
        BaseIOHandler.__init__(self, ds)
        self.ds = ds
        self.dmpfile = ds.parameter_filename

    def _read_fluid_selection(self, chunks, selector, fields, size):
        # THIS IS STRAIGHT FROM netCDF FRONTEND
        data = {}
        offset = 0

        for field in fields:
            ftype, fname = field
            data[ftype, fname] = np.empty(size, dtype="float64")
            for chunk in chunks:
                for grid in chunk.objs:
                    _fprops, fdata = read_idefix_dmpfile(self.dmpfile)
                    values = fdata[fname]
                    offset += grid.select(selector, values, data[field], offset)
        return data

    def _read_particle_coords(self, chunks, ptf):
        # This needs to *yield* a series of tuples of (ptype, (x, y, z)).
        # chunks is a list of chunks, and ptf is a dict where the keys are
        # ptypes and the values are lists of fields.
        raise NotImplementedError

    def _read_particle_fields(self, chunks, ptf, selector):
        # idefix doesn't have particles (yet)
        raise NotImplementedError

    def _read_chunk_data(self, chunk, fields):
        # This reads the data from a single chunk without doing any selection,
        # and is only used for caching data that might be used by multiple
        # different selectors later. For instance, this can speed up ghost zone
        # computation.
        raise NotImplementedError
