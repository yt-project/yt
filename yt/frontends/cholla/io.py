from yt.utilities.io_handler import BaseIOHandler
from yt.utilities.on_demand_imports import _h5py as h5py


class ChollaIOHandler(BaseIOHandler):
    _particle_reader = False
    _dataset_type = "cholla"

    def _read_particle_coords(self, chunks, ptf):
        raise NotImplementedError

    def _read_particle_fields(self, chunks, ptf, selector):
        raise NotImplementedError

    def io_iter(self, chunks, fields):
        # this is loosely inspired by the implementation used for Enzo/Enzo-E
        # - those other options use the lower-level hdf5 interface. Unclear
        #   whether that affords any advantages...
        fh, filename = None, None
        for chunk in chunks:
            for obj in chunk.objs:
                if obj.filename is None:  # unclear when this case arises...
                    continue
                elif obj.filename != filename:
                    if fh is not None:
                        fh.close()
                    fh, filename = h5py.File(obj.filename, "r"), obj.filename
                for field in fields:
                    ftype, fname = field
                    yield field, obj, fh[fname][:].astype("=f8")
        if fh is not None:
            fh.close()

    def _read_chunk_data(self, chunk, fields):
        raise NotImplementedError
