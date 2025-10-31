from yt.utilities.io_handler import BaseIOHandler

from .misc import _CachedH5Openner


class ChollaIOHandler(BaseIOHandler):
    _particle_reader = False
    _dataset_type = "cholla"

    def _read_particle_coords(self, chunks, ptf):
        raise NotImplementedError

    def _read_particle_fields(self, chunks, ptf, selector):
        raise NotImplementedError

    def io_iter(self, chunks, fields):
        # this is loosely inspired by the implementation used for Enzo/Enzo-E
        # - those other implementations use the lower-level hdf5 interface. Unclear
        #   whether that affords any advantages...
        mapper = self.ds.index._block_mapping
        with _CachedH5Openner(mode="r") as h5_context_manager:
            for chunk in chunks:
                for obj in chunk.objs:
                    if obj.filename is None:  # unclear when this case arises...
                        continue

                    # ensure the file containing data for obj is open
                    fh = h5_context_manager.open_fh(obj.filename)

                    # access the HDF5 group containing the datasets of field values
                    grp = fh[mapper.h5_group]
                    # get the indices in a generic dataset that correspond to obj.id
                    idx = mapper.idx_map[obj.id]

                    for field in fields:
                        ftype, fname = field
                        yield field, obj, grp[fname][idx].astype("=f8")

    def _read_chunk_data(self, chunk, fields):
        raise NotImplementedError
