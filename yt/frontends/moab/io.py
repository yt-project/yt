import numpy as np

from yt.funcs import mylog
from yt.utilities.io_handler import BaseIOHandler


def field_dname(field_name):
    return f"/tstt/elements/Hex8/tags/{field_name}"


# TODO all particle bits were removed
class IOHandlerMoabH5MHex8(BaseIOHandler):
    _dataset_type = "moab_hex8"

    def __init__(self, ds):
        super().__init__(ds)
        self._handle = ds._handle

    def _read_fluid_selection(self, chunks, selector, fields, size):
        chunks = list(chunks)
        assert len(chunks) == 1
        fhandle = self._handle
        rv = {}
        for field in fields:
            ftype, fname = field
            rv[field] = np.empty(size, dtype=fhandle[field_dname(fname)].dtype)
        ngrids = sum(len(chunk.objs) for chunk in chunks)
        mylog.debug(
            "Reading %s cells of %s fields in %s blocks",
            size,
            [fname for ft, fn in fields],
            ngrids,
        )
        for field in fields:
            ftype, fname = field
            ds = np.array(fhandle[field_dname(fname)][:], dtype="float64")
            ind = 0
            for chunk in chunks:
                for g in chunk.objs:
                    ind += g.select(selector, ds, rv[field], ind)  # caches
        return rv


class IOHandlerMoabPyneHex8(BaseIOHandler):
    _dataset_type = "moab_hex8_pyne"

    def _read_fluid_selection(self, chunks, selector, fields, size):
        chunks = list(chunks)
        assert len(chunks) == 1
        rv = {}
        pyne_mesh = self.ds.pyne_mesh
        for field in fields:
            rv[field] = np.empty(size, dtype="float64")
        ngrids = sum(len(chunk.objs) for chunk in chunks)
        mylog.debug(
            "Reading %s cells of %s fields in %s blocks",
            size,
            [fname for ftype, fname in fields],
            ngrids,
        )
        for field in fields:
            ftype, fname = field
            if pyne_mesh.structured:
                tag = pyne_mesh.mesh.tag_get_handle("idx")
                hex_list = [ent for ent in pyne_mesh.structured_iterate_hex()]
                indices = pyne_mesh.mesh.tag_get_data(tag, hex_list).flatten()
            else:
                indices = slice(None)
            ds = np.asarray(getattr(pyne_mesh, fname)[indices], "float64")

            ind = 0
            for chunk in chunks:
                for g in chunk.objs:
                    ind += g.select(selector, ds, rv[field], ind)  # caches
        return rv
