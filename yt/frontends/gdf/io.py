import numpy as np

from yt.funcs import mylog
from yt.geometry.selection_routines import GridSelector
from yt.utilities.io_handler import BaseParticleIOHandler
from yt.utilities.on_demand_imports import _h5py as h5py


def _grid_dname(grid_id):
    return "/data/grid_%010i" % grid_id


def _field_dname(grid_id, field_name):
    return f"{_grid_dname(grid_id)}/{field_name}"


# TODO all particle bits were removed
class IOHandlerGDFHDF5(BaseParticleIOHandler):
    _dataset_type = "grid_data_format"
    _offset_string = "data:offsets=0"
    _data_string = "data:datatype=0"

    def _read_fluid_selection(self, chunks, selector, fields, size):

        rv = {}
        chunks = list(chunks)

        if isinstance(selector, GridSelector):
            if not (len(chunks) == len(chunks[0].objs) == 1):
                raise RuntimeError
            grid = chunks[0].objs[0]
            h5f = h5py.File(grid.filename, mode="r")
            gds = h5f.get(_grid_dname(grid.id))
            for ftype, fname in fields:
                if self.ds.field_ordering == 1:
                    rv[(ftype, fname)] = gds.get(fname)[()].swapaxes(0, 2)
                else:
                    rv[(ftype, fname)] = gds.get(fname)[()]
            h5f.close()
            return rv
        if size is None:
            size = sum(grid.count(selector) for chunk in chunks for grid in chunk.objs)

        if any((ftype != "gdf" for ftype, fname in fields)):
            raise NotImplementedError

        for field in fields:
            ftype, fname = field
            fsize = size
            # check the dtype instead
            rv[field] = np.empty(fsize, dtype="float64")
        ngrids = sum(len(chunk.objs) for chunk in chunks)
        mylog.debug(
            "Reading %s cells of %s fields in %s blocks",
            size,
            [fn for ft, fn in fields],
            ngrids,
        )
        ind = 0
        for chunk in chunks:
            fid = None
            for grid in chunk.objs:
                if grid.filename is None:
                    continue
                if fid is None:
                    fid = h5py.h5f.open(
                        bytes(grid.filename, "utf-8"), h5py.h5f.ACC_RDONLY
                    )
                if self.ds.field_ordering == 1:
                    # check the dtype instead
                    data = np.empty(grid.ActiveDimensions[::-1], dtype="float64")
                    data_view = data.swapaxes(0, 2)
                else:
                    # check the dtype instead
                    data_view = data = np.empty(grid.ActiveDimensions, dtype="float64")
                for field in fields:
                    ftype, fname = field
                    dg = h5py.h5d.open(
                        fid, bytes(_field_dname(grid.id, fname), "utf-8")
                    )
                    dg.read(h5py.h5s.ALL, h5py.h5s.ALL, data)
                    # caches
                    nd = grid.select(selector, data_view, rv[field], ind)
                ind += nd  # I don't get that part, only last nd is added
            if fid is not None:
                fid.close()
        return rv
