import numpy as np

from yt.utilities.file_handler import NetCDF4FileHandler
from yt.utilities.io_handler import BaseIOHandler


class IOHandlerExodusII(BaseIOHandler):
    _particle_reader = False
    _dataset_type = "exodus_ii"
    _INDEX_OFFSET = 1

    def __init__(self, ds):
        self.filename = ds.index_filename
        exodus_ii_handler = NetCDF4FileHandler(self.filename)
        self.handler = exodus_ii_handler
        super().__init__(ds)
        self.node_fields = ds._get_nod_names()
        self.elem_fields = ds._get_elem_names()

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
        with self.handler.open_ds() as ds:
            chunks = list(chunks)
            rv = {}
            for field in fields:
                ftype, fname = field
                if ftype == "all":
                    ci = np.concatenate(
                        [
                            mesh.connectivity_indices - self._INDEX_OFFSET
                            for mesh in self.ds.index.mesh_union
                        ]
                    )
                else:
                    ci = ds.variables[ftype][:] - self._INDEX_OFFSET
                num_elem = ci.shape[0]
                if fname in self.node_fields:
                    nodes_per_element = ci.shape[1]
                    rv[field] = np.zeros((num_elem, nodes_per_element), dtype="float64")
                elif fname in self.elem_fields:
                    rv[field] = np.zeros(num_elem, dtype="float64")
            for field in fields:
                ind = 0
                ftype, fname = field
                if ftype == "all":
                    mesh_ids = [mesh.mesh_id + 1 for mesh in self.ds.index.mesh_union]
                    objs = [mesh for mesh in self.ds.index.mesh_union]
                else:
                    mesh_ids = [int(ftype.replace("connect", ""))]
                    chunk = chunks[mesh_ids[0] - 1]
                    objs = chunk.objs
                if fname in self.node_fields:
                    field_ind = self.node_fields.index(fname)
                    fdata = ds.variables["vals_nod_var%d" % (field_ind + 1)]
                    for g in objs:
                        ci = g.connectivity_indices - self._INDEX_OFFSET
                        data = fdata[self.ds.step][ci]
                        ind += g.select(selector, data, rv[field], ind)  # caches
                if fname in self.elem_fields:
                    field_ind = self.elem_fields.index(fname)
                    for g, mesh_id in zip(objs, mesh_ids):
                        fdata = ds.variables[
                            "vals_elem_var%deb%s" % (field_ind + 1, mesh_id)
                        ][:]
                        data = fdata[self.ds.step, :]
                        ind += g.select(selector, data, rv[field], ind)  # caches
                rv[field] = rv[field][:ind]
            return rv

    def _read_chunk_data(self, chunk, fields):
        # This reads the data from a single chunk, and is only used for
        # caching.
        pass
