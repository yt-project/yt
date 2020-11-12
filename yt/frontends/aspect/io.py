import os

import meshio  # should be on demand
import numpy as np

from yt.utilities.io_handler import BaseIOHandler


class IOHandlerASPECT(BaseIOHandler):
    _particle_reader = False
    _dataset_type = "aspect"
    _INDEX_OFFSET = 1

    def __init__(self, ds):
        self.filename = ds.index_filename
        # exodus_ii_handler = NetCDF4FileHandler(self.filename)
        # self.handler = exodus_ii_handler
        super(IOHandlerASPECT, self).__init__(ds)
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

        # always select all since these are chunks of a single mesh...
        rv = {}
        for field in fields:
            rv[field] = []

        for chunk in chunks:
            mesh = chunk.objs[0]
            mask = selector.fill_mesh_cell_mask(mesh)
            masked_conn = mesh.connectivity_indices[mask, :].ravel()
            vtu_file = os.path.join(self.ds.data_dir, mesh.filename)
            print(vtu_file)
            meshPiece = meshio.read(vtu_file)
            for field in fields:
                ftype, fname = field
                rv[field].append(meshPiece.point_data[fname][masked_conn])

        for field in fields:
            rv[field] = np.concatenate(rv[field])

        return rv

    def _read_chunk_data(self, chunk, fields):
        # This reads the data from a single chunk, and is only used for
        # caching.
        pass
