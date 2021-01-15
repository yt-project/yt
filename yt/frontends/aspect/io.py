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
        ftype_list = []
        rv = {}
        for field in fields:
            ftype, fname = field
            rv[field] = []
            ftype_list.append(ftype)

        for chunk_id, chunk in enumerate(chunks):
            mesh = chunk.objs[0]
            el_counts = np.array(mesh.element_count)
            # operations here should not offset mesh indeces as the indices for
            # a chunk's mesh are not global. need to temporarily
            # zero the mesh index for use within the selector objects, then
            # reset after for places where global concatentation of the indices
            # is required.
            # orig_offset = mesh._index_offset
            # mesh._index_offset = 0

            # for now, we have all file info in single chunk
            mesh_name = f"connect{chunk_id+1}"
            if "all" in ftype_list or mesh_name in ftype_list:
                # mask here is the **element** mask. i.e., shape(mask) = (n_elments,)
                # rather than (n_elements, n_verts). These apply to all fields, so
                # pull them out here:
                mask = selector.fill_mesh_cell_mask(mesh)
                if mask is not None:
                    # each mesh piece will pull the connectivity and mask values
                    # for the indices corresponding to each vtu file
                    for vtu_id, vtu_filename in enumerate(mesh.filenames):
                        # load this vtu's part of the mesh
                        vtu_file = os.path.join(self.ds.data_dir, vtu_filename)
                        meshPiece = meshio.read(vtu_file)

                        # the connectivty for this piece
                        cell_types = list(meshPiece.cells_dict.keys())
                        vtu_con = np.array(meshPiece.cells_dict[cell_types[0]])
                        con_shape = vtu_con.shape
                        conn1d = vtu_con.ravel()

                        # the element mask for this piece
                        element_offset_start = el_counts[0:vtu_id].sum()
                        element_offset_end = element_offset_start + el_counts[vtu_id]
                        vtu_mask = mask[element_offset_start:element_offset_end]

                        for field in fields:
                            ftype, fname = field
                            if ftype == "all" or ftype == mesh_name:
                                # meshio returns a 1d data array, so we need to:
                                # 1. select with 1d connectivity to ensure proper order
                                # 2. reshape to expected (element, n_verts) shape
                                # 3. select the masked data
                                # field_data = meshPiece.point_data[fname]
                                if "velocity" in fname:
                                    vdim = fname.split("_")[-1]
                                    vdim_to_dim = {"x": 0, "y": 1, "z": 2}
                                    dim = vdim_to_dim[vdim]
                                    raw_data = meshPiece.point_data["velocity"][:, dim]
                                else:
                                    raw_data = meshPiece.point_data[fname]

                                raw_data = raw_data[conn1d].reshape(con_shape)
                                rv[field].append(raw_data[vtu_mask, :])

        for field in fields:
            rv[field] = np.concatenate(rv[field]).astype(np.float64)

        return rv

    def _read_chunk_data(self, chunk, fields):
        # This reads the data from a single chunk, and is only used for
        # caching.
        pass
