from yt.utilities.lib.cyoctree_raytracing import CythonOctreeRayTracing

import numpy as np
from itertools import product

class OctreeRayTracing(object):
    octree = None
    data_source = None
    log_fields = None
    fields = None

    # Internal data
    _cell_index = None
    _tvalues = None

    def __init__(self, data_source):
        self.data_source = data_source
        LE = np.array([0, 0, 0], dtype=np.float64)
        RE = np.array([1, 1, 1], dtype=np.float64)
        depth = data_source.ds.parameters['levelmax']

        self.octree = CythonOctreeRayTracing(LE, RE, depth)
        ds = data_source.ds

        xyz = np.stack([data_source[_].to('unitary').value for _ in 'x y z'.split()], axis=-1)
        lvl = data_source['grid_level'].astype(int).value + ds.parameters['levelmin']

        ipos = np.floor(xyz * (1<<(ds.parameters['levelmax']))).astype(int)
        self.octree.add_nodes(ipos.astype(np.int32), lvl.astype(np.int32), np.arange(len(ipos), dtype=np.int32))

    def vertex_centered_data(self, field):
        data_source = self.data_source
        chunks = data_source.index._chunk(data_source, "spatial", ngz=1)

        finfo = data_source.ds._get_field_info(*field)
        units = finfo.units
        rv = data_source.ds.arr(np.zeros((2, 2, 2, data_source.ires.size), dtype="float64"), units)
        ind = {(i, j, k): 0 for i, j, k in product(*[range(2)]*3)}
        for chunk in chunks:
            with data_source._chunked_read(chunk):
                gz = data_source._current_chunk.objs[0]
                gz.field_parameters = data_source.field_parameters
                wogz = gz._base_grid
                vertex_data = gz.get_vertex_centered_data([field])[field]

                for i, j, k in product(*[range(2)]*3):
                    ind[i, j, k] += wogz.select(
                        data_source.selector,
                        vertex_data[i:i+2, j:j+2, k:k+2, ...],
                        rv[i, j, k, :], ind[i, j, k])
        return rv

    def set_fields(self, fields, log_fields, no_ghost, force=False):
        if no_ghost:
            raise NotImplementedError('Ghost zones are required with Octree datasets')

        assert len(fields) == 1
        field = self.data_source._determine_fields(fields)[0]
        take_log = log_fields[0]
        vertex_data = self.vertex_centered_data(field)

        if take_log:
            vertex_data = np.log10(vertex_data)

        # Vertex_data has shape (2, 2, 2, ...)
        self.data = vertex_data.reshape(8, -1)

    def cast_rays(self, vp_pos, vp_dir):
        """Cast the rays through the oct.

        Parameters
        ----------
        vp_pos, vp_dir : float arrays (Nrays, Ndim)
            The position (unitary) and direction of each ray

        Returns
        -------
        cell_index : list of integer arrays of shape (Ncell)
            For each ray, contains an ordered array of cell ids
            that it intersects with
        tvalues : list of float arrays of shape (Ncell, 2)
            The t value at entry and exit for each cell.
        """
        if not self._cell_index:
            # TODO: cache indices of cells
            self._cell_index, self._tvalues = \
                self.octree.cast_rays(vp_pos, vp_dir)
        return self._cell_index, self._tvalues
