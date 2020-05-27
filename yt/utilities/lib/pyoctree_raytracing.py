from yt.funcs import iterable

from yt.utilities.lib.cyoctree_raytracing import CythonOctreeRayTracing
from yt.utilities.amr_kdtree.amr_kdtree import _apply_log
import numpy as np

import operator



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

    def set_fields(self, fields, log_fields, no_ghost, force=False):
        if no_ghost:
            raise NotImplementedError('Ghost zones are required with Octree datasets')

        assert len(fields) == 1
        fields = self.data_source._determine_fields(fields)

        for field, take_log in zip(fields, log_fields):
            if take_log:
                tmp = np.log10(self.data_source[field])
            else:
                tmp = self.data_source[field]
            self.data = [tmp]*8

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

    # def sample(self, sampler):
    #     # TODO: Apply to sampler to each oct encountered by all rays.
    #     self.octree.sample(sampler, self._cell_index, self._tvalues)

    def traverse(self, viewpoint):
        raise Exception()
        self.octree.cast_rays()