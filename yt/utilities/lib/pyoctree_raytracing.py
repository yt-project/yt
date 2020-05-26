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
            raise NotImplementedError('Cannot use no ghost with Octree datasets')
        new_fields = self.data_source._determine_fields(fields)
        regenerate_data = self.fields is None or \
                          len(self.fields) != len(new_fields) or \
                          self.fields != new_fields or force
        if not iterable(log_fields):
            log_fields = [log_fields]
        new_log_fields = list(log_fields)
        self.fields = new_fields

        if self.log_fields is not None and not regenerate_data:
            flip_log = list(map(operator.ne, self.log_fields, new_log_fields))
        else:
            flip_log = [False] * len(new_log_fields)
        self.log_fields = new_log_fields

        # TODO: cache data in the 3x3x3 neighbouring cells

    def cast_rays(self, vp_pos, vp_dir):
        # TODO: cache indices of cells
        self.cell_index, self.tvalues = \
            self.octree.cast_rays(vp_pos, vp_dir)

    def sample(self, sampler):
        # TODO: Apply to sampler to each oct encountered by all rays.
        pass

    def traverse(self, viewpoint):
        raise Exception()
        self.octree.cast_rays()