"""
Unstructured-mesh geometry handler




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
import os
import weakref

from yt.utilities.logger import ytLogger as mylog
from yt.geometry.geometry_handler import Index, YTDataChunk
from yt.utilities.lib.mesh_utilities import smallest_fwidth

class UnstructuredIndex(Index):
    """The Index subclass for unstructured and hexahedral mesh datasets. """
    _global_mesh = False
    _unsupported_objects = ('proj', 'covering_grid', 'smoothed_covering_grid')

    def __init__(self, ds, dataset_type):
        self.dataset_type = dataset_type
        self.dataset = weakref.proxy(ds)
        self.index_filename = self.dataset.parameter_filename
        self.directory = os.path.dirname(self.index_filename)
        self.float_type = np.float64
        super(UnstructuredIndex, self).__init__(ds, dataset_type)

    def _setup_geometry(self):
        mylog.debug("Initializing Unstructured Mesh Geometry Handler.")
        self._initialize_mesh()

    def get_smallest_dx(self):
        """
        Returns (in code units) the smallest cell size in the simulation.
        """
        dx = min(smallest_fwidth(mesh.connectivity_coords,
                                 mesh.connectivity_indices,
                                 mesh._index_offset)
                 for mesh in self.meshes)
        return dx

    def convert(self, unit):
        return self.dataset.conversion_factors[unit]

    def _initialize_mesh(self):
        raise NotImplementedError

    def _identify_base_chunk(self, dobj):
        if getattr(dobj, "_chunk_info", None) is None:
            dobj._chunk_info = self.meshes
        if getattr(dobj, "size", None) is None:
            dobj.size = self._count_selection(dobj)
        dobj._current_chunk = list(self._chunk_all(dobj))[0]

    def _count_selection(self, dobj, meshes = None):
        if meshes is None: meshes = dobj._chunk_info
        count = sum((m.count(dobj.selector) for m in meshes))
        return count

    def _chunk_all(self, dobj, cache = True):
        oobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        yield YTDataChunk(dobj, "all", oobjs, dobj.size, cache)

    def _chunk_spatial(self, dobj, ngz, sort = None, preload_fields = None):
        sobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        # We actually do not really use the data files except as input to the
        # ParticleOctreeSubset.
        # This is where we will perform cutting of the Octree and
        # load-balancing.  That may require a specialized selector object to
        # cut based on some space-filling curve index.
        for i,og in enumerate(sobjs):
            if ngz > 0:
                g = og.retrieve_ghost_zones(ngz, [], smoothed=True)
            else:
                g = og
            size = self._count_selection(dobj, [og])
            if size == 0: continue
            yield YTDataChunk(dobj, "spatial", [g], size)

    def _chunk_io(self, dobj, cache = True, local_only = False):
        oobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        for subset in oobjs:
            s = self._count_selection(dobj, oobjs)
            yield YTDataChunk(dobj, "io", [subset], s, cache = cache)
