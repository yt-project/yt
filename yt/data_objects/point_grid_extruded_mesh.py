"""
Unstructured mesh base container.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.funcs import mylog
from yt.utilities.exceptions import \
    YTParticleDepositionNotImplemented
from yt.utilities.lib.mesh_utilities import \
    fill_fcoords_extruded

from yt.data_objects.data_containers import \
    YTSelectionContainer
import yt.geometry.particle_deposit as particle_deposit

class PointGridExtrudedMesh(YTSelectionContainer):
    # This particular format expects to have a grid of x, y points that define
    # the vertices, coupled with a "distance" that defines the ranges to the
    # different grid points.  We assume that we can connect along these grids.
    # For instance, this might be the form of an elevation-azimuth-range
    # dataset.
    _connectivity_length = 8
    _type_name = "point_grid_extruded_mesh"
    _container_fields = ("dx", "dy", "dz")
    _spatial = False
    _skip_add = True
    _index_offset = 0
    _con_args = ('mesh_id', 'filename', '_grid_locations', '_grid_extrusion')
    def __init__(self, mesh_id, filename, xy_grid, extruded_grid, index):
        super(PointGridExtrudedMesh, self).__init__(index.dataset, None)
        self.filename = filename
        self.mesh_id = mesh_id
        # This is where we set up the connectivity information
        self._grid_locations = xy_grid.copy()
        self._grid_extrusion = extruded_grid.copy()
        self.ds = index.dataset
        self._index = index
        self._last_mask = None
        self._last_count = -1
        self._last_selector_id = None

    def get_global_startindex(self):
        """
        Return the integer starting index for each dimension at the current
        level.

        """
        raise NotImplementedError

    def convert(self, datatype):
        """
        This will attempt to convert a given unit to cgs from code units. It
        either returns the multiplicative factor or throws a KeyError.

        """
        return self.ds[datatype]

    @property
    def shape(self):
        raise NotImplementedError

    def __repr__(self):
        return "PointGridExtrudedMesh_%04i" % (self.mesh_id)

    def select_fcoords(self, dobj = None):
        # This computes centroids!
        mask = self._get_selector_mask(dobj.selector)
        if mask is None: return np.empty((0,3), dtype='float64')
        centers = fill_fcoords_extruded(self._grid_locations,
                                        self._grid_extrusion)
        return centers[mask, :]

    def select_fwidth(self, dobj):
        # This is a placeholder; it's not immediately obvious to me if we will
        # be able to implement a meaningful fwidth, since the data points will
        # be irregular in some dimensions.  A width is not necessarily
        # "meaningful" when discussing the types of polyhedra.
        raise NotImplementedError
        mask = self._get_selector_mask(dobj.selector)
        if mask is None: return np.empty((0,3), dtype='float64')
        widths = fill_fwidths(self.connectivity_coords,
                              self.connectivity_indices,
                              self._index_offset)
        return widths[mask, :]

    def select_ires(self, dobj):
        ind = np.zeros(self.connectivity_indices.shape[0])
        mask = self._get_selector_mask(dobj.selector)
        if mask is None: return np.empty(0, dtype='int32')
        return ind[mask]

    def select_tcoords(self, dobj):
        mask = self._get_selector_mask(dobj.selector)
        if mask is None: return np.empty(0, dtype='float64')
        dt, t = dobj.selector.get_dt_mesh(self, mask.sum(), self._index_offset)
        return dt, t

    def _generate_container_field(self, field):
        if self._current_chunk is None:
            self.index._identify_base_chunk(self)
        if field == "dx":
            return self._current_chunk.fwidth[:,0]
        elif field == "dy":
            return self._current_chunk.fwidth[:,1]
        elif field == "dz":
            return self._current_chunk.fwidth[:,2]

    def select(self, selector, source, dest, offset):
        mask = self._get_selector_mask(selector)
        count = self.count(selector)
        if count == 0: return 0
        # Note: this likely will not work with vector fields.
        dest[offset:offset+count] = source.flat[mask]
        return count

