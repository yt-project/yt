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
    fill_fcoords, fill_fwidths

from yt.data_objects.data_containers import \
    YTFieldData, \
    YTSelectionContainer
import yt.geometry.particle_deposit as particle_deposit

class UnstructuredMesh(YTSelectionContainer):
    # This is a base class, not meant to be used directly.
    _spatial = False
    _connectivity_length = -1
    _type_name = 'unstructured_mesh'
    _skip_add = True
    _index_offset = 0
    _con_args = ('mesh_id', 'filename', 'connectivity_indices',
                 'connectivity_coords')

    def __init__(self, mesh_id, filename, connectivity_indices,
                 connectivity_coords, index):
        self.field_data = YTFieldData()
        self.filename = filename
        self.field_parameters = {}
        self.mesh_id = mesh_id
        # This is where we set up the connectivity information
        self.connectivity_indices = connectivity_indices
        if connectivity_indices.shape[1] != self._connectivity_length:
            if self._connectivity_length == -1:
                self._connectivity_length = connectivity_indices.shape[1]
            else:
                raise RuntimeError
        self.connectivity_coords = connectivity_coords
        self.ds = index.dataset
        self._index = index
        self._last_mask = None
        self._last_count = -1
        self._last_selector_id = None
        self._current_particle_type = 'all'
        self._current_fluid_type = self.ds.default_fluid_type

    def _check_consistency(self):
        if self.connectivity_indices.shape[1] != self._connectivity_length:
            raise RuntimeError

        for gi in range(self.connectivity_indices.shape[0]):
            ind = self.connectivity_indices[gi, :] - self._index_offset
            coords = self.connectivity_coords[ind, :]
            for i in range(3):
                assert(np.unique(coords[:,i]).size == 2)
        mylog.debug("Connectivity is consistent.")

    def __repr__(self):
        return "UnstructuredMesh_%04i" % (self.mesh_id)

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

    def _generate_container_field(self, field):
        raise NotImplementedError

    def select_fcoords(self, dobj = None):
        # This computes centroids!
        mask = self._get_selector_mask(dobj.selector)
        if mask is None: return np.empty((0,3), dtype='float64')
        centers = fill_fcoords(self.connectivity_coords,
                               self.connectivity_indices,
                               self._index_offset)
        return centers[mask, :]

    def select_fwidth(self, dobj):
        raise NotImplementedError

    def select_icoords(self, dobj):
        raise NotImplementedError

    def select_ires(self, dobj):
        raise NotImplementedError

    def select_tcoords(self, dobj):
        raise NotImplementedError

    def deposit(self, positions, fields = None, method = None,
                kernel_name = 'cubic'):
        raise NotImplementedError
        # Here we perform our particle deposition.
        cls = getattr(particle_deposit, "deposit_%s" % method, None)
        if cls is None:
            raise YTParticleDepositionNotImplemented(method)
        # We allocate number of zones, not number of octs
        op = cls(self.ActiveDimensions.prod(), kernel_name)
        op.initialize()
        op.process_grid(self, positions, fields)
        vals = op.finalize()
        if vals is None: return
        return vals.reshape(self.ActiveDimensions, order="C")

    def select_blocks(self, selector):
        mask = self._get_selector_mask(selector)
        yield self, mask

    def select(self, selector, source, dest, offset):
        mask = self._get_selector_mask(selector)
        count = self.count(selector)
        if count == 0: return 0
        dest[offset:offset+count] = source[mask, ...]
        return count

    def count(self, selector):
        mask = self._get_selector_mask(selector)
        if mask is None: return 0
        return self._last_count

    def count_particles(self, selector, x, y, z):
        # We don't cache the selector results
        count = selector.count_points(x,y,z, 0.0)
        return count

    def select_particles(self, selector, x, y, z):
        mask = selector.select_points(x,y,z, 0.0)
        return mask

    def _get_selector_mask(self, selector):
        if hash(selector) == self._last_selector_id:
            mask = self._last_mask
        else:
            self._last_mask = mask = selector.fill_mesh_cell_mask(self)
            self._last_selector_id = hash(selector)
            if mask is None:
                self._last_count = 0
            else:
                self._last_count = mask.sum()
        return mask

    def select_fcoords_vertex(self, dobj = None):
        mask = self._get_selector_mask(dobj.selector)
        if mask is None: return np.empty((0, self._connectivity_length, 3), dtype='float64')
        vertices = self.connectivity_coords[
                self.connectivity_indices - 1]
        return vertices[mask, :, :]


class SemiStructuredMesh(UnstructuredMesh):
    _connectivity_length = 8
    _type_name = 'semi_structured_mesh'
    _container_fields = ("dx", "dy", "dz")

    def __repr__(self):
        return "SemiStructuredMesh_%04i" % (self.mesh_id)

    def _generate_container_field(self, field):
        if self._current_chunk is None:
            self.index._identify_base_chunk(self)
        if field == "dx":
            return self._current_chunk.fwidth[:,0]
        elif field == "dy":
            return self._current_chunk.fwidth[:,1]
        elif field == "dz":
            return self._current_chunk.fwidth[:,2]

    def select_fwidth(self, dobj):
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

    def _get_selector_mask(self, selector):
        if hash(selector) == self._last_selector_id:
            mask = self._last_mask
        else:
            self._last_mask = mask = selector.fill_mesh_cell_mask(self)
            self._last_selector_id = hash(selector)
            if mask is None:
                self._last_count = 0
            else:
                self._last_count = mask.sum()
        return mask

    def select(self, selector, source, dest, offset):
        mask = self._get_selector_mask(selector)
        count = self.count(selector)
        if count == 0: return 0
        # Note: this likely will not work with vector fields.
        dest[offset:offset+count] = source.flat[mask]
        return count

