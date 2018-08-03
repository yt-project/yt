"""
Data structures for Cholla.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2018, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
import os
import weakref
import glob
import h5py

from yt.funcs import \
    mylog, \
    ensure_tuple
from yt.data_objects.grid_patch import \
    AMRGridPatch
from yt.geometry.grid_geometry_handler import \
    GridIndex
from yt.data_objects.static_output import \
    Dataset
from yt.utilities.lib.misc_utilities import \
    get_box_grids_level
from yt.geometry.geometry_handler import \
    YTDataChunk
from yt.utilities.file_handler import \
    HDF5FileHandler

from .fields import ChollaFieldInfo
from yt.utilities.decompose import \
    decompose_array, get_psize

class ChollaGrid(AMRGridPatch):
    _id_offset = 0

    def __init__(self, id, index, level):
        super(ChollaGrid, self).__init__(self, id, \
                                         filename=index.index_filename, \
                                         index=index)
        self.filename = index.index_filename
        self.Parent = None
        self.Children = []
        self.Level = level

    def _setup_dx(self):
        # So first we figure out what the index is.  We don't assume
        # that dx=dy=dz , at least here.  We probably do elsewhere.
        id = self.id - self._id_offset
        if len(self.Parent) > 0:
            self.dds = self.Parent[0].dds / self.ds.refine_by
        else:
            LE, RE = self.index.grid_left_edge[id,:], \
                     self.index.grid_right_edge[id,:]
            self.dds = self.ds.arr((RE-LE)/self.ActiveDimensions, "code_length")
        if self.ds.dimensionality < 2: self.dds[1] = 1.0
        if self.ds.dimensionality < 3: self.dds[2] = 1.0
        self.field_data['dx'], self.field_data['dy'], self.field_data['dz'] = self.dds

    def __repr__(self):
        return "ChollaGrid_%04i (%s)" % (self.id, self.ActiveDimensions)

class ChollaHierarchy(GridIndex):

    grid = ChollaGrid
    _dataset_type='cholla'
    _data_file = None

    def __init__(self, ds, dataset_type='cholla'):
        self.dataset = weakref.proxy(ds)
        self.directory = os.path.dirname(self.dataset.filename)
        self.dataset_type = dataset_type
        # for now, the index file is the dataset!
        self.index_filename = self.dataset.filename
        self._handle = ds._handle
        GridIndex.__init__(self, ds, dataset_type)

    def _parse_index(self):
        self.num_grids = 1
        self.grid_left_edge = self.dataset.domain_left_edge
        self.grid_right_edge = self.dataset.domain_right_edge
        self.grid_dimensions = self.dataset.dimensions
        self.grid_particle_count = np.ones([self.num_grids, 1], dtype='int64')
        self.grid_levels = np.ones([self.num_grids, 1], dtype='int64')
        self.max_level = np.ones([self.num_grids, 1], dtype='int64')
        self.grids = np.empty(self.num_grids, dtype='object')
        for i in range(num_grids):
            self.grids[i] = self.grid(i, self, levels[i])

    def _populate_grid_objects(self):
        for g in self.grids:
            g._prepare_grid()
            g._setup_dx()
        self.max_level = self.grid_levels.max()

    def _chunk_io(self, dobj, cache = True, local_only = False):
        gobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        for subset in gobjs:
            yield YTDataChunk(dobj, "io", [subset],
                              self._count_selection(dobj, [subset]),
                              cache = cache)

class ChollaDataset(Dataset):
    _index_class = ChollaHierarchy
    _field_info_class = ChollaFieldInfo
    _dataset_type = "cholla"

    def __init__(self, filename, dataset_type='cholla',
                 storage_filename=None, parameters=None,
                 units_override=None, nprocs=1, unit_system="code"):
        self.fluid_types += ("cholla",)
        if parameters is None:
            parameters = {}
        if units_override is None:
            units_override = {}
        self._handle = HDF5FileHandler(filename)
        Dataset.__init__(self, filename, dataset_type, units_override=units_override,
                         unit_system=unit_system)
        self.filename = filename
        if storage_filename is None:
            storage_filename = '%s.yt' % filename.split('/')[-1]
        self.storage_filename = storage_filename
        self.backup_filename = self.filename[:-4] + "_backup.gdf"

    def _set_code_unit_attributes(self):
        """
        Sets the code units
        """
        self.length_unit = self.quan(1.0, 'kpc')
        self.mass_unit = self.quan(1.0, 'Msun')
        self.time_unit = self.quan(1.0, 'kyr')
        self.velocity_unit = self.length_unit / self.time_unit

    def _parse_parameter_file(self):
        self.unique_identifier = self.parameter_filename.__hash__()
        self.parameters = self._handle.attrs
        self.domain_left_edge = self._handle.attrs['offset']
        self.domain_right_edge = self.domain_left_edge + \
                                 self._handle.attrs['domain']
        self.domain_dimensions = self._handle.attrs["dims"]
        self.dimensionality = len(self.domain_dimensions)
        self.current_time = self._handle.attrs["t"]
        self.cosmological_simulation = False
        self.periodicity = (False, False, False)
        self.gamma = self._handle.attrs['gamma']
        self._handle.close()

    @classmethod
    def _is_valid(self, *args, **kwargs):
        # Tag as Cholla if 'n_step' attribute in h5 file root level
        # Need something better than this
        try:
            fh = h5py.File(args[0], mode='r')
            if 'n_step' in fh.attrs:
                return True
        except:
            pass
        return False

    @property
    def _skip_cache(self):
        return True

    def __repr__(self):
        return self.basename.rsplit(".", 1)[0]
