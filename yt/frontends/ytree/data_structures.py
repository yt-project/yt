"""
Data structures for ytree frontend.




"""

#-----------------------------------------------------------------------------
# Copyright (c) yt Development Team. All rights reserved.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.utilities.on_demand_imports import _h5py as h5py
import numpy as np
import json
import os
import stat

from .fields import \
    YTreeFieldInfo

from yt.frontends.halo_catalog.data_structures import \
    HaloCatalogFile
from yt.frontends.ytdata.data_structures import \
    SavedDataset
from yt.geometry.particle_geometry_handler import \
    ParticleIndex

class YTreeParticleIndex(ParticleIndex):
    def _setup_filenames(self):
        template = self.dataset.filename_template
        cls = self.dataset._file_class
        self.data_files = \
          [cls(self.dataset, self.io, template % i, i)
           for i in range(self.dataset.file_count)]

class YTreeHDF5File(HaloCatalogFile):
    def __init__(self, ds, io, filename, file_id):
        with h5py.File(filename, "r") as f:
            self.particle_count = f['data'].attrs['num_elements']
        super(YTreeHDF5File, self).__init__(
            ds, io, filename, file_id)

    def _read_particle_positions(self, ptype, f=None):
        """
        Read all particle positions in this file.
        """

        if f is None:
            close = True
            f = h5py.File(self.filename, "r")
        else:
            close = False

        pos = np.empty((self.particle_count, 3), dtype="float64")
        for i, ax in enumerate('xyz'):
            pos[:, i] = f["data/position_%s" % ax].value

        if close:
            f.close()

        return pos

class YTreeDataset(SavedDataset):
    _index_class = YTreeParticleIndex
    _file_class = YTreeHDF5File
    _field_info_class = YTreeFieldInfo
    _suffix = ".h5"
    _con_attrs = ("hubble_constant", "omega_matter", "omega_lambda")

    def __init__(self, filename, dataset_type="ytree_arbor",
                 n_ref = 16, over_refine_factor = 1, units_override=None,
                 unit_system="cgs"):
        self.n_ref = n_ref
        self.over_refine_factor = over_refine_factor
        super(YTreeDataset, self).__init__(filename, dataset_type,
                                                 units_override=units_override,
                                                 unit_system=unit_system)

    def _set_derived_attrs(self):
        self.domain_center = 0.5 * (self.domain_right_edge +
                                    self.domain_left_edge)
        self.domain_width = self.domain_right_edge - self.domain_left_edge

    def _with_parameter_file_open(self, f):
        self.file_count = f.attrs['total_files']
        self.particle_count = f.attrs['total_nodes']
        self.box_size = self.quan(
            f.attrs['box_size'], f.attrs['box_size_units'])
        self._field_dict = json.loads(f.attrs['field_info'])

    def _parse_parameter_file(self):
        self.current_redshift = None
        self.current_time = None
        self.cosmological_simulation = 1
        self.refine_by = 2
        self.dimensionality = 3
        nz = 1 << self.over_refine_factor
        self.domain_dimensions = np.ones(self.dimensionality, "int32") * nz
        self.periodicity = (True, True, True)
        self.unique_identifier = \
          int(os.stat(self.parameter_filename)[stat.ST_CTIME])
        prefix = self.parameter_filename[:-len(self._suffix)]
        self.filename_template = "%s_%%04d%s" % (prefix, self._suffix)
        self.particle_types = ("halos")
        self.particle_types_raw = ("halos")
        super(YTreeDataset, self)._parse_parameter_file()
        # Set this again because unit registry has been updated.
        self.box_size = self.quan(self.box_size)
        self.domain_left_edge = self.box_size * \
          np.zeros(self.dimensionality)
        self.domain_right_edge = self.box_size * \
          np.ones(self.dimensionality)

    @classmethod
    def _is_valid(self, *args, **kwargs):
        if not args[0].endswith(".h5"): return False
        with h5py.File(args[0], "r") as f:
            if "data_type" in f.attrs and "arbor_type" in f.attrs:
                return True
        return False
