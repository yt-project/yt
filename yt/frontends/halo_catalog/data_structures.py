"""
Data structures for HaloCatalog frontend.




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.utilities.on_demand_imports import _h5py as h5py
from numbers import \
    Number as numeric_type
import numpy as np
import stat
import glob
import os

from .fields import \
    HaloCatalogFieldInfo

from yt.frontends.ytdata.data_structures import \
    SavedDataset
from yt.extern.six import \
    string_types
from yt.funcs import \
    parse_h5_attr
from yt.geometry.particle_geometry_handler import \
    ParticleIndex
from yt.data_objects.static_output import \
    Dataset, \
    ParticleFile
from yt.units.unit_registry import \
    UnitRegistry
from yt.units.yt_array import \
    YTQuantity

class HaloCatalogHDF5File(ParticleFile):
    def __init__(self, ds, io, filename, file_id):
        with h5py.File(filename, "r") as f:
            self.header = dict((field, parse_h5_attr(f, field)) \
                               for field in f.attrs.keys())

        super(HaloCatalogHDF5File, self).__init__(ds, io, filename, file_id)
    
class HaloCatalogDataset(SavedDataset):
    _index_class = ParticleIndex
    _file_class = HaloCatalogHDF5File
    _field_info_class = HaloCatalogFieldInfo
    _suffix = ".h5"
    _con_attrs = ("cosmological_simulation",
                  "current_time", "current_redshift",
                  "hubble_constant", "omega_matter", "omega_lambda",
                  "domain_left_edge", "domain_right_edge")

    def __init__(self, filename, dataset_type="halocatalog_hdf5",
                 n_ref = 16, over_refine_factor = 1, units_override=None,
                 unit_system="cgs"):
        self.n_ref = n_ref
        self.over_refine_factor = over_refine_factor
        super(HaloCatalogDataset, self).__init__(filename, dataset_type,
                                                 units_override=units_override,
                                                 unit_system=unit_system)

    def _parse_parameter_file(self):
        self.refine_by = 2
        self.dimensionality = 3
        nz = 1 << self.over_refine_factor
        self.domain_dimensions = np.ones(self.dimensionality, "int32") * nz
        self.periodicity = (True, True, True)
        prefix = ".".join(self.parameter_filename.rsplit(".", 2)[:-2])
        self.filename_template = "%s.%%(num)s%s" % (prefix, self._suffix)
        self.file_count = len(glob.glob(prefix + "*" + self._suffix))
        self.particle_types = ("halos")
        self.particle_types_raw = ("halos")
        super(HaloCatalogDataset, self)._parse_parameter_file()

    @classmethod
    def _is_valid(self, *args, **kwargs):
        if not args[0].endswith(".h5"): return False
        with h5py.File(args[0], "r") as f:
            if "data_type" in f.attrs and \
              parse_h5_attr(f, "data_type") == "halo_catalog":
                return True
        return False
