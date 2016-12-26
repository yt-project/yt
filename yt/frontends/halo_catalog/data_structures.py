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
    
class HaloCatalogDataset(Dataset):
    _index_class = ParticleIndex
    _file_class = HaloCatalogHDF5File
    _field_info_class = HaloCatalogFieldInfo
    _suffix = ".h5"

    def __init__(self, filename, dataset_type="halocatalog_hdf5",
                 n_ref = 16, over_refine_factor = 1, units_override=None,
                 unit_system="cgs"):
        self.n_ref = n_ref
        self.over_refine_factor = over_refine_factor
        super(HaloCatalogDataset, self).__init__(filename, dataset_type,
                                                 units_override=units_override,
                                                 unit_system=unit_system)

    def _parse_parameter_file(self):
        with h5py.File(self.parameter_filename, "r") as f:
            hvals = dict((key, parse_h5_attr(f, key)) for key in f.attrs.keys())
        self.parameters.update(hvals)
        self.dimensionality = 3
        self.refine_by = 2
        self.unique_identifier = \
            int(os.stat(self.parameter_filename)[stat.ST_CTIME])
        prefix = ".".join(self.parameter_filename.rsplit(".", 2)[:-2])
        self.filename_template = "%s.%%(num)s%s" % (prefix, self._suffix)
        self.file_count = len(glob.glob(prefix + "*" + self._suffix))

        # if saved, restore unit registry from the json string
        if "unit_registry_json" in self.parameters:
            self.unit_registry = UnitRegistry.from_json(
                self.parameters["unit_registry_json"])
            del self.parameters["unit_registry_json"]
            # reset self.arr and self.quan to use new unit_registry
            self._arr = None
            self._quan = None
        self._assign_unit_system("cgs")

        # assign units to parameters that have associated unit string
        del_pars = []
        for par in self.parameters:
            ustr = "%s_units" % par
            if ustr in self.parameters:
                if isinstance(self.parameters[par], np.ndarray):
                    to_u = self.arr
                else:
                    to_u = self.quan
                self.parameters[par] = to_u(
                    self.parameters[par], self.parameters[ustr])
                del_pars.append(ustr)
        for par in del_pars:
            del self.parameters[par]

        for attr in ["cosmological_simulation", "current_time", "current_redshift",
                     "hubble_constant", "omega_matter", "omega_lambda",
                     "domain_left_edge", "domain_right_edge"]:
            setattr(self, attr, self.parameters[attr])
        self.periodicity = (True, True, True)
        self.particle_types = ("halos")
        self.particle_types_raw = ("halos")

        nz = 1 << self.over_refine_factor
        self.domain_dimensions = np.ones(3, "int32") * nz

    def set_units(self):
        if "unit_registry_json" in self.parameters:
            self._set_code_unit_attributes()
        else:
            super(HaloCatalogDataset, self).set_units()

    def _set_code_unit_attributes(self):
        attrs = ('length_unit', 'mass_unit', 'time_unit',
                 'velocity_unit', 'magnetic_unit')
        cgs_units = ('cm', 'g', 's', 'cm/s', 'gauss')
        base_units = np.ones(len(attrs))
        for unit, attr, cgs_unit in zip(base_units, attrs, cgs_units):
            if attr in self.parameters and \
              isinstance(self.parameters[attr], YTQuantity):
                uq = self.parameters[attr]
            elif attr in self.parameters and \
              "%s_units" % attr in self.parameters:
                uq = self.quan(self.parameters[attr],
                               self.parameters["%s_units" % attr])
                del self.parameters[attr]
                del self.parameters["%s_units" % attr]
            elif isinstance(unit, string_types):
                uq = self.quan(1.0, unit)
            elif isinstance(unit, numeric_type):
                uq = self.quan(unit, cgs_unit)
            elif isinstance(unit, YTQuantity):
                uq = unit
            elif isinstance(unit, tuple):
                uq = self.quan(unit[0], unit[1])
            else:
                raise RuntimeError("%s (%s) is invalid." % (attr, unit))
            setattr(self, attr, uq)

    @classmethod
    def _is_valid(self, *args, **kwargs):
        if not args[0].endswith(".h5"): return False
        with h5py.File(args[0], "r") as f:
            if "data_type" in f.attrs and \
              parse_h5_attr(f, "data_type") == "halo_catalog":
                return True
        return False
