"""
FLASH-specific data structures

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: UCSD
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2010 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import h5py
import stat
import numpy as na
import weakref

from yt.funcs import *
from yt.data_objects.grid_patch import \
    AMRGridPatch
from yt.data_objects.hierarchy import \
    AMRHierarchy
from yt.data_objects.static_output import \
    StaticOutput
from yt.utilities.definitions import \
    mpc_conversion
from yt.utilities.io_handler import \
    io_registry

from .fields import \
    FLASHFieldContainer, \
    add_field

class FLASHGrid(AMRGridPatch):
    _id_offset = 1
    #__slots__ = ["_level_id", "stop_index"]
    def __init__(self, id, hierarchy, level):
        AMRGridPatch.__init__(self, id, filename = hierarchy.hierarchy_filename,
                              hierarchy = hierarchy)
        self.Parent = None
        self.Children = []
        self.Level = level

    def __repr__(self):
        return "FLASHGrid_%04i (%s)" % (self.id, self.ActiveDimensions)

class FLASHHierarchy(AMRHierarchy):

    grid = FLASHGrid
    _handle = None
    
    def __init__(self,pf,data_style='chombo_hdf5'):
        self.data_style = data_style
        self.field_info = FLASHFieldContainer()
        self.field_indexes = {}
        self.parameter_file = weakref.proxy(pf)
        # for now, the hierarchy file is the parameter file!
        self.hierarchy_filename = self.parameter_file.parameter_filename
        self.directory = os.path.dirname(self.hierarchy_filename)
        self._handle = h5py.File(self.hierarchy_filename)

        self.float_type = na.float64
        AMRHierarchy.__init__(self,pf,data_style)

        self._handle.close()
        self._handle = None

    def _initialize_data_storage(self):
        pass

    def _detect_fields(self):
        ncomp = self._handle["/unknown names"].shape[0]
        self.field_list = [s for s in self._handle["/unknown names"][:].flat]
        if ("/particle names" in self._handle) :
            self.field_list += ["particle_" + s[0].strip() for s
                                in self._handle["/particle names"][:]]
    
    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        AMRHierarchy._setup_classes(self, dd)
        self.object_types.sort()

    def _count_grids(self):
        try:
            self.num_grids = self.parameter_file._find_parameter(
                "integer", "globalnumblocks", True, self._handle)
        except KeyError:
            self.num_grids = self._handle["/simulation parameters"][0][0]
        
    def _parse_hierarchy(self):
        f = self._handle # shortcut
        pf = self.parameter_file # shortcut
        
        self.grid_left_edge[:] = f["/bounding box"][:,:,0]
        self.grid_right_edge[:] = f["/bounding box"][:,:,1]
        # Move this to the parameter file
        try:
            nxb = pf._find_parameter("integer", "nxb", True, f)
            nyb = pf._find_parameter("integer", "nyb", True, f)
            nzb = pf._find_parameter("integer", "nzb", True, f)
        except KeyError:
            nxb, nyb, nzb = [int(f["/simulation parameters"]['n%sb' % ax])
                              for ax in 'xyz']
        self.grid_dimensions[:] *= (nxb, nyb, nzb)
        try:
            self.grid_particle_count[:] = f["/localnp"][:][:,None]
        except KeyError:
            self.grid_particle_count[:] = 0.0
        self._particle_indices = na.zeros(self.num_grids + 1, dtype='int64')
        na.add.accumulate(self.grid_particle_count, out=self._particle_indices[1:])
        # This will become redundant, as _prepare_grid will reset it to its
        # current value.  Note that FLASH uses 1-based indexing for refinement
        # levels, but we do not, so we reduce the level by 1.
        self.grid_levels.flat[:] = f["/refine level"][:][:] - 1
        g = [self.grid(i+1, self, self.grid_levels[i,0])
                for i in xrange(self.num_grids)]
        self.grids = na.array(g, dtype='object')

    def _populate_grid_objects(self):
        # We only handle 3D data, so offset is 7 (nfaces+1)
        
        offset = 7
        ii = na.argsort(self.grid_levels.flat)
        gid = self._handle["/gid"][:]
        for g in self.grids[ii].flat:
            gi = g.id - g._id_offset
            # FLASH uses 1-indexed group info
            g.Children = [self.grids[i - 1] for i in gid[gi,7:] if i > -1]
            for g1 in g.Children:
                g1.Parent = g
            g._prepare_grid()
            g._setup_dx()
        self.max_level = self.grid_levels.max()

    def _setup_unknown_fields(self):
        for field in self.field_list:
            if field in self.parameter_file.field_info: continue
            pfield = field.startswith("particle_")
            mylog.info("Adding %s to list of fields", field)
            cf = None
            if self.parameter_file.has_key(field):
                def external_wrapper(f):
                    def _convert_function(data):
                        return data.convert(f)
                    return _convert_function
                cf = external_wrapper(field)
            add_field(field, lambda a, b: None,
                      convert_function=cf, take_log=False,
                      particle_type=pfield)

    def _setup_derived_fields(self):
        self.derived_field_list = []
        for field in self.parameter_file.field_info:
            try:
                fd = self.parameter_file.field_info[field].get_dependencies(
                            pf = self.parameter_file)
            except:
                continue
            available = na.all([f in self.field_list for f in fd.requested])
            if available: self.derived_field_list.append(field)
        for field in self.field_list:
            if field not in self.derived_field_list:
                self.derived_field_list.append(field)

    def _setup_data_io(self):
        self.io = io_registry[self.data_style](self.parameter_file)

class FLASHStaticOutput(StaticOutput):
    _hierarchy_class = FLASHHierarchy
    _fieldinfo_class = FLASHFieldContainer
    _handle = None
    
    def __init__(self, filename, data_style='flash_hdf5',
                 storage_filename = None,
                 conversion_override = None):

        if conversion_override is None: conversion_override = {}
        self._conversion_override = conversion_override

        StaticOutput.__init__(self, filename, data_style)
        self.storage_filename = storage_filename

        self.field_info = self._fieldinfo_class()
        # These should be explicitly obtained from the file, but for now that
        # will wait until a reorganization of the source tree and better
        # generalization.
        self.dimensionality = 3
        self.refine_by = 2
        self.parameters["HydroMethod"] = 'flash' # always PPM DE
        self.parameters["Time"] = 1. # default unit is 1...
        self._set_units()
        
    def _set_units(self):
        """
        Generates the conversion to various physical _units based on the parameter file
        """
        self.units = {}
        self.time_units = {}
        if len(self.parameters) == 0:
            self._parse_parameter_file()
        self.conversion_factors = defaultdict(lambda: 1.0)
        if self.cosmological_simulation == 1:
            self._setup_comoving_units()
        else:
            self._setup_nounits_units()
        self.time_units['1'] = 1
        self.units['1'] = 1.0
        self.units['unitary'] = 1.0 / \
            (self.domain_right_edge - self.domain_left_edge).max()
        seconds = 1 #self["Time"]
        self.time_units['years'] = seconds / (365*3600*24.0)
        self.time_units['days']  = seconds / (3600*24.0)
        for p, v in self._conversion_override.items():
            self.conversion_factors[p] = v

    def _setup_comoving_units(self):
        self.conversion_factors['dens'] = (1.0 + self.current_redshift)**3.0
        self.conversion_factors['pres'] = (1.0 + self.current_redshift)**1.0
        self.conversion_factors['eint'] = (1.0 + self.current_redshift)**-2.0
        self.conversion_factors['ener'] = (1.0 + self.current_redshift)**-2.0
        self.conversion_factors['temp'] = (1.0 + self.current_redshift)**-2.0
        self.conversion_factors['velx'] = (1.0 + self.current_redshift)
        self.conversion_factors['vely'] = self.conversion_factors['velx']
        self.conversion_factors['velz'] = self.conversion_factors['velx']
        self.conversion_factors['particle_velx'] = (1.0 + self.current_redshift)
        self.conversion_factors['particle_vely'] = \
            self.conversion_factors['particle_velx']
        self.conversion_factors['particle_velz'] = \
            self.conversion_factors['particle_velx']
        if not self.has_key("TimeUnits"):
            mylog.warning("No time units.  Setting 1.0 = 1 second.")
            self.conversion_factors["Time"] = 1.0
        for unit in mpc_conversion.keys():
            self.units[unit] = mpc_conversion[unit] / mpc_conversion["cm"]

    def _setup_nounits_units(self):
        self.conversion_factors['dens'] = 1.0
        self.conversion_factors['pres'] = 1.0
        self.conversion_factors['eint'] = 1.0
        self.conversion_factors['ener'] = 1.0
        self.conversion_factors['temp'] = 1.0
        self.conversion_factors['velx'] = 1.0
        self.conversion_factors['vely'] = 1.0
        self.conversion_factors['velz'] = 1.0
        self.conversion_factors['particle_velx'] = 1.0
        self.conversion_factors['particle_vely'] = 1.0
        self.conversion_factors['particle_velz'] = 1.0
        z = 0
        mylog.warning("Setting 1.0 in code units to be 1.0 cm")
        if not self.has_key("TimeUnits"):
            mylog.warning("No time units.  Setting 1.0 = 1 second.")
            self.conversion_factors["Time"] = 1.0
        for unit in mpc_conversion.keys():
            self.units[unit] = mpc_conversion[unit] / mpc_conversion["cm"]

    def _find_parameter(self, ptype, pname, scalar = False, handle = None):
        # We're going to implement handle caching eventually
        if handle is None: handle = self._handle
        if handle is None:
            handle = h5py.File(self.parameter_filename, "r")
        nn = "/%s %s" % (ptype,
                {False: "runtime parameters", True: "scalars"}[scalar])
        for tpname, pval in handle[nn][:]:
            if tpname.strip() == pname:
                return pval
        raise KeyError(pname)

    def _parse_parameter_file(self):
        self.unique_identifier = \
            int(os.stat(self.parameter_filename)[stat.ST_CTIME])
        self._handle = h5py.File(self.parameter_filename, "r")
        if "file format version" in self._handle:
            self._flash_version = int(
                self._handle["file format version"][:])
        elif "sim info" in self._handle:
            self._flash_version = int(
                self._handle["sim info"][:]["file format version"])
        else:
            raise RuntimeError("Can't figure out FLASH file version.")
        self.domain_left_edge = na.array(
            [self._find_parameter("real", "%smin" % ax) for ax in 'xyz'])
        self.domain_right_edge = na.array(
            [self._find_parameter("real", "%smax" % ax) for ax in 'xyz'])
        if self._flash_version == 7:
            self.current_time = float(
                self._handle["simulation parameters"][:]["time"])
        else:
            self.current_time = \
                float(self._find_parameter("real", "time", scalar=True))

        try:
            use_cosmo = self._find_parameter("logical", "usecosmology") 
        except KeyError:
            use_cosmo = 0

        if use_cosmo == 1:
            self.cosmological_simulation = 1
            self.current_redshift = self._find_parameter("real", "redshift",
                                        scalar = True)
            self.omega_lambda = self._find_parameter("real", "cosmologicalconstant")
            self.omega_matter = self._find_parameter("real", "omegamatter")
            self.hubble_constant = self._find_parameter("real", "hubbleconstant")
        else:
            self.current_redshift = self.omega_lambda = self.omega_matter = \
                self.hubble_constant = self.cosmological_simulation = 0.0
        self._handle.close()

    @classmethod
    def _is_valid(self, *args, **kwargs):
        try:
            fileh = h5py.File(args[0],'r')
            if "bounding box" in fileh["/"].keys():
                return True
        except:
            pass
        return False


