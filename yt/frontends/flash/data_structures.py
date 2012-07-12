"""
FLASH-specific data structures

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: UCSD
Homepage: http://yt-project.org/
License:
  Copyright (C) 2010-2011 Matthew Turk.  All Rights Reserved.

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
    mpc_conversion, sec_conversion
from yt.utilities.io_handler import \
    io_registry

from .fields import FLASHFieldInfo, add_flash_field, KnownFLASHFields
from yt.data_objects.field_info_container import FieldInfoContainer, NullFunc, \
     ValidateDataField

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
    
    def __init__(self,pf,data_style='flash_hdf5'):
        self.data_style = data_style
        self.field_indexes = {}
        self.parameter_file = weakref.proxy(pf)
        # for now, the hierarchy file is the parameter file!
        self.hierarchy_filename = self.parameter_file.parameter_filename
        self.directory = os.path.dirname(self.hierarchy_filename)
        self._handle = pf._handle

        self.float_type = na.float64
        AMRHierarchy.__init__(self,pf,data_style)

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
                "integer", "globalnumblocks", True)
        except KeyError:
            self.num_grids = self._handle["/simulation parameters"][0][0]
        
    def _parse_hierarchy(self):
        f = self._handle # shortcut
        pf = self.parameter_file # shortcut
        
        # Initialize to the domain left / domain right
        ND = self.parameter_file.dimensionality
        DLE = self.parameter_file.domain_left_edge
        DRE = self.parameter_file.domain_right_edge
        for i in range(3):
            self.grid_left_edge[:,i] = DLE[i]
            self.grid_right_edge[:,i] = DRE[i]
        # We only go up to ND for 2D datasets
        if (f["/bounding box"][:,:,0].shape[1] == ND) :
            #FLASH 2/3 2D data
            self.grid_left_edge[:,:ND] = f["/bounding box"][:,:,0]
            self.grid_right_edge[:,:ND] = f["/bounding box"][:,:,1]
        else:
            self.grid_left_edge[:,:] = f["/bounding box"][:,:,0]
            self.grid_right_edge[:,:] = f["/bounding box"][:,:,1]
            

        # Move this to the parameter file
        try:
            nxb = pf.parameters['nxb']
            nyb = pf.parameters['nyb']
            nzb = pf.parameters['nzb']
        except KeyError:
            nxb, nyb, nzb = [int(f["/simulation parameters"]['n%sb' % ax])
                              for ax in 'xyz']
        self.grid_dimensions[:] *= (nxb, nyb, nzb)
        try:
            self.grid_particle_count[:] = f["/localnp"][:][:,None]
        except KeyError:
            self.grid_particle_count[:] = 0.0
        self._particle_indices = na.zeros(self.num_grids + 1, dtype='int64')
        na.add.accumulate(self.grid_particle_count.squeeze(),
                          out=self._particle_indices[1:])
        # This will become redundant, as _prepare_grid will reset it to its
        # current value.  Note that FLASH uses 1-based indexing for refinement
        # levels, but we do not, so we reduce the level by 1.
        self.grid_levels.flat[:] = f["/refine level"][:][:] - 1
        self.grids = na.empty(self.num_grids, dtype='object')
        for i in xrange(self.num_grids):
            self.grids[i] = self.grid(i+1, self, self.grid_levels[i,0])
        

        # This is a possibly slow and verbose fix, and should be re-examined!
        rdx = (self.parameter_file.domain_right_edge -
                self.parameter_file.domain_left_edge)/self.parameter_file.domain_dimensions
        nlevels = self.grid_levels.max()
        dxs = na.zeros((nlevels+1,3),dtype='float64')
        for i in range(nlevels+1):
            dxs[i] = rdx/self.parameter_file.refine_by**i
       
        for i in xrange(self.num_grids):
            dx = dxs[self.grid_levels[i],:]
            self.grid_left_edge[i] = na.rint(self.grid_left_edge[i]/dx)*dx
            self.grid_right_edge[i] = na.rint(self.grid_right_edge[i]/dx)*dx
                        
    def _populate_grid_objects(self):
        # We only handle 3D data, so offset is 7 (nfaces+1)
        
        offset = 7
        ii = na.argsort(self.grid_levels.flat)
        gid = self._handle["/gid"][:]
        first_ind = -(self.parameter_file.refine_by**self.parameter_file.dimensionality)
        for g in self.grids[ii].flat:
            gi = g.id - g._id_offset
            # FLASH uses 1-indexed group info
            g.Children = [self.grids[i - 1] for i in gid[gi,first_ind:] if i > -1]
            for g1 in g.Children:
                g1.Parent = g
            g._prepare_grid()
            g._setup_dx()
        if self.parameter_file.dimensionality < 3:
            DD = (self.parameter_file.domain_right_edge[2] -
                  self.parameter_file.domain_left_edge[2])
            for g in self.grids:
                g.dds[2] = DD
        if self.parameter_file.dimensionality < 2:
            DD = (self.parameter_file.domain_right_edge[1] -
                  self.parameter_file.domain_left_edge[1])
            for g in self.grids:
                g.dds[1] = DD
        self.max_level = self.grid_levels.max()

    def _setup_derived_fields(self):
        AMRHierarchy._setup_derived_fields(self)
        [self.parameter_file.conversion_factors[field] 
         for field in self.field_list]
        for field in self.field_list:
            if field not in self.derived_field_list:
                self.derived_field_list.append(field)
            if (field not in KnownFLASHFields and
                field.startswith("particle")) :
                self.parameter_file.field_info.add_field(field,
                                                         function=NullFunc,
                                                         take_log=False,
                                                         validators = [ValidateDataField(field)],
                                                         particle_type=True)
                
    def _setup_data_io(self):
        self.io = io_registry[self.data_style](self.parameter_file)

class FLASHStaticOutput(StaticOutput):
    _hierarchy_class = FLASHHierarchy
    _fieldinfo_fallback = FLASHFieldInfo
    _fieldinfo_known = KnownFLASHFields
    _handle = None
    
    def __init__(self, filename, data_style='flash_hdf5',
                 storage_filename = None,
                 conversion_override = None):

        self._handle = h5py.File(filename, "r")
        if conversion_override is None: conversion_override = {}
        self._conversion_override = conversion_override

        StaticOutput.__init__(self, filename, data_style)
        self.storage_filename = storage_filename

        # These should be explicitly obtained from the file, but for now that
        # will wait until a reorganization of the source tree and better
        # generalization.
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
        if "EOSType" not in self.parameters:
            self.parameters["EOSType"] = -1
        if self.cosmological_simulation == 1:
            self._setup_comoving_units()
        if "pc_unitsbase" in self.parameters:
            if self.parameters["pc_unitsbase"] == "CGS":
                self._setup_cgs_units()
        else:
            self._setup_nounits_units()
        self.time_units['1'] = 1
        self.units['1'] = 1.0
        self.units['unitary'] = 1.0 / \
            (self.domain_right_edge - self.domain_left_edge).max()
        for unit in sec_conversion.keys():
            self.time_units[unit] = 1.0 / sec_conversion[unit]

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

    def _setup_cgs_units(self):
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
        mylog.warning("Setting 1.0 in code units to be 1.0 cm")
        if not self.has_key("TimeUnits"):
            mylog.warning("No time units.  Setting 1.0 = 1 second.")
            self.conversion_factors["Time"] = 1.0
        for unit in mpc_conversion.keys():
            self.units[unit] = mpc_conversion[unit] / mpc_conversion["cm"]

    def _find_parameter(self, ptype, pname, scalar = False):
        nn = "/%s %s" % (ptype,
                {False: "runtime parameters", True: "scalars"}[scalar])
        if nn not in self._handle: raise KeyError(nn)
        for tpname, pval in zip(self._handle[nn][:,'name'],
                                self._handle[nn][:,'value']):
            if tpname.strip() == pname:
                if ptype == "string" :
                    return pval.strip()
                else :
                    return pval
        raise KeyError(pname)

    def _parse_parameter_file(self):
        self.unique_identifier = \
            int(os.stat(self.parameter_filename)[stat.ST_CTIME])
        if "file format version" in self._handle:
            self._flash_version = int(
                self._handle["file format version"][:])
        elif "sim info" in self._handle:
            self._flash_version = int(
                self._handle["sim info"][:]["file format version"])
        else:
            raise RuntimeError("Can't figure out FLASH file version.")
        # First we load all of the parameters
        hns = ["simulation parameters"]
        # note the ordering here is important: runtime parameters should
        # ovewrite scalars with the same name.
        for ptype in ['scalars', 'runtime parameters']:
            for vtype in ['integer', 'real', 'logical', 'string']:
                hns.append("%s %s" % (vtype, ptype))
        if self._flash_version > 7:
            for hn in hns:
                if hn not in self._handle:
                    continue
                for varname, val in zip(self._handle[hn][:,'name'],
                                        self._handle[hn][:,'value']):
                    vn = varname.strip()
                    if hn.startswith("string") :
                        pval = val.strip()
                    else :
                        pval = val
                    if vn in self.parameters and self.parameters[vn] != pval:
                        mylog.warning("{0} {1} overwrites a simulation scalar of the same name".format(hn[:-1],vn)) 
                    self.parameters[vn] = pval
        if self._flash_version == 7:
            for hn in hns:
                if hn not in self._handle:
                    continue
                if hn is 'simulation parameters':
                    zipover = zip(self._handle[hn].dtype.names,self._handle[hn][0])
                else:
                    zipover = zip(self._handle[hn][:,'name'],self._handle[hn][:,'value'])
                for varname, val in zipover:
                    vn = varname.strip()
                    if hn.startswith("string") :
                        pval = val.strip()
                    else :
                        pval = val
                    if vn in self.parameters and self.parameters[vn] != pval:
                        mylog.warning("{0} {1} overwrites a simulation scalar of the same name".format(hn[:-1],vn))
                    self.parameters[vn] = pval
        self.domain_left_edge = na.array(
            [self.parameters["%smin" % ax] for ax in 'xyz']).astype("float64")
        self.domain_right_edge = na.array(
            [self.parameters["%smax" % ax] for ax in 'xyz']).astype("float64")
        self.min_level = self.parameters["lrefine_min"] -1

        # Determine domain dimensions
        try:
            nxb = self.parameters["nxb"]
            nyb = self.parameters["nyb"]
            nzb = self.parameters["nzb"]
        except KeyError:
            nxb, nyb, nzb = [int(self._handle["/simulation parameters"]['n%sb' % ax])
                              for ax in 'xyz'] # FLASH2 only!
        try:
            dimensionality = self.parameters["dimensionality"]
        except KeyError:
            dimensionality = 3
            if nzb == 1: dimensionality = 2
            if nyb == 1: dimensionality = 1
            if dimensionality < 3:
                mylog.warning("Guessing dimensionality as %s", dimensionality)

        nblockx = self.parameters["nblockx"]
        nblocky = self.parameters["nblocky"]
        nblockz = self.parameters["nblockz"]
        self.dimensionality = dimensionality
        self.domain_dimensions = \
            na.array([nblockx*nxb,nblocky*nyb,nblockz*nzb])
        try:
            self.parameters["Gamma"] = self.parameters["gamma"]
        except:
            mylog.warning("Cannot find Gamma")
            pass

        self.current_time = self.parameters["time"]

        try: 
            self.parameters["usecosmology"]
            self.cosmological_simulation = 1
            self.current_redshift = self.parameters['redshift']
            self.omega_lambda = self.parameters['cosmologicalconstant']
            self.omega_matter = self.parameters['omegamatter']
            self.hubble_constant = self.parameters['hubbleconstant']
        except:
            self.current_redshift = self.omega_lambda = self.omega_matter = \
                self.hubble_constant = self.cosmological_simulation = 0.0

    def __del__(self):
        self._handle.close()

    @classmethod
    def _is_valid(self, *args, **kwargs):
        try:
            fileh = h5py.File(args[0],'r')
            if "bounding box" in fileh["/"].keys():
                fileh.close()
                return True
            fileh.close()
        except:
            pass
        return False


