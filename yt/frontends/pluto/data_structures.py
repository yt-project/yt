"""
Data structures for Pluto.

Author: Matthew Turk <matthewturk@gmail.com>
Author: J. S. Oishi <jsoishi@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2008-2011 Matthew Turk, J. S. Oishi.  All Rights Reserved.

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
import re
import os
import weakref
import numpy as np

from collections import \
     defaultdict
from string import \
     strip, \
     rstrip
from stat import \
     ST_CTIME

from .definitions import \
     pluto2enzoDict, \
     yt2plutoFieldsDict, \
     parameterDict \

from yt.funcs import *
from yt.data_objects.grid_patch import \
     AMRGridPatch
from yt.data_objects.hierarchy import \
     AMRHierarchy
from yt.data_objects.static_output import \
     StaticOutput
from yt.utilities.definitions import \
     mpc_conversion, sec_conversion
from yt.utilities.parallel_tools.parallel_analysis_interface import \
     parallel_root_only
from yt.utilities.io_handler import \
    io_registry

from yt.data_objects.field_info_container import \
    FieldInfoContainer, NullFunc
from .fields import PlutoFieldInfo, KnownPlutoFields

class PlutoGrid(AMRGridPatch):
    _id_offset = 0
    __slots__ = ["_level_id", "stop_index"]
    def __init__(self, id, hierarchy, level, start, stop):
        AMRGridPatch.__init__(self, id, filename = hierarchy.hierarchy_filename,
                              hierarchy = hierarchy)
        self.Parent = []
        self.Children = []
        self.Level = level
        self.ActiveDimensions = stop - start + 1

    def get_global_startindex(self):
        """
        Return the integer starting index for each dimension at the current
        level.

        """
        if self.start_index != None:
            return self.start_index
        if self.Parent == []:
            iLE = self.LeftEdge - self.pf.domain_left_edge
            start_index = iLE / self.dds
            return np.rint(start_index).astype('int64').ravel()
        pdx = self.Parent[0].dds
        start_index = (self.Parent[0].get_global_startindex()) + \
            np.rint((self.LeftEdge - self.Parent[0].LeftEdge)/pdx)
        self.start_index = (start_index*self.pf.refine_by).astype('int64').ravel()
        return self.start_index

    def _setup_dx(self):
        # has already been read in and stored in hierarchy
        self.dds = self.hierarchy.dds_list[self.Level]
        self.field_data['dx'], self.field_data['dy'], self.field_data['dz'] = self.dds

class PlutoHierarchy(AMRHierarchy):

    grid = PlutoGrid

    def __init__(self,pf,data_style='chombo_hdf5'):
        self.domain_left_edge = pf.domain_left_edge
        self.domain_right_edge = pf.domain_right_edge
        self.data_style = data_style
        self.field_indexes = {}
        self.parameter_file = weakref.proxy(pf)
        # for now, the hierarchy file is the parameter file!
        self.hierarchy_filename = os.path.abspath(
            self.parameter_file.parameter_filename)
        self.directory = pf.fullpath
        self._handle = pf._handle

        self.float_type = self._handle['/level_0']['data:datatype=0'].dtype.name
        self._levels = self._handle.keys()[2:]
        AMRHierarchy.__init__(self,pf,data_style)

    def _detect_fields(self):
        ncomp = int(self._handle['/'].attrs['num_components'])
        self.field_list = [c[1] for c in self._handle['/'].attrs.items()[-ncomp:]]
          
    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        AMRHierarchy._setup_classes(self, dd)
        self.object_types.sort()

    def _count_grids(self):
        self.num_grids = 0
        for lev in self._levels:
            self.num_grids += self._handle[lev]['Processors'].len()

    def _parse_hierarchy(self):
        f = self._handle # shortcut

        # this relies on the first Group in the H5 file being
        # 'Chombo_global' and the second 'Expressions'
        levels = f.keys()[2:]
        grids = []
        self.dds_list = []
        i = 0
        for lev in levels:
            level_number = int(re.match('level_(\d+)',lev).groups()[0])
            boxes = f[lev]['boxes'].value
            dx = f[lev].attrs['dx']
            self.dds_list.append(dx * np.ones(3))
            for level_id, box in enumerate(boxes):
                si = np.array([box['lo_%s' % ax] for ax in 'ijk'])
                ei = np.array([box['hi_%s' % ax] for ax in 'ijk'])
                pg = self.grid(len(grids),self,level=level_number,
                               start = si, stop = ei)
                grids.append(pg)
                grids[-1]._level_id = level_id
                self.grid_left_edge[i] = dx*si.astype(self.float_type) + self.domain_left_edge
                self.grid_right_edge[i] = dx*(ei.astype(self.float_type)+1) + self.domain_left_edge
                self.grid_particle_count[i] = 0
                self.grid_dimensions[i] = ei - si + 1
                i += 1
        self.grids = np.empty(len(grids), dtype='object')
        for gi, g in enumerate(grids): self.grids[gi] = g
#        self.grids = np.array(self.grids, dtype='object')

    def _populate_grid_objects(self):
        for g in self.grids:
            g._prepare_grid()
            g._setup_dx()

        for g in self.grids:
            g.Children = self._get_grid_children(g)
            for g1 in g.Children:
                g1.Parent.append(g)
        self.max_level = self.grid_levels.max()

    def _setup_derived_fields(self):
        self.derived_field_list = []

    def _get_grid_children(self, grid):
        mask = np.zeros(self.num_grids, dtype='bool')
        grids, grid_ind = self.get_box_grids(grid.LeftEdge, grid.RightEdge)
        mask[grid_ind] = True
        return [g for g in self.grids[mask] if g.Level == grid.Level + 1]

    def _setup_data_io(self):
        self.io = io_registry[self.data_style](self.parameter_file)

class PlutoStaticOutput(StaticOutput):
    _hierarchy_class = PlutoHierarchy
    _fieldinfo_fallback = PlutoFieldInfo
    _fieldinfo_known = KnownPlutoFields

    def __init__(self, filename, data_style='chombo_hdf5',
                 storage_filename = None, ini_filename = None):
        self._handle = h5py.File(filename,'r')
        self.current_time = self._handle.attrs['time']
        self.ini_filename = ini_filename
        self.fullplotdir = os.path.abspath(filename)
        StaticOutput.__init__(self,filename,data_style)
        self.storage_filename = storage_filename
        self.cosmological_simulation = False

        # These are parameters that I very much wish to get rid of.
        self.parameters["HydroMethod"] = 'chombo' # always PPM DE
        self.parameters["DualEnergyFormalism"] = 0 
        self.parameters["EOSType"] = -1 # default

    def __del__(self):
        self._handle.close()

    def _set_units(self):
        """
        Generates the conversion to various physical _units based on the parameter file
        """
        self.units = {}
        self.time_units = {}
        if len(self.parameters) == 0:
            self._parse_parameter_file()
        self._setup_nounits_units()
        self.conversion_factors = defaultdict(lambda: 1.0)
        self.time_units['1'] = 1
        self.units['1'] = 1.0
        self.units['unitary'] = 1.0 / (self.domain_right_edge - self.domain_left_edge).max()
        seconds = 1 #self["Time"]
        for unit in sec_conversion.keys():
            self.time_units[unit] = seconds / sec_conversion[unit]
        for key in yt2plutoFieldsDict:
            self.conversion_factors[key] = 1.0

    def _setup_nounits_units(self):
        z = 0
        mylog.warning("Setting 1.0 in code units to be 1.0 cm")
        if not self.has_key("TimeUnits"):
            mylog.warning("No time units.  Setting 1.0 = 1 second.")
            self.conversion_factors["Time"] = 1.0
        for unit in mpc_conversion.keys():
            self.units[unit] = mpc_conversion[unit] / mpc_conversion["cm"]


    def _localize(self, f, default):
        if f is None:
            return os.path.join(self.directory, default)
        return f

    def _parse_parameter_file(self):
        """
        Reads in an inputs file in the 'pluto.ini' format. Probably not
        especially robust at the moment.
        """

        ini_filename = 'pluto.ini'
        self.fullplotdir = os.path.abspath(self.parameter_filename)
        self.ini_filename = self._localize( \
            self.ini_filename, ini_filename)
        self.unique_identifier = \
                               int(os.stat(self.parameter_filename)[ST_CTIME])
        lines = open(self.ini_filename).readlines()
        # read the file line by line, storing important parameters
        for lineI, line in enumerate(lines):
            try:
                param, sep, vals = map(rstrip,line.partition(' '))
            except ValueError:
                mylog.error("ValueError: '%s'", line)
            if pluto2enzoDict.has_key(param):
                paramName = pluto2enzoDict[param]
                t = map(parameterDict[paramName], vals.split())
                if len(t) == 1:
                    self.parameters[paramName] = t[0]
                else:
                    if paramName == "RefineBy":
                        self.parameters[paramName] = t[0]
                    else:
                        self.parameters[paramName] = t

            # assumes 3D for now
            elif param.startswith("X1-grid"):
                t = vals.split()
                low1 = float(t[1])
                high1 = float(t[4])
                N1 = int(t[2])
            elif param.startswith("X2-grid"):
                t = vals.split()
                low2 = float(t[1])
                high2 = float(t[4])
                N2 = int(t[2])
            elif param.startswith("X3-grid"):
                t = vals.split()
                low3 = float(t[1])
                high3 = float(t[4])
                N3 = int(t[2])
            
        self.dimensionality = 3
        self.domain_left_edge = np.array([low1,low2,low3])
        self.domain_right_edge = np.array([high1,high2,high3])
        self.domain_dimensions = np.array([N1,N2,N3])
        self.refine_by = self.parameters["RefineBy"]
            
    @classmethod
    def _is_valid(self, *args, **kwargs):
        return os.path.isfile('pluto.ini')

    @parallel_root_only
    def print_key_parameters(self):
        for a in ["current_time", "domain_dimensions", "domain_left_edge",
                  "domain_right_edge"]:
            if not hasattr(self, a):
                mylog.error("Missing %s in parameter file definition!", a)
                continue
            v = getattr(self, a)
            mylog.info("Parameters: %-25s = %s", a, v)
