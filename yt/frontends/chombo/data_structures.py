"""
Data structures for Chombo.

Author: Matthew Turk <matthewturk@gmail.com>
Author: J. S. Oishi <jsoishi@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2008-2010 Matthew Turk, J. S. Oishi.  All Rights Reserved.

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
import numpy as na

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
     mpc_conversion

from .fields import ChomboFieldContainer

class ChomboGrid(AMRGridPatch):
    _id_offset = 0
    __slots__ = ["_level_id", "stop_index"]
    def __init__(self, id, hierarchy, level, start, stop):
        AMRGridPatch.__init__(self, id, filename = hierarchy.hierarchy_filename,
                              hierarchy = hierarchy)
        self.Parent = []
        self.Children = []
        self.Level = level
        self.start_index = start.copy()#.transpose()
        self.stop_index = stop.copy()#.transpose()
        self.ActiveDimensions = stop - start + 1

    def _setup_dx(self):
        # So first we figure out what the index is.  We don't assume
        # that dx=dy=dz , at least here.  We probably do elsewhere.
        id = self.id - self._id_offset
        if len(self.Parent) > 0:
            self.dds = self.Parent[0].dds / self.pf.refine_by
        else:
            LE, RE = self.hierarchy.grid_left_edge[id,:], \
                     self.hierarchy.grid_right_edge[id,:]
            self.dds = na.array((RE-LE)/self.ActiveDimensions)
        if self.pf.dimensionality < 2: self.dds[1] = 1.0
        if self.pf.dimensionality < 3: self.dds[2] = 1.0
        self.data['dx'], self.data['dy'], self.data['dz'] = self.dds

class ChomboHierarchy(AMRHierarchy):

    grid = ChomboGrid
    
    def __init__(self,pf,data_style='chombo_hdf5'):
        self.domain_left_edge = pf.domain_left_edge # need these to determine absolute grid locations
        self.domain_right_edge = pf.domain_right_edge # need these to determine absolute grid locations
        self.data_style = data_style
        self.field_info = ChomboFieldContainer()
        self.field_indexes = {}
        self.parameter_file = weakref.proxy(pf)
        # for now, the hierarchy file is the parameter file!
        self.hierarchy_filename = self.parameter_file.parameter_filename
        self.hierarchy = os.path.abspath(self.hierarchy_filename)
        self.directory = os.path.dirname(self.hierarchy_filename)
        self._fhandle = h5py.File(self.hierarchy_filename)

        self.float_type = self._fhandle['/level_0']['data:datatype=0'].dtype.name
        self._levels = self._fhandle.listnames()[1:]
        AMRHierarchy.__init__(self,pf,data_style)

        self._fhandle.close()

    def _initialize_data_storage(self):
        pass

    def _detect_fields(self):
        ncomp = int(self._fhandle['/'].attrs['num_components'])
        self.field_list = [c[1] for c in self._fhandle['/'].attrs.listitems()[-ncomp:]]
    
    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        AMRHierarchy._setup_classes(self, dd)
        self.object_types.sort()

    def _count_grids(self):
        self.num_grids = 0
        for lev in self._levels:
            self.num_grids += self._fhandle[lev]['Processors'].len()
        
    def _parse_hierarchy(self):
        f = self._fhandle # shortcut
        
        # this relies on the first Group in the H5 file being
        # 'Chombo_global'
        levels = f.listnames()[1:]
        self.grids = []
        i = 0
        for lev in levels:
            level_number = int(re.match('level_(\d+)',lev).groups()[0])
            boxes = f[lev]['boxes'].value
            dx = f[lev].attrs['dx']
            for level_id, box in enumerate(boxes):
                si = na.array([box['lo_%s' % ax] for ax in 'ijk'])
                ei = na.array([box['hi_%s' % ax] for ax in 'ijk'])
                pg = self.grid(len(self.grids),self,level=level_number,
                               start = si, stop = ei)
                self.grids.append(pg)
                self.grids[-1]._level_id = level_id
                self.grid_left_edge[i] = dx*si.astype(self.float_type) + self.domain_left_edge
                self.grid_right_edge[i] = dx*(ei.astype(self.float_type)+1) + self.domain_left_edge
                self.grid_particle_count[i] = 0
                self.grid_dimensions[i] = ei - si + 1
                i += 1
        self.grids = na.array(self.grids, dtype='object')

    def _populate_grid_objects(self):
        for g in self.grids:
            g._prepare_grid()
            g._setup_dx()

        for g in self.grids:
            g.Children = self._get_grid_children(g)
            for g1 in g.Children:
                g1.Parent.append(g)
        self.max_level = self.grid_levels.max()

    def _setup_unknown_fields(self):
        pass

    def _setup_derived_fields(self):
        self.derived_field_list = []

    def _get_grid_children(self, grid):
        mask = na.zeros(self.num_grids, dtype='bool')
        grids, grid_ind = self.get_box_grids(grid.LeftEdge, grid.RightEdge)
        mask[grid_ind] = True
        return [g for g in self.grids[mask] if g.Level == grid.Level + 1]

class ChomboStaticOutput(StaticOutput):
    _hierarchy_class = ChomboHierarchy
    _fieldinfo_class = ChomboFieldContainer
    
    def __init__(self, filename, data_style='chombo_hdf5',
                 storage_filename = None, ini_filename = None):
        # hardcoded for now 
        self.current_time = 0.0
        self.ini_filename = ini_filename
        StaticOutput.__init__(self,filename,data_style)
        self.storage_filename = storage_filename
        self.field_info = self._fieldinfo_class()
        
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
        self.time_units['years'] = seconds / (365*3600*24.0)
        self.time_units['days']  = seconds / (3600*24.0)
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
        Check to see whether a 'pluto.ini' or 'orion2.ini' file
        exists in the plot file directory. If one does, attempt to parse it.
        Otherwise, assume the left edge starts at 0 and get the right edge
        from the hdf5 file.
        """
        if os.path.isfile('pluto.ini'):
            self._parse_pluto_file('pluto.ini')
        elif os.path.isfile('orion2.ini'):
            self._parse_pluto_file('orion2.ini')
        else:
            self.unique_identifier = \
                                   int(os.stat(self.parameter_filename)[ST_CTIME])
            self.domain_left_edge = na.array([0.,0.,0.])
            self.domain_right_edge = self.__calc_right_edge()
            self.dimensionality = 3
            self.refine_by = 2

    def _parse_pluto_file(self, ini_filename):
        """
        Reads in an inputs file in the 'pluto.ini' format. Probably not
        especially robust at the moment.
        """
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
        self.domain_left_edge = na.array([low1,low2,low3])
        self.domain_right_edge = na.array([high1,high2,high3])
        self.domain_dimensions = na.array([N1,N2,N3])
        self.refine_by = self.parameters["RefineBy"]
            
    def __calc_right_edge(self):
        fileh = h5py.File(self.parameter_filename,'r')
        dx0 = fileh['/level_0'].attrs['dx']
        RE = dx0*((na.array(fileh['/level_0'].attrs['prob_domain']))[3:] + 1)
        fileh.close()
        return RE
                   
    @classmethod
    def _is_valid(self, *args, **kwargs):
        try:
            fileh = h5py.File(args[0],'r')
            if (fileh.listnames())[0] == 'Chombo_global':
                return True
        except:
            pass
        return False


