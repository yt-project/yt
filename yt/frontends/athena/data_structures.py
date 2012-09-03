"""
Data structures for Athena.

Author: Samuel W. Skillman <samskillman@gmail.com>
Affiliation: University of Colorado at Boulder
Author: Matthew Turk <matthewturk@gmail.com>
Author: J. S. Oishi <jsoishi@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2008-2011 Samuel W. Skillman, Matthew Turk, J. S. Oishi.  
  All Rights Reserved.

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

from .fields import AthenaFieldInfo, KnownAthenaFields
from yt.data_objects.field_info_container import \
    FieldInfoContainer, NullFunc
import pdb

def _get_convert(fname):
    def _conv(data):
        return data.convert(fname)
    return _conv

class AthenaGrid(AMRGridPatch):
    _id_offset = 0
    def __init__(self, id, hierarchy, level, start, dimensions):
        df = hierarchy.storage_filename
        if 'id0' not in hierarchy.parameter_file.filename:
            gname = hierarchy.parameter_file.filename
        else:
            if id == 0:
                gname = 'id0/' + df + '.vtk'
            else:
                gname = 'id%i/' % id + df[:-5] + '-id%i'%id + df[-5:] + '.vtk'
        AMRGridPatch.__init__(self, id, filename = gname,
                              hierarchy = hierarchy)
        self.filename = gname
        self.Parent = []
        self.Children = []
        self.Level = level
        self.start_index = start.copy()
        self.stop_index = self.start_index + dimensions
        self.ActiveDimensions = dimensions.copy()

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
        self.field_data['dx'], self.field_data['dy'], self.field_data['dz'] = self.dds

def parse_line(line, grid):
    # grid is a dictionary
    splitup = line.strip().split()
    if "vtk" in splitup:
        grid['vtk_version'] = splitup[-1]
    elif "Really" in splitup:
        grid['time'] = splitup[-1]
    elif 'PRIMITIVE' in splitup:
        grid['time'] = float(splitup[4].rstrip(','))
        grid['level'] = int(splitup[6].rstrip(','))
        grid['domain'] = int(splitup[8].rstrip(','))
    elif "DIMENSIONS" in splitup:
        grid['dimensions'] = na.array(splitup[-3:]).astype('int')
    elif "ORIGIN" in splitup:
        grid['left_edge'] = na.array(splitup[-3:]).astype('float64')
    elif "SPACING" in splitup:
        grid['dds'] = na.array(splitup[-3:]).astype('float64')
    elif "CELL_DATA" in splitup:
        grid["ncells"] = int(splitup[-1])
    elif "SCALARS" in splitup:
        field = splitup[1]
        grid['read_field'] = field
        grid['read_type'] = 'scalar'
    elif "VECTORS" in splitup:
        field = splitup[1]
        grid['read_field'] = field
        grid['read_type'] = 'vector'



class AthenaHierarchy(AMRHierarchy):

    grid = AthenaGrid
    _data_style='athena'
    
    def __init__(self, pf, data_style='athena'):
        self.parameter_file = weakref.proxy(pf)
        self.data_style = data_style
        # for now, the hierarchy file is the parameter file!
        self.storage_filename = self.parameter_file.storage_filename
        self.hierarchy_filename = self.parameter_file.filename
        #self.directory = os.path.dirname(self.hierarchy_filename)
        self._fhandle = file(self.hierarchy_filename,'rb')
        AMRHierarchy.__init__(self,pf,data_style)

        self._fhandle.close()

    def _initialize_data_storage(self):
        pass

    def _detect_fields(self):
        field_map = {}
        f = open(self.hierarchy_filename,'rb')
        line = f.readline()
        while line != '':
            splitup = line.strip().split()
            if "DIMENSIONS" in splitup:
                grid_dims = na.array(splitup[-3:]).astype('int')
                line = f.readline()
                continue
            elif "CELL_DATA" in splitup:
                grid_ncells = int(splitup[-1])
                line = f.readline()
                if na.prod(grid_dims) != grid_ncells:
                    grid_dims -= 1
                    grid_dims[grid_dims==0]=1
                if na.prod(grid_dims) != grid_ncells:
                    mylog.error('product of dimensions %i not equal to number of cells %i' % 
                          (na.prod(grid_dims), grid_ncells))
                    raise TypeError
                break
            else:
                del line
                line = f.readline()
        read_table = False
        read_table_offset = f.tell()
        while line != '':
            if len(line) == 0: break
            splitup = line.strip().split()
            if 'SCALARS' in splitup:
                field = splitup[1]
                if not read_table:
                    line = f.readline() # Read the lookup table line
                    read_table = True
                field_map[field] = 'scalar',f.tell() - read_table_offset
                read_table=False

            elif 'VECTORS' in splitup:
                field = splitup[1]
                vfield = field+'_x'
                field_map[vfield] = 'vector',f.tell() - read_table_offset
                vfield = field+'_y'
                field_map[vfield] = 'vector',f.tell() - read_table_offset
                vfield = field+'_z'
                field_map[vfield] = 'vector',f.tell() - read_table_offset
            del line
            line = f.readline()

        f.close()
        del f

        self.field_list = field_map.keys()
        self._field_map = field_map

    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        AMRHierarchy._setup_classes(self, dd)
        self.object_types.sort()

    def _count_grids(self):
        self.num_grids = self.parameter_file.nvtk

    def _parse_hierarchy(self):
        f = open(self.hierarchy_filename,'rb')
        grid = {}
        grid['read_field'] = None
        grid['read_type'] = None
        table_read=False
        line = f.readline()
        while grid['read_field'] is None:
            parse_line(line, grid)
            if "SCALAR" in line.strip().split():
                break
            if "VECTOR" in line.strip().split():
                break
            if 'TABLE' in line.strip().split():
                break
            if len(line) == 0: break
            del line
            line = f.readline()
        f.close()
        del f

        if na.prod(grid['dimensions']) != grid['ncells']:
            grid['dimensions'] -= 1
            grid['dimensions'][grid['dimensions']==0]=1
        if na.prod(grid['dimensions']) != grid['ncells']:
            mylog.error('product of dimensions %i not equal to number of cells %i' % 
                  (na.prod(grid['dimensions']), grid['ncells']))
            raise TypeError

        dxs=[]
        self.grids = na.empty(self.num_grids, dtype='object')
        levels = na.zeros(self.num_grids, dtype='int32')
        single_grid_width = grid['dds']*grid['dimensions']
        grids_per_dim = (self.parameter_file.domain_width/single_grid_width).astype('int32')
        glis = na.empty((self.num_grids,3), dtype='int64')
        for i in range(self.num_grids):
            procz = i/(grids_per_dim[0]*grids_per_dim[1])
            procy = (i - procz*(grids_per_dim[0]*grids_per_dim[1]))/grids_per_dim[0]
            procx = i - procz*(grids_per_dim[0]*grids_per_dim[1]) - procy*grids_per_dim[1]
            glis[i, 0] = procx*grid['dimensions'][0]
            glis[i, 1] = procy*grid['dimensions'][1]
            glis[i, 2] = procz*grid['dimensions'][2]
        gdims = na.ones_like(glis)
        gdims[:] = grid['dimensions']
        for i in range(levels.shape[0]):
            self.grids[i] = self.grid(i, self, levels[i],
                                      glis[i],
                                      gdims[i])
            self.grids[i]._level_id = levels[i]

            dx = (self.parameter_file.domain_right_edge-
                  self.parameter_file.domain_left_edge)/self.parameter_file.domain_dimensions
            dx = dx/self.parameter_file.refine_by**(levels[i])
            dxs.append(grid['dds'])
        dx = na.array(dxs)
        self.grid_left_edge = self.parameter_file.domain_left_edge + dx*glis
        self.grid_dimensions = gdims.astype("int32")
        self.grid_right_edge = self.grid_left_edge + dx*self.grid_dimensions
        self.grid_particle_count = na.zeros([self.num_grids, 1], dtype='int64')
        del levels, glis, gdims

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
        mask = na.zeros(self.num_grids, dtype='bool')
        grids, grid_ind = self.get_box_grids(grid.LeftEdge, grid.RightEdge)
        mask[grid_ind] = True
        return [g for g in self.grids[mask] if g.Level == grid.Level + 1]

class AthenaStaticOutput(StaticOutput):
    _hierarchy_class = AthenaHierarchy
    _fieldinfo_fallback = AthenaFieldInfo
    _fieldinfo_known = KnownAthenaFields
    _data_style = "athena"

    def __init__(self, filename, data_style='athena',
                 storage_filename = None, parameters = {}):
        StaticOutput.__init__(self, filename, data_style)
        self.filename = filename
        self.storage_filename = filename[4:-4]
        self.specified_parameters = parameters

    def _set_units(self):
        """
        Generates the conversion to various physical _units based on the parameter file
        """
        self.units = {}
        self.time_units = {}
        if len(self.parameters) == 0:
            self._parse_parameter_file()
        self.time_units['1'] = 1
        self.units['1'] = 1.0
        self.units['cm'] = 1.0
        self.units['unitary'] = 1.0 / (self.domain_right_edge - self.domain_left_edge).max()
        for unit in mpc_conversion.keys():
            self.units[unit] = 1.0 * mpc_conversion[unit] / mpc_conversion["cm"]
        for unit in sec_conversion.keys():
            self.time_units[unit] = 1.0 / sec_conversion[unit]

        # Here should read through and add fields.

        #default_fields=['density']
        # for field in self.field_list:
        #     self.units[field] = 1.0
        #     self._fieldinfo_known.add_field(field, function=NullFunc, take_log=False,
        #             units="", projected_units="",
        #             convert_function=None)

        # This should be improved.
        # self._handle = h5py.File(self.parameter_filename, "r")
        # for field_name in self._handle["/field_types"]:
        #     current_field = self._handle["/field_types/%s" % field_name]
        #     try:
        #         self.units[field_name] = current_field.attrs['field_to_cgs']
        #     except:
        #         self.units[field_name] = 1.0
        #     try:
        #         current_fields_unit = current_field.attrs['field_units'][0]
        #     except:
        #         current_fields_unit = ""
        #     self._fieldinfo_known.add_field(field_name, function=NullFunc, take_log=False,
        #            units=current_fields_unit, projected_units="", 
        #            convert_function=_get_convert(field_name))

        # self._handle.close()
        # del self._handle

    def _parse_parameter_file(self):
        self._handle = open(self.parameter_filename, "rb")
        # Read the start of a grid to get simulation parameters.
        grid = {}
        grid['read_field'] = None
        line = self._handle.readline()
        while grid['read_field'] is None:
            parse_line(line, grid)
            if "SCALAR" in line.strip().split():
                break
            if "VECTOR" in line.strip().split():
                break
            if 'TABLE' in line.strip().split():
                break
            if len(line) == 0: break
            del line
            line = self._handle.readline()

        self.domain_left_edge = grid['left_edge']
        try:
            self.domain_right_edge = na.array(self.specified_parameters['domain_right_edge'])
        except:
            mylog.info("Please set 'domain_right_edge' in parameters dictionary argument " +
                    "if it is not equal to -domain_left_edge.")
        self.domain_width = self.domain_right_edge-self.domain_left_edge
        self.domain_dimensions = self.domain_width/grid['dds']
        refine_by = None
        if refine_by is None: refine_by = 2
        self.refine_by = refine_by
        self.dimensionality = 3
        self.current_time = grid["time"]
        self.unique_identifier = None
        self.cosmological_simulation = False
        self.num_ghost_zones = 0
        self.field_ordering = 'fortran'
        self.boundary_conditions = [1]*6

        self.nvtk = int(na.product(self.domain_dimensions/(grid['dimensions']-1)))

        # if self.cosmological_simulation:
        #     self.current_redshift = sp["current_redshift"]
        #     self.omega_lambda = sp["omega_lambda"]
        #     self.omega_matter = sp["omega_matter"]
        #     self.hubble_constant = sp["hubble_constant"]
        # else:
        self.current_redshift = self.omega_lambda = self.omega_matter = \
            self.hubble_constant = self.cosmological_simulation = 0.0
        self.parameters['Time'] = self.current_time # Hardcode time conversion for now.
        self.parameters["HydroMethod"] = 0 # Hardcode for now until field staggering is supported.
        self._handle.close()
        del self._handle

    @classmethod
    def _is_valid(self, *args, **kwargs):
        try:
            if 'vtk' in args[0]:
                return True
        except:
            pass
        return False

    def __repr__(self):
        return self.basename.rsplit(".", 1)[0]

