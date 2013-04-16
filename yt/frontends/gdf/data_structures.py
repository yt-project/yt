"""
Data structures for GDF.

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
import numpy as np
import weakref
from yt.funcs import *
from yt.data_objects.grid_patch import \
           AMRGridPatch
from yt.geometry.grid_geometry_handler import \
           GridGeometryHandler
from yt.data_objects.static_output import \
           StaticOutput
from yt.utilities.lib import \
    get_box_grids_level
from yt.utilities.io_handler import \
    io_registry
from yt.utilities.definitions import \
    mpc_conversion, sec_conversion

from .fields import GDFFieldInfo, KnownGDFFields
from yt.data_objects.field_info_container import \
    FieldInfoContainer, NullFunc
import pdb

def _get_convert(fname):
    def _conv(data):
        return data.convert(fname)
    return _conv

class GDFGrid(AMRGridPatch):
    _id_offset = 0
    def __init__(self, id, hierarchy, level, start, dimensions):
        AMRGridPatch.__init__(self, id, filename = hierarchy.hierarchy_filename,
                              hierarchy = hierarchy)
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
            self.dds = np.array((RE-LE)/self.ActiveDimensions)
        if self.pf.dimensionality < 2: self.dds[1] = 1.0
        if self.pf.dimensionality < 3: self.dds[2] = 1.0
        self.field_data['dx'], self.field_data['dy'], self.field_data['dz'] = self.dds

    @property
    def filename(self):
        return None

class GDFHierarchy(GridGeometryHandler):

    grid = GDFGrid

    def __init__(self, pf, data_style='grid_data_format'):
        self.parameter_file = weakref.proxy(pf)
        self.data_style = data_style
        self.max_level = 10  # FIXME
        # for now, the hierarchy file is the parameter file!
        self.hierarchy_filename = self.parameter_file.parameter_filename
        self.directory = os.path.dirname(self.hierarchy_filename)
#        self._handle = h5py.File(self.hierarchy_filename, 'r')
        self._handle = pf._handle
#        import pudb; pudb.set_trace()
        GridGeometryHandler.__init__(self, pf, data_style)
        print "!!!!"

#        self._handle.close()

    def _initialize_data_storage(self):
        pass

    def _detect_fields(self):
        self.field_list = self._handle['field_types'].keys()

    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        GridGeometryHandler._setup_classes(self, dd)
        self.object_types.sort()

    def _count_grids(self):
        self.num_grids = self._handle['/grid_parent_id'].shape[0]

    def _parse_hierarchy(self):
        f = self._handle
        dxs = []
        self.grids = np.empty(self.num_grids, dtype='object')
        levels = (f['grid_level'][:]).copy()
        glis = (f['grid_left_index'][:]).copy()
        gdims = (f['grid_dimensions'][:]).copy()
        active_dims = ~((np.max(gdims, axis=0) == 1) &
                        (self.parameter_file.domain_dimensions == 1))

        for i in range(levels.shape[0]):
            self.grids[i] = self.grid(i, self, levels[i],
                                      glis[i],
                                      gdims[i])
            self.grids[i]._level_id = levels[i]

            dx = (self.parameter_file.domain_right_edge-
                  self.parameter_file.domain_left_edge)/self.parameter_file.domain_dimensions
            dx[active_dims] = dx[active_dims]/self.parameter_file.refine_by**(levels[i])
            dxs.append(dx)
        dx = np.array(dxs)
        self.grid_left_edge = self.parameter_file.domain_left_edge + dx*glis
        self.grid_dimensions = gdims.astype("int32")
        self.grid_right_edge = self.grid_left_edge + dx*self.grid_dimensions
        self.grid_particle_count = f['grid_particle_count'][:]
        del levels, glis, gdims

    def _populate_grid_objects(self):
        mask = np.empty(self.grids.size, dtype='int32')
        for gi, g in enumerate(self.grids):
            g._prepare_grid()
            g._setup_dx()
        return
        for gi, g in enumerate(self.grids):
            g.Children = self._get_grid_children(g)
            for g1 in g.Children:
                g1.Parent.append(g)
            get_box_grids_level(self.grid_left_edge[gi,:],
                                self.grid_right_edge[gi,:],
                                self.grid_levels[gi],
                                self.grid_left_edge, self.grid_right_edge,
                                self.grid_levels, mask)
            m = mask.astype("bool")
            m[gi] = False
            siblings = self.grids[gi:][m[gi:]]
            if len(siblings) > 0:
                g.OverlappingSiblings = siblings.tolist()
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

class GDFStaticOutput(StaticOutput):
    _hierarchy_class = GDFHierarchy
    _fieldinfo_fallback = GDFFieldInfo
    _fieldinfo_known = KnownGDFFields
    _handle = None

    def __init__(self, filename, data_style='grid_data_format',
                 storage_filename = None):
        if self._handle is not None: return
        self._handle = h5py.File(filename, "r")
        self.storage_filename = storage_filename
        self.filename = filename
        StaticOutput.__init__(self, filename, data_style)

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

        # This should be improved.
        self._handle = h5py.File(self.parameter_filename, "r")
        for field_name in self._handle["/field_types"]:
            current_field = self._handle["/field_types/%s" % field_name]
            try:
                self.units[field_name] = current_field.attrs['field_to_cgs']
            except:
                self.units[field_name] = 1.0
            try:
                current_fields_unit = current_field.attrs['field_units'][0]
            except:
                current_fields_unit = ""
            self._fieldinfo_known.add_field(field_name, function=NullFunc, take_log=False,
                   units=current_fields_unit, projected_units="",
                   convert_function=_get_convert(field_name))
        for p, v in self.units.items():
            self.conversion_factors[p] = v
#        self._handle.close()
#        del self._handle

    def _parse_parameter_file(self):
        self._handle = h5py.File(self.parameter_filename, "r")
        sp = self._handle["/simulation_parameters"].attrs
        self.domain_left_edge = sp["domain_left_edge"][:]
        self.domain_right_edge = sp["domain_right_edge"][:]
        self.domain_dimensions = sp["domain_dimensions"][:]
        refine_by = sp["refine_by"]
        if refine_by is None: refine_by = 2
        self.refine_by = refine_by
        self.dimensionality = sp["dimensionality"]
        self.current_time = sp["current_time"]
        self.unique_identifier = sp["unique_identifier"]
        self.cosmological_simulation = sp["cosmological_simulation"]
        if sp["num_ghost_zones"] != 0: raise RuntimeError
        self.num_ghost_zones = sp["num_ghost_zones"]
        self.field_ordering = sp["field_ordering"]
        self.boundary_conditions = sp["boundary_conditions"][:]
        p = [bnd == 0 for bnd in self.boundary_conditions[::2]]
        self.periodicity = ensure_tuple(p)
        if self.cosmological_simulation:
            self.current_redshift = sp["current_redshift"]
            self.omega_lambda = sp["omega_lambda"]
            self.omega_matter = sp["omega_matter"]
            self.hubble_constant = sp["hubble_constant"]
        else:
            self.current_redshift = self.omega_lambda = self.omega_matter = \
                self.hubble_constant = self.cosmological_simulation = 0.0
        self.parameters['Time'] = 1.0 # Hardcode time conversion for now.
        self.parameters["HydroMethod"] = 0 # Hardcode for now until field staggering is supported.
#        self._handle.close()
#        del self._handle

    @classmethod
    def _is_valid(self, *args, **kwargs):
        try:
            fileh = h5py.File(args[0],'r')
            if "gridded_data_format" in fileh:
                fileh.close()
                return True
            fileh.close()
        except:
            pass
        return False

    def __repr__(self):
        return self.basename.rsplit(".", 1)[0]

    def __del__(self):
        self._handle.close()
