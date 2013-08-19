"""
A means of extracting a subset of the hierarchy

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2008-2011 Matthew Turk.  All Rights Reserved.

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

import h5py, os.path
import numpy as np

from yt.funcs import *
from yt.data_objects.data_containers import YTFieldData
from yt.data_objects.grid_patch import \
    AMRGridPatch
from yt.data_objects.static_output import \
    StaticOutput
from yt.data_objects.hierarchy import \
    AMRHierarchy

class DummyHierarchy(object):
    pass

class ConstructedRootGrid(AMRGridPatch):
    __slots__ = ['base_grid', 'id', 'base_pf']
    _id_offset = 1
    def __init__(self, base_pf, pf, hierarchy, level, left_edge, right_edge):
        """
        This is a fake root grid, constructed by creating a
        :class:`yt.data_objects.api.CoveringGridBase` at a given *level* between
        *left_edge* and *right_edge*.
        """
        self.pf = pf
        self.base_pf = base_pf
        self.field_parameters = {}
        self.NumberOfParticles = 0
        self.id = 1
        self.hierarchy = hierarchy
        self._child_mask = self._child_indices = self._child_index_mask = None
        self.Level = level
        self.LeftEdge = left_edge
        self.RightEdge = right_edge
        self.start_index = np.min([grid.get_global_startindex() for grid in
                             base_pf.h.select_grids(level)], axis=0).astype('int64')
        self.dds = base_pf.h.select_grids(level)[0].dds.copy()
        dims = (self.RightEdge-self.LeftEdge)/self.dds
        self.ActiveDimensions = dims
        print "Constructing base grid of size %s" % (self.ActiveDimensions)
        self.base_grid = base_pf.h.smoothed_covering_grid(level, self.LeftEdge,
                        self.RightEdge, dims=dims)
        self.base_grid.Level = self.base_grid.level
        self.field_data = YTFieldData()
        #self._calculate_child_masks()
        self.Parent = None
        self.Children = []

    def get_vertex_centered_data(self, field, smoothed=True):
        vc = self.base_pf.h.smoothed_covering_grid(self.base_grid.Level,
                self.base_grid.LeftEdge - self.base_grid.dds*0.5,
                self.base_grid.RightEdge + self.base_grid.dds*0.5,
                dims = self.ActiveDimensions + 1)
        return vc[field]

class AMRExtractedGridProxy(AMRGridPatch):
    __slots__ = ['base_grid']
    _id_offset = 1
    def __init__(self, grid_id, base_grid, hierarchy):
        # We make a little birdhouse in our soul for the base_grid
        # (they're the only bee in our bonnet!)
        self.base_grid = base_grid
        AMRGridPatch.__init__(self, grid_id, filename = None, hierarchy=hierarchy)
        self.Parent = None
        self.Children = []
        self.Level = -1

    def get_vertex_centered_data(self, *args, **kwargs):
        return self.base_grid.get_vertex_centered_data(*args, **kwargs)

class OldExtractedHierarchy(object):

    def __init__(self, pf, min_level, max_level = -1, offset = None,
                 always_copy=False):
        """
        This is a class that extracts a hierarchy from another hierarchy,
        filling in regions as necessary.  It accepts a parameter file (*pf*), a
        *min_level*, a *max_level*, and alternately an *offset*.  This class is
        typically or exclusively used to extract for the purposes of visualization.
        """
        self.pf = pf
        self.always_copy = always_copy
        self.min_level = min_level
        self.int_offset = np.min([grid.get_global_startindex() for grid in
                             pf.h.select_grids(min_level)], axis=0).astype('float64')
        min_left = np.min([grid.LeftEdge for grid in
                           pf.h.select_grids(min_level)], axis=0).astype('float64')
        max_right = np.max([grid.RightEdge for grid in 
                                   pf.h.select_grids(min_level)], axis=0).astype('float64')
        if offset is None: offset = (max_right + min_left)/2.0
        self.left_edge_offset = offset
        self.mult_factor = 2**min_level
        self.min_left_edge = self._convert_coords(min_left)
        self.max_right_edge = self._convert_coords(max_right)
        if max_level == -1: max_level = pf.h.max_level
        self.max_level = min(max_level, pf.h.max_level)
        self.final_level = self.max_level - self.min_level
        if len(self.pf.h.select_grids(self.min_level)) > 0:
            self._base_grid = ConstructedRootGrid(self.pf, self.min_level,
                               min_left, max_right)
        else: self._base_grid = None
        
    def select_level(self, level):
        if level == 0 and self._base_grid is not None:
            return [self._base_grid]
        return self.pf.h.select_grids(self.min_level + level)

    def export_output(self, afile, n, field):
        # I prefer dict access, but tables doesn't.
        # But h5py does!
        time_node = afile.create_group("/time-%s" % n)
        time_node.attrs['time'] = self.pf.current_time
        time_node.attrs['numLevels'] = self.pf.h.max_level+1-self.min_level
        # Can take a while, so let's get a progressbar
        self._export_all_levels(afile, time_node, field)

    def _export_all_levels(self, afile, time_node, field):
        pbar = yt.funcs.get_pbar("Exporting levels", self.final_level+1)
        for i,grid_set in enumerate(self.get_levels()):
            pbar.update(i)
            self.export_level(afile, time_node, i, field, grid_set)
        pbar.finish()

    def export_level(self, afile, time_node, level, field, grids = None):
        level_node = afile.create_group("%s/level-%s" % (time_node,level))
        # Grid objects on this level...
        if grids is None: grids = self.pf.h.select_grids(level+self.min_level)
        level_node.attrs['delta'] = grids[0].dds*self.mult_factor
        level_node.attrs['relativeRefinementFactor'] = np.array([2]*3, dtype='int32')
        level_node.attrs['numGrids'] = len(grids)
        for i,g in enumerate(grids):
            self.export_grid(afile, level_node, g, i, field)

    def _convert_grid(self, grid):
        int_origin = (grid.get_global_startindex() \
                    - self.int_offset*2**(grid.Level-self.min_level)).astype('int64')
        level_int_origin = (grid.LeftEdge - self.left_edge_offset)/grid.dds
        origin = self._convert_coords(grid.LeftEdge)
        dds = grid.dds * self.mult_factor
        return int_origin, level_int_origin, origin, dds

    def export_grid(self, afile, level_node, grid, i, field):
        grid_node = afile.create_group("%s/grid-%s" % (level_node,i))
        int_origin, lint, origin, dds = self._convert_grid(grid)
        grid_node.attrs['integerOrigin'] = int_origin
        grid_node.attrs['origin'] = origin
        grid_node.attrs['ghostzoneFlags'] = np.zeros(6, dtype='int32')
        grid_node.attrs['numGhostzones'] = np.zeros(3, dtype='int32')
        grid_node.attrs['dims'] = grid.ActiveDimensions[::-1].astype('int32')
        if not self.always_copy and self.pf.h.data_style == 6 \
           and field in self.pf.h.field_list:
            if grid.hierarchy.data_style == -1: # constructed grid
                # if we can get conversion in amira we won't need to do this
                ff = grid[field].astype('float32')
                ff /= self.pf.conversion_factors.get(field, 1.0)
                afile.create_dataset("%s/grid-data" % grid_node, data=ff.swapaxes(0,2))
            else:
                tfn = os.path.abspath(afile.filename)
                gfn = os.path.abspath(grid.filename)
                fpn = os.path.commonprefix([tfn, grid.filename])
                fn = grid.filename[len(os.path.commonprefix([tfn, grid.filename])):]
                grid_node.attrs['referenceFileName'] = fn
                grid_node.attrs['referenceDataPath'] = \
                    "/Grid%08i/%s" % (grid.id, field)
        else:
            # Export our array
            afile.create_dataset("%s/grid-data" % grid_node, 
                                 data = grid[field].astype('float32').swapaxes(0,2))

    def _convert_coords(self, val):
        return (val - self.left_edge_offset)*self.mult_factor

class ExtractedHierarchy(AMRHierarchy):

    grid = AMRExtractedGridProxy

    def __init__(self, pf, data_style):
        # First we set up our translation between original and extracted
        self.data_style = data_style
        self.min_level = pf.min_level
        self.int_offset = np.min([grid.get_global_startindex() for grid in
                           pf.base_pf.h.select_grids(pf.min_level)], axis=0).astype('float64')
        min_left = np.min([grid.LeftEdge for grid in
                           pf.base_pf.h.select_grids(pf.min_level)], axis=0).astype('float64')
        max_right = np.max([grid.RightEdge for grid in 
                           pf.base_pf.h.select_grids(pf.min_level)], axis=0).astype('float64')
        level_dx = pf.base_pf.h.select_grids(pf.min_level)[0].dds[0]
        dims = ((max_right-min_left)/level_dx)
        max_right += (dims.max() - dims) * level_dx
        offset = pf.offset
        if offset is None: offset = min_left
        self.left_edge_offset = offset
        pf.offset = offset
        self.mult_factor = 2**pf.min_level
        self.min_left_edge = self._convert_coords(min_left)
        self.max_right_edge = self._convert_coords(max_right)
        self.min_left, self.max_right = min_left, max_right
        max_level = pf.max_level
        if max_level == -1: max_level = pf.base_pf.h.max_level
        self.max_level = min(max_level, pf.base_pf.h.max_level)
        self.final_level = self.max_level - self.min_level

        # Now we utilize the existing machinery for generating the appropriate
        # arrays of grids, etc etc.
        self.base_pf = pf.base_pf
        AMRHierarchy.__init__(self, pf, data_style)

        # Now a few cleanups
        self.pf.override["DomainRightEdge"] = self.max_right_edge
        self.pf.override["DomainLeftEdge"] = self.min_left_edge
        for u,v in self.base_pf.units.items():
            self.pf.override[u] = v / self.mult_factor
        self.pf.override['unitary'] = 1.0 / (self.pf.domain_right_edge -
                                             self.pf.domain_left_edge).max()

    def _count_grids(self):
        self.num_grids = 1 + sum( ( # 1 is the base grid
            len(self.base_pf.h.select_grids(level)) 
                for level in range(self.min_level+1, self.max_level)) )

    def _parse_hierarchy(self):
        # Here we need to set up the grid info, which for the Enzo hierarchy
        # is done like:
        # self.grid_dimensions.flat[:] = ei
        # self.grid_dimensions -= np.array(si, self.float_type)
        # self.grid_dimensions += 1
        # self.grid_left_edge.flat[:] = LE
        # self.grid_right_edge.flat[:] = RE
        # self.grid_particle_count.flat[:] = np
        # self.grids = np.array(self.grids, dtype='object')
        #
        # For now, we make the presupposition that all of our grids are
        # strictly nested and we are not doing any cuts.  However, we do
        # construct a root grid!
        root_level_grids = self.base_pf.h.select_grids(self.min_level)
        base_grid = ConstructedRootGrid(self.base_pf, self.pf, self,
                        self.min_level, self.min_left, self.max_right)
        self._fill_grid_arrays(base_grid, 0)
        grids = [base_grid]
        # We need to ensure we have the correct parentage relationships
        # However, we want the parent/child to be to the new proxy grids
        # so we need to map between the old ids and the new ids
        self.id_map = {}
        grid_id = 2 # id 0 is the base grid
        for level in range(self.min_level+1, self.max_level):
            for grid in self.base_pf.h.select_grids(level):
                # This next little bit will have to be changed if we ever move to
                # not-strictly-nested AMR hierarchies
                parent = self.id_map.get(grid.Parent.id, base_grid)
                grids.append(self.grid(grid_id, grid, self))
                parent.Children.append(grids[-1])
                grids[-1].Parent = parent
                self.id_map[grid.id] = grids[-1]
                # Now we fill in our arrays of values -- note that we
                # are filling in values from the base grids, not the newly
                # extracted grids.  We will perform bulk changes after we
                # finish.
                self._fill_grid_arrays(grid, grid_id-1)
                grid_id += 1

        self.grid_left_edge = self._convert_coords(self.grid_left_edge)
        self.grid_right_edge = self._convert_coords(self.grid_right_edge)
        self.grids = np.array(grids, dtype='object')

    def _fill_grid_arrays(self, grid, i):
        # This just fills in the grid arrays for a single grid --
        # note that we presuppose here that we are being handed a grid object
        # that has these defined; this means we are being handed the *base*
        # grid, not the newly extracted one
        self.grid_dimensions[i,:] = grid.ActiveDimensions
        self.grid_left_edge[i,:] = grid.LeftEdge
        self.grid_right_edge[i,:] = grid.RightEdge
        self.grid_particle_count[i] = grid.NumberOfParticles

    def _populate_grid_objects(self):
        for grid in self.grids:
            grid.Level = grid.base_grid.Level - self.pf.min_level
            grid._prepare_grid()
            grid._setup_dx()
            grid.start_index = None
        self.max_level -= self.pf.min_level
        print "New max level:", self.max_level

    def _convert_coords(self, val):
        return (val - self.left_edge_offset)*self.mult_factor

    def _detect_fields(self):
        self.field_list = self.base_pf.h.field_list[:]

    def _setup_unknown_fields(self):
        pass # Done in the base_h

    def _setup_derived_fields(self):
        self.derived_field_list = self.base_pf.h.derived_field_list[:]

    def _initialize_data_storage(self):
        self._data_file = None

    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        AMRHierarchy._setup_classes(self, dd)
        self.object_types.sort()

class ExtractedParameterFile(StaticOutput):
    _hierarchy_class = ExtractedHierarchy
    data_style = "extracted"
    
    def __init__(self, base_pf, min_level, max_level = -1, offset = None):
        self.base_pf = base_pf
        self.min_level = min_level
        self.max_level = max_level
        self.offset = offset
        self.override = {}

    def __repr__(self):
        return "extracted_%s" % self.base_pf

    def __getattr__(self, name):
        # This won't get called if 'name' is found already
        # and we'd like it to raise AttributeError if it's not anywhere
        if name in ['h', 'hierarchy']:
            return StaticOutput._get_hierarchy(self)
        return getattr(self.base_pf, name)

    def __getitem__(self, key):
        if key not in self.override:
            return self.base_pf[key]
        return self.override[key]

