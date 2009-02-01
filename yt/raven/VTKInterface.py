"""
This is the preliminary interface to VTK.  Note that as of VTK 5.2, it still
requires a patchset prepared here:
http://yt.enzotools.org/files/vtk_composite_data.zip

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2007-2009 Matthew Turk.  All Rights Reserved.

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

from enthought.tvtk.tools import ivtk
from enthought.tvtk.api import tvtk 
from enthought.traits.api import Float, HasTraits, Instance, Range, Any, \
                                 Delegate, Tuple, File, Int, Str, CArray

#from yt.reason import *
import sys
import numpy as na
import yt.lagos as lagos
from yt.funcs import *
from yt.logger import ravenLogger as mylog
from yt.extensions.HierarchySubset import ExtractedHierarchy

from enthought.tvtk.pyface.ui.wx.wxVTKRenderWindowInteractor \
     import wxVTKRenderWindowInteractor

#wxVTKRenderWindowInteractor.USE_STEREO = 1

class TVTKMapperWidget(HasTraits):
    lookup_table = Instance(tvtk.LookupTable)
    alpha_range = Tuple(Float(1.0), Float(1.0))
    post_call = Any

    def _alpha_range_changed(self, old, new):
        self.lookup_table.alpha_range = new
        self.post_call()

class MappingPlane(TVTKMapperWidget):
    plane = Instance(tvtk.Plane)

    def __init__(self, vmin, vmax, vdefault, **traits):
        HasTraits.__init__(self, **traits)
        trait = Range(float(vmin), float(vmax), value=vdefault)
        self.add_trait("coord", trait)
        self.coord = vdefault

    def _coord_changed(self, old, new):
        orig = self.plane.origin[:]
        orig[self.axis] = new
        self.plane.origin = orig
        self.post_call()

class MappingMarchingCubes(TVTKMapperWidget):
    cubes = Instance(tvtk.MarchingCubes)
    mapper = Instance(tvtk.HierarchicalPolyDataMapper)

    def __init__(self, vmin, vmax, vdefault, **traits):
        HasTraits.__init__(self, **traits)
        trait = Range(float(vmin), float(vmax), value=vdefault)
        self.add_trait("value", trait)
        self.value = vdefault

    def _value_changed(self, old, new):
        self.cubes.set_value(0, new)
        self.post_call()

class YTVTKScene(object):
    _grid_boundaries_shown = False
    _grid_boundaries_actor = None

    def __init__(self):
        window = ivtk.IVTKWithCrustAndBrowser(size=(800,600), stereo=1)
        window.open()
        self.window = window
        self.scene = window.scene

class YTScene(HasTraits):

    parameter_fn = File(filter=["*.hierarchy"])
    min_grid_level = Int(0)
    max_grid_level = Int(-1)
    field = Str("Density")
    scene_frame = Instance(YTVTKScene)
    center = CArray(shape = (3,), dtype = 'float64')
    _grid_boundaries_actor = None
    
    def _center_default(self):
        return [0.5,0.5,0.5]

    def _scene_frame_default(self):
        return YTVTKScene()

    def _parameter_fn_changed(self, scene_frame = None, center = None):
        self.pf = lagos.EnzoStaticOutput(self.parameter_fn[:-10])
        self.scene = self.scene_frame.scene
        self.extracted_hierarchy = ExtractedHierarchy(
                        self.pf, self.min_grid_level, self.max_grid_level,
                        offset=None)
        self._take_log = True
        self._grids = []
        self._ugs = []
        self._vtk_objs = []
        self.operators = []
        self._hdata_set = tvtk.HierarchicalBoxDataSet()
        self._min_val = 1e60
        self._max_val = -1e60
        self.center = self.extracted_hierarchy._convert_coords(
            self.pf.h.find_max("Density")[1])
        gid = 0
        for l, grid_set in enumerate(self.extracted_hierarchy.get_levels()):
            gid = self._add_level(grid_set, l, gid)
        self._hdata_set.generate_visibility_arrays()
        self.toggle_grid_boundaries()
        self.cubs = self.add_contour()
        self.cubs.configure_traits()
        self.scene.camera.focal_point = self.center

    def _add_level(self, grid_set, level, gid):
        for grid in grid_set:
            self._hdata_set.set_refinement_ratio(level, 2)
            gid = self._add_grid(grid, gid, level)
        return gid

    def _add_grid(self, grid, gid, level=0):
        mylog.debug("Adding grid %s on level %s (%s)",
                    grid.id, level, grid.Level)
        if grid in self._grids: return
        self._grids.append(grid)

        scalars = grid.get_vertex_centered_data(self.field)

        io, left_index, origin, dds = \
            self.extracted_hierarchy._convert_grid(grid)
        right_index = left_index + scalars.shape - 2
        ug = tvtk.UniformGrid(origin=origin, spacing=dds,
                              dimensions=grid.ActiveDimensions+1)
        if self.field not in self.pf.field_info or \
            self.pf.field_info[self.field].take_log:
            scalars = na.log10(scalars)
        ug.point_data.scalars = scalars.transpose().ravel()
        ug.point_data.scalars.name = self.field
        self._ugs.append(ug)
        self._hdata_set.set_data_set(level, gid, left_index, right_index, ug)

        self._min_val = min(self._min_val, scalars.min())
        self._max_val = max(self._max_val, scalars.max())

        gid += 1
        return gid

    def zoom(self, dist, unit='1'):
        vec = self.scene.camera.focal_point - \
              self.scene.camera.position
        self.scene.camera.position += \
            vec * dist/self._grids[0].pf[unit]
        self.scene.render()

    def toggle_grid_boundaries(self):
        if self._grid_boundaries_actor is None:
            # We don't need to track this stuff right now.
            ocf = tvtk.OutlineCornerFilter(
                    executive=tvtk.CompositeDataPipeline(),
                    corner_factor = 0.5)
            ocf.input = self._hdata_set
            ocm = tvtk.HierarchicalPolyDataMapper(
                    input_connection = ocf.output_port)
            self._grid_boundaries_actor = tvtk.Actor(mapper = ocm)
            self.scene.add_actor(self._grid_boundaries_actor)
        else:
            self._grid_boundaries_actor.visibility = \
            (not self._grid_boundaries_actor.visibility)

    def _add_plane(self, origin=(0.0,0.0,0.0), normal=(0,1,0)):
        plane = tvtk.Plane(origin=origin, normal=normal)
        cutter = tvtk.Cutter(executive = tvtk.CompositeDataPipeline(),
                             cut_function = plane)
        cutter.input = self._hdata_set
        clut = tvtk.LookupTable(hue_range=(0.0,1.00))
        clut.build()
        smap = tvtk.HierarchicalPolyDataMapper(
                        scalar_range=(self._min_val, self._max_val),
                        lookup_table=clut,
                        input_connection = cutter.output_port)
        sactor = tvtk.Actor(mapper=smap)
        self.scene.add_actors(sactor)
        return plane, clut

    def add_plane(self, origin=(0.0,0.0,0.0), normal=(0,1,0)):
        self.operators.append(self._add_plane(origin, normal))
        return self.operators[-1]

    def _add_axis_plane(self, axis):
        normal = [0,0,0]
        normal[axis] = 1
        np, lookup_table = self._add_plane(self.center, normal=normal)
        LE = self.extracted_hierarchy._convert_coords(
                self.extracted_hierarchy.left_edge_offset)
        RE = self.extracted_hierarchy._convert_coords(
                self.extracted_hierarchy.right_edge_offset)
        self.operators.append(MappingPlane(
                vmin=LE[axis], vmax=RE[axis],
                vdefault = self.center[axis],
                post_call = self.scene.render,
                plane = np, axis=axis, coord=0.0,
                lookup_table = lookup_table))

    def add_x_plane(self):
        self._add_axis_plane(0)
        return self.operators[-1]

    def add_y_plane(self):
        self._add_axis_plane(1)
        return self.operators[-1]

    def add_z_plane(self):
        self._add_axis_plane(2)
        return self.operators[-1]

    def add_contour(self, val=None):
        if val is None: val = (self._max_val+self._min_val) * 0.5
        cubes = tvtk.MarchingCubes(
                    executive = tvtk.CompositeDataPipeline())
        cubes.input = self._hdata_set
        cubes.set_value(0, val)
        cube_lut = tvtk.LookupTable(hue_range=(0.0,1.0))
        cube_lut.build()
        cube_mapper = tvtk.HierarchicalPolyDataMapper(
                                input_connection = cubes.output_port,
                                lookup_table=cube_lut)
        cube_mapper.color_mode = 'map_scalars'
        cube_mapper.scalar_range = (self._min_val, self._max_val)
        cube_actor = tvtk.Actor(mapper=cube_mapper)
        self.scene.add_actors(cube_actor)
        self.operators.append(MappingMarchingCubes(cubes=cubes,
                    vmin=self._min_val, vmax=self._max_val,
                    vdefault=val,
                    mapper = cube_mapper,
                    post_call = self.scene.render,
                    lookup_table = cube_lut))
        return self.operators[-1]

def get_all_parents(grid):
    parents = []
    if len(grid.Parents) == 0: return grid
    for parent in grid.Parents: parents.append(get_all_parents(parent))
    return list(set(parents))

if __name__=="__main__":
    print "This code probably won't work.  You need to install the patchset to VTK,"
    print "the Enthought Tool Suite, and then run th TVTK test."
    print
    print "If you've done all that, remove the 'sys.exit()' line and try again."
    sys.exit()
    import yt.lagos as lagos

    from enthought.pyface.api import GUI
    gui = GUI()
    ehds = YTScene()
    ehds.configure_traits()
    ehds.toggle_grid_boundaries()
    gui.start_event_loop()
