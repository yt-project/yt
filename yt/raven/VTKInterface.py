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

from enthought.tvtk.pyface.ui.wx.wxVTKRenderWindowInteractor \
     import wxVTKRenderWindowInteractor

wxVTKRenderWindowInteractor.USE_STEREO = 1

class IVTKScene(object):
    def __init__(self):
        window = ivtk.IVTKWithCrustAndBrowser(size=(800,600), stereo=True)
        window.open()
        self.window = window
        self.scene = window.scene

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

    def __init__(self, vmin, vmax, vdefault, **traits):
        HasTraits.__init__(self, **traits)
        trait = Range(float(vmin), float(vmax), value=vdefault)
        self.add_trait("value", trait)
        self.value = vdefault

    def _value_changed(self, old, new):
        self.cubes.set_value(0, new)
        self.post_call()


class ExtractedVTKHierarchicalDataSet(HasTraits):

    parameter_fn = File(filter=["*.hierarchy"])
    base_grid_level = Int(22)
    field = Str("Density")
    scene_frame = Instance(IVTKScene)
    center = CArray(shape = (3,), dtype = 'float64')
    
    
    _gid = 0 # AMRBoxes require unique ids
    _grid_boundaries_shown = False
    _grid_boundaries_actor = None

    def _scene_frame_default(self):
        return IVTKScene()

    def _center_default(self):
        return [0.5,0.5,0.5]

    def _parameter_fn_changed(self, 
                 scene_frame = None, center = None):
        pf = lagos.EnzoStaticOutput(self.parameter_fn[:-10])
        base_grid = pf.h.select_grids(self.base_grid_level).tolist()
        center = pf.h.find_max("Density")[1]
        self.scene = self.scene_frame.scene
        self._take_log = True
        self._grids = []
        self._vtk_objs = []
        self._ugs = []
        self._xs = []
        self._ys = []
        self._zs = []
        self._vals = []
        self.operators = []
        self._oleft_edge = na.min([grid.LeftEdge for grid in base_grid], axis=0)
        self._base_level = base_grid[0].Level
        self._hdata_set = tvtk.HierarchicalBoxDataSet()
        self._mult_factor = 2**base_grid[0].Level
        self._min_val = 1e60
        self._max_val = -1e60
        self.left_edge = na.zeros(3, dtype='float64')
        #self.right_edge = (base_grid.RightEdge - base_grid.LeftEdge)*self._mult_factor
        self.right_edge = (na.max([grid.RightEdge for grid in base_grid], axis=0) -
                                self._oleft_edge) * self._mult_factor
        if center is None: center = (base_grid.RightEdge - base_grid.LeftEdge)/2.0
        print center, self._oleft_edge, self._mult_factor
        self.center = (center - self._oleft_edge)*self._mult_factor
        for grid in base_grid:
            self._add_grid(grid)
        self._hdata_set.generate_visibility_arrays()
        self.scene.camera.focal_point = self.center

    def _add_grid(self, grid):
        # Some of this will get wrapped into a get-vertex-centered-data routine
        # on the grid object
        # We recalculate here; easier than munging get_global_startindex
        if grid in self._grids: return
        self._grids.append(grid)
        left_index = (grid.LeftEdge - self._oleft_edge)/grid.dx
        right_index = left_index + grid.ActiveDimensions - 1
        # We need vertex-centered, so we need smoothed ghost zones and then we
        # interpolate them to the vertices
        cg = grid.retrieve_ghost_zones(2, self.field, smoothed=True)
        # Bounds should be cell-centered
        bds = na.array(zip(cg.left_edge+cg.dx/2.0, cg.right_edge-cg.dx/2.0)).ravel()
        interp = lagos.TrilinearFieldInterpolator(cg[self.field], bds, ['x','y','z'])
        dx = grid.dx * self._mult_factor
        origin=(grid.LeftEdge - self._oleft_edge)*self._mult_factor
        ad = grid.ActiveDimensions + 1
        ug = tvtk.UniformGrid(origin=origin, spacing=[dx]*3, dimensions=ad)
        # LE and RE are vertex locations
        x,y,z = na.mgrid[grid.LeftEdge[0]:grid.RightEdge[0]:ad[0]*1j,
                         grid.LeftEdge[1]:grid.RightEdge[1]:ad[1]*1j,
                         grid.LeftEdge[2]:grid.RightEdge[2]:ad[2]*1j]
        dd = {'x':x,'y':y,'z':z}
        if self._take_log: scalars = na.log10(interp(dd)).transpose().copy()
        else: scalars = interp(dd).transpose().copy()
        self._xs.append(x.ravel())
        self._ys.append(y.ravel())
        self._zs.append(z.ravel())
        self._vals.append(scalars.ravel())
        self._min_val = min(self._min_val, scalars.min())
        self._max_val = max(self._max_val, scalars.max())
        ug.point_data.scalars = scalars.ravel()
        ug.point_data.scalars.name = self.field
        self._ugs.append(ug)
        self._hdata_set.set_data_set(grid.Level-self._base_level, self._gid,
                                     left_index, right_index, ug)
        self._gid +=1
        # This is cheap, so we can set it every time
        self._hdata_set.set_refinement_ratio(grid.Level-self._base_level, 2)
        # Now we recurse
        # We should have a set of masks for grids that have been added, to
        # allow support for multiple children of a single parent
        for child in grid.Children:
            self._add_grid(child)

    def zoom(self, dist, unit='1'):
        vec = self.scene.camera.focal_point - \
              self.scene.camera.position
        self.scene.camera.position += \
            vec * dist/self._grids[0].pf[unit]
        self.scene.render()

    def toggle_grid_boundaries(self):
        if self._grid_boundaries_shown:
            self._grid_boundaries_actor.visibility = False
            return
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
            self._grid_boundaries_actor.visibility = True

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

    def add_x_plane(self):
        np, lookup_table = self._add_plane(self.center, normal=(1,0,0))
        self.operators.append(MappingPlane(
                vmin=self.left_edge[0], vmax=self.right_edge[0],
                vdefault = self.center[0],
                post_call = self.scene.render,
                plane = np, axis=0, coord=0.0,
                lookup_table = lookup_table))
        return self.operators[-1]

    def add_y_plane(self):
        np, lookup_table = self._add_plane(self.center, normal=(0,1,0))
        self.operators.append(MappingPlane(
                vmin=self.left_edge[1], vmax=self.right_edge[1],
                vdefault = self.center[1],
                post_call = self.scene.render,
                plane = np, axis=1, coord=0.0,
                lookup_table = lookup_table))
        return self.operators[-1]

    def add_z_plane(self):
        np, lookup_table = self._add_plane(self.center, normal=(0,0,1))
        self.operators.append(MappingPlane(
                vmin=self.left_edge[2], vmax=self.right_edge[2],
                vdefault = self.center[2],
                post_call = self.scene.render,
                plane = np, axis=2, coord=0.0,
                lookup_table = lookup_table))
        return self.operators[-1]

    def add_contour(self, val):
        cubes = tvtk.MarchingCubes(
                    executive = tvtk.CompositeDataPipeline())
        cubes.input = self._hdata_set
        cubes.set_value(0, val)
        cube_lut = tvtk.LookupTable(hue_range=(0.0,1.0),
                                    range=(self._min_val,
                                           self._max_val),
                                    table_range=(self._min_val,
                                                 self._max_val))
        print cube_lut.table_range, cube_lut.range
        cube_lut.build()
        cube_mapper = tvtk.HierarchicalPolyDataMapper(
                                input_connection = cubes.output_port,
                                lookup_table=cube_lut)

        cube_actor = tvtk.Actor(mapper=cube_mapper)
        self.scene.add_actors(cube_actor)
        self.operators.append(MappingMarchingCubes(cubes=cubes,
                    vmin=self._min_val, vmax=self._max_val,
                    vdefault=val,
                    post_call = self.scene.render,
                    lookup_table = cube_lut))
        return self.operators[-1]

    def add_volren(self):
        otf = tvtk.PiecewiseFunction()
        otf.add_point(self._min_val,0.0)
        otf.add_point(self._max_val,1.0)
        ctf = tvtk.ColorTransferFunction()
        vs = na.mgrid[self._min_val:self._max_val:5j]
        ctf.add_rgb_point(vs[0], 0.0, 0.0, 0.0)
        ctf.add_rgb_point(vs[1], 1.0, 0.0, 0.0)
        ctf.add_rgb_point(vs[2], 0.0, 0.0, 1.0)
        ctf.add_rgb_point(vs[3], 0.0, 1.0, 0.0)
        ctf.add_rgb_point(vs[4], 0.0, 0.2, 0.0)
        vp = tvtk.VolumeProperty()
        vp.set_scalar_opacity(otf)
        vp.set_color(ctf)
        cf = tvtk.VolumeRayCastCompositeFunction()
        self.cf = cf
        self.ctf = ctf
        self.otf = otf
        vals = na.concatenate(self._vals)
        xs = na.concatenate(self._xs)
        ys = na.concatenate(self._ys)
        zs = na.concatenate(self._zs)
        uug = tvtk.UnstructuredGrid()
        uug.point_data.t_coords = na.array([xs,ys,zs]).transpose()
        uug.point_data.scalars = vals
        #vmap = tvtk.FixedPointVolumeRayCastMapper(input=uug)
        vmap = tvtk.UnstructuredGridVolumeRayCastMapper(input=uug)
        #vmap.set_volume_raycastFunction
        #vmap = tvtk.VolumeTextureMapper3D()
        #vmap.sample_distance = self._grids[i].dx
        #vmap.volume_ray_cast_function = cf
        #vmap.input = uug
        volume = tvtk.Volume(mapper=vmap, property=vp)
        self.scene.add_actor(volume)

    def _convert_coords(self, val):
        return (self.val - self._oleft_edge)*self._mult_factor

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
    ehds = ExtractedVTKHierarchicalDataSet()
    ehds.configure_traits()
    ehds.toggle_grid_boundaries()
    gui.start_event_loop()
