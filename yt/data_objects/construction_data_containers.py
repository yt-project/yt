"""
Data containers that require processing before they can be utilized.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Author: Britton Smith <Britton.Smith@colorado.edu>
Affiliation: University of Colorado at Boulder
Homepage: http://yt-project.org/
License:
  Copyright (C) 2007-2011 Matthew Turk.  All Rights Reserved.

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

import numpy as np
import math
import weakref
import exceptions
import itertools
import shelve
from exceptions import ValueError, KeyError
from functools import wraps

from yt.funcs import *
from yt.utilities.logger import ytLogger
from .data_containers import \
    YTSelectionContainer1D, YTSelectionContainer2D, YTSelectionContainer3D, \
    restore_field_information_state
from .field_info_container import \
    NeedsOriginalGrid
from yt.utilities.lib import \
    QuadTree, ghost_zone_interpolate, fill_region
from yt.utilities.data_point_utilities import CombineGrids,\
    DataCubeRefine, DataCubeReplace, FillRegion, FillBuffer
from yt.utilities.definitions import axis_names, x_dict, y_dict
from yt.utilities.minimal_representation import \
    MinimalProjectionData
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    parallel_objects, parallel_root_only, ParallelAnalysisInterface

from .field_info_container import\
    NeedsGridType,\
    NeedsOriginalGrid,\
    NeedsDataField,\
    NeedsProperty,\
    NeedsParameter

class YTStreamlineBase(YTSelectionContainer1D):
    _type_name = "streamline"
    _con_args = ('positions')
    sort_by = 't'
    def __init__(self, positions, length = 1.0, fields=None, pf=None, **kwargs):
        """
        This is a streamline, which is a set of points defined as
        being parallel to some vector field.

        This object is typically accessed through the Streamlines.path
        function.  The resulting arrays have their dimensionality
        reduced to one, and an ordered list of points at an (x,y)
        tuple along `axis` are available, as is the `t` field, which
        corresponds to a unitless measurement along the ray from start
        to end.

        Parameters
        ----------
        positions : array-like
            List of streamline positions
        length : float
            The magnitude of the distance; dts will be divided by this
        fields : list of strings, optional
            If you want the object to pre-retrieve a set of fields, supply them
            here.  This is not necessary.
        pf : Parameter file object
            Passed in to access the hierarchy
        kwargs : dict of items
            Any additional values are passed as field parameters that can be
            accessed by generated fields.

        Examples
        --------

        >>> from yt.visualization.api import Streamlines
        >>> streamlines = Streamlines(pf, [0.5]*3) 
        >>> streamlines.integrate_through_volume()
        >>> stream = streamlines.path(0)
        >>> matplotlib.pylab.semilogy(stream['t'], stream['Density'], '-x')
        
        """
        YTSelectionContainer1D.__init__(self, pf, fields, **kwargs)
        self.positions = positions
        self.dts = np.empty_like(positions[:,0])
        self.dts[:-1] = np.sqrt(np.sum((self.positions[1:]-
                                        self.positions[:-1])**2,axis=1))
        self.dts[-1] = self.dts[-1]
        self.length = length
        self.dts /= length
        self.ts = np.add.accumulate(self.dts)
        self._set_center(self.positions[0])
        self.set_field_parameter('center', self.positions[0])
        self._dts, self._ts = {}, {}
        #self._refresh_data()

    def _get_list_of_grids(self):
        # Get the value of the line at each LeftEdge and RightEdge
        LE = self.pf.h.grid_left_edge
        RE = self.pf.h.grid_right_edge
        # Check left faces first
        min_streampoint = np.min(self.positions, axis=0)
        max_streampoint = np.max(self.positions, axis=0)
        p = np.all((min_streampoint <= RE) & (max_streampoint > LE), axis=1)
        self._grids = self.hierarchy.grids[p]

    def _get_data_from_grid(self, grid, field):
        # No child masking here; it happens inside the mask cut
        mask = self._get_cut_mask(grid) 
        if field == 'dts': return self._dts[grid.id]
        if field == 't': return self._ts[grid.id]
        return grid[field].flat[mask]
        

    def _get_cut_mask(self, grid):
        #pdb.set_trace()
        points_in_grid = np.all(self.positions > grid.LeftEdge, axis=1) & \
                         np.all(self.positions <= grid.RightEdge, axis=1) 
        pids = np.where(points_in_grid)[0]
        mask = np.zeros(points_in_grid.sum(), dtype='int')
        dts = np.zeros(points_in_grid.sum(), dtype='float64')
        ts = np.zeros(points_in_grid.sum(), dtype='float64')
        for mi, (i, pos) in enumerate(zip(pids, self.positions[points_in_grid])):
            if not points_in_grid[i]: continue
            ci = ((pos - grid.LeftEdge)/grid.dds).astype('int')
            if grid.child_mask[ci[0], ci[1], ci[2]] == 0: continue
            for j in range(3):
                ci[j] = min(ci[j], grid.ActiveDimensions[j]-1)
            mask[mi] = np.ravel_multi_index(ci, grid.ActiveDimensions)
            dts[mi] = self.dts[i]
            ts[mi] = self.ts[i]
        self._dts[grid.id] = dts
        self._ts[grid.id] = ts
        return mask


class YTQuadTreeProjBase(YTSelectionContainer2D):
    _key_fields = YTSelectionContainer2D._key_fields + ['weight_field']
    _type_name = "proj"
    _con_args = ('axis', 'weight_field')
    _container_fields = ('px', 'py', 'pdx', 'pdy')
    def __init__(self, field, axis, weight_field = None,
                 center = None, pf = None, data_source=None, 
                 style = "integrate", field_parameters = None):
        """
        This is a data object corresponding to a line integral through the
        simulation domain.

        This object is typically accessed through the `proj` object that
        hangs off of hierarchy objects.  AMRQuadProj is a projection of a
        `field` along an `axis`.  The field can have an associated
        `weight_field`, in which case the values are multiplied by a weight
        before being summed, and then divided by the sum of that weight; the
        two fundamental modes of operating are direct line integral (no
        weighting) and average along a line of sight (weighting.)  What makes
        `proj` different from the standard projection mechanism is that it
        utilizes a quadtree data structure, rather than the old mechanism for
        projections.  It will not run in parallel, but serial runs should be
        substantially faster.  Note also that lines of sight are integrated at
        every projected finest-level cell.

        Parameters
        ----------
        field : string
            This is the field which will be "projected" along the axis.  If
            multiple are specified (in a list) they will all be projected in
            the first pass.
        axis : int
            The axis along which to slice.  Can be 0, 1, or 2 for x, y, z.
        weight_field : string
            If supplied, the field being projected will be multiplied by this
            weight value before being integrated, and at the conclusion of the
            projection the resultant values will be divided by the projected
            `weight_field`.
        max_level : int
            If supplied, only cells at or below this level will be projected.
        center : array_like, optional
            The 'center' supplied to fields that use it.  Note that this does
            not have to have `coord` as one value.  Strictly optional.
        source : `yt.data_objects.api.YTDataContainer`, optional
            If specified, this will be the data source used for selecting
            regions to project.
        node_name: string, optional
            The node in the .yt file to find or store this slice at.  Should
            probably not be used.
        field_cuts : list of strings, optional
            If supplied, each of these strings will be evaluated to cut a
            region of a grid out.  They can be of the form "grid['Temperature']
            > 100" for instance.
        preload_style : string
            Either 'level' (default) or 'all'.  Defines how grids are loaded --
            either level by level, or all at once.  Only applicable during
            parallel runs.
        serialize : bool, optional
            Whether we should store this projection in the .yt file or not.
        kwargs : dict of items
            Any additional values are passed as field parameters that can be
            accessed by generated fields.

        Examples
        --------

        >>> pf = load("RedshiftOutput0005")
        >>> qproj = pf.h.proj("Density", 0)
        >>> print qproj["Density"]
        """
        YTSelectionContainer2D.__init__(self, axis, pf, field_parameters)
        self.proj_style = style
        if style == "mip":
            self.func = np.max
            op = "max"
        elif style == "integrate":
            self.func = np.sum # for the future
            op = "sum"
        else:
            raise NotImplementedError(style)
        self.weight_field = weight_field
        self._set_center(center)
        if center is not None: self.set_field_parameter('center',center)
        if data_source is None: data_source = self.pf.h.all_data()
        self.data_source = data_source
        self.weight_field = weight_field
        self.get_data(field)

    @property
    def _mrep(self):
        return MinimalProjectionData(self)

    def hub_upload(self):
        self._mrep.upload()

    def _get_tree(self, nvals):
        xax = x_dict[self.axis]
        yax = y_dict[self.axis]
        xd = self.pf.domain_dimensions[xax]
        yd = self.pf.domain_dimensions[yax]
        bounds = (self.pf.domain_left_edge[xax],
                  self.pf.domain_right_edge[yax],
                  self.pf.domain_left_edge[xax],
                  self.pf.domain_right_edge[yax])
        return QuadTree(np.array([xd,yd], dtype='int64'), nvals,
                        bounds, style = self.proj_style)

    def _get_conv(self, fields):
        # Place holder for a time when maybe we will not be doing just
        # a single dx for every field.
        convs = np.empty(len(fields), dtype="float64")
        fields = self._determine_fields(fields)
        for i, field in enumerate(fields):
            fi = self.pf._get_field_info(*field)
            convs[i] = (self.pf.units[fi.projection_conversion])
        return convs

    def get_data(self, fields = None):
        fields = self._determine_fields(ensure_list(fields))
        # We need a new tree for every single set of fields we add
        if len(fields) == 0: return
        chunk_fields = fields[:]
        if self.weight_field is not None:
            chunk_fields.append(self.weight_field)
        tree = self._get_tree(len(fields))
        # We do this once
        for chunk in self.data_source.chunks(None, "io"):
            self._initialize_chunk(chunk, tree)
        # This needs to be parallel_objects-ified
        for chunk in parallel_objects(self.data_source.chunks(
                chunk_fields, "io")): 
            mylog.debug("Adding chunk (%s) to tree", chunk.size)
            self._handle_chunk(chunk, fields, tree)
        # Note that this will briefly double RAM usage
        if self.proj_style == "mip":
            merge_style = -1
            op = "max"
        elif self.proj_style == "integrate":
            merge_style = 1
            op = "sum"
        else:
            raise NotImplementedError
        # TODO: Add the combine operation
        ox = self.pf.domain_left_edge[x_dict[self.axis]]
        oy = self.pf.domain_left_edge[y_dict[self.axis]]
        px, py, pdx, pdy, nvals, nwvals = tree.get_all(False)
        nvals = self.comm.mpi_allreduce(nvals, op=op)
        nwvals = self.comm.mpi_allreduce(nwvals, op=op)
        np.multiply(px, self.pf.domain_width[x_dict[self.axis]], px)
        np.add(px, ox, px)
        np.multiply(pdx, self.pf.domain_width[x_dict[self.axis]], pdx)

        np.multiply(py, self.pf.domain_width[y_dict[self.axis]], py)
        np.add(py, oy, py)
        np.multiply(pdy, self.pf.domain_width[y_dict[self.axis]], pdy)

        if self.weight_field is not None:
            np.divide(nvals, nwvals[:,None], nvals)
        if self.weight_field is None:
            convs = self._get_conv(fields)
            nvals *= convs[None,:]
        # We now convert to half-widths and center-points
        data = {}
        data['px'] = px
        data['py'] = py
        data['weight_field'] = nwvals
        data['pdx'] = pdx
        data['pdy'] = pdy
        data['fields'] = nvals
        # Now we run the finalizer, which is ignored if we don't need it
        fd = data['fields']
        field_data = np.hsplit(data.pop('fields'), len(fields))
        for fi, field in enumerate(fields):
            mylog.debug("Setting field %s", field)
            self[field] = field_data[fi].ravel()
        for i in data.keys(): self[i] = data.pop(i)
        mylog.info("Projection completed")

    def _initialize_chunk(self, chunk, tree):
        icoords = chunk.icoords
        i1 = icoords[:,0]
        i2 = icoords[:,1]
        ilevel = chunk.ires
        tree.initialize_chunk(i1, i2, ilevel)

    def _handle_chunk(self, chunk, fields, tree):
        if self.proj_style == "mip":
            dl = 1.0
        else:
            dl = chunk.fwidth[:, self.axis]
        v = np.empty((chunk.size, len(fields)), dtype="float64")
        for i in range(len(fields)):
            v[:,i] = chunk[fields[i]] * dl
        if self.weight_field is not None:
            w = chunk[self.weight_field]
            np.multiply(v, w[:,None], v)
            np.multiply(w, dl, w)
        else:
            w = np.ones(chunk.size, dtype="float64")
        icoords = chunk.icoords
        i1 = icoords[:,x_dict[self.axis]]
        i2 = icoords[:,y_dict[self.axis]]
        ilevel = chunk.ires
        tree.add_chunk_to_tree(i1, i2, ilevel, v, w)

    def to_pw(self, fields=None, center='c', width=None, axes_unit=None, 
               origin='center-window'):
        r"""Create a :class:`~yt.visualization.plot_window.PWViewerMPL` from this
        object.

        This is a bare-bones mechanism of creating a plot window from this
        object, which can then be moved around, zoomed, and on and on.  All
        behavior of the plot window is relegated to that routine.
        """
        pw = self._get_pw(fields, center, width, origin, axes_unit, 'Projection')
        return pw

class YTCoveringGridBase(YTSelectionContainer3D):
    _spatial = True
    _type_name = "covering_grid"
    _con_args = ('level', 'left_edge', 'ActiveDimensions')
    _base_grid = None
    def __init__(self, level, left_edge, dims, fields = None,
                 pf = None, num_ghost_zones = 0, use_pbar = True, 
                 field_parameters = None):
        """A 3D region with all data extracted to a single, specified
        resolution.

        Parameters
        ----------
        level : int
            The resolution level data is uniformly gridded at
        left_edge : array_like
            The left edge of the region to be extracted
        dims : array_like
            Number of cells along each axis of resulting covering_grid
        fields : array_like, optional
            A list of fields that you'd like pre-generated for your object

        Examples
        --------
        cube = pf.h.covering_grid(2, left_edge=[0.0, 0.0, 0.0], \
                                  dims=[128, 128, 128])

        """
        if field_parameters is None:
            center = None
        else:
            center = field_parameters.get("center", None)
        YTSelectionContainer3D.__init__(self,
            center, pf, field_parameters)
        self.left_edge = np.array(left_edge)
        self.level = level
        rdx = self.pf.domain_dimensions*self.pf.refine_by**level
        self.dds = self.pf.domain_width/rdx.astype("float64")
        self.ActiveDimensions = np.array(dims, dtype='int32')
        self.right_edge = self.left_edge + self.ActiveDimensions*self.dds
        self._num_ghost_zones = num_ghost_zones
        self._use_pbar = use_pbar
        self.global_startindex = np.rint((self.left_edge-self.pf.domain_left_edge)/self.dds).astype('int64')
        self.domain_width = np.rint((self.pf.domain_right_edge -
                    self.pf.domain_left_edge)/self.dds).astype('int64')
        self._setup_data_source()

    def _setup_data_source(self):
        self._data_source = self.pf.h.region(
            self.center, self.left_edge, self.right_edge)
        self._data_source.min_level = 0
        self._data_source.max_level = self.level

    def get_data(self, fields = None):
        fields = self._determine_fields(ensure_list(fields))
        fields_to_get = [f for f in fields if f not in self.field_data]
        fields_to_get = self._identify_dependencies(fields_to_get)
        if len(fields_to_get) == 0: return
        fill, gen = self._split_fields(fields_to_get)
        if len(fill) > 0: self._fill_fields(fill)
        if len(gen) > 0: self._generate_fields(gen)

    def _split_fields(self, fields_to_get):
        fill, gen = self.pf.h._split_fields(fields_to_get)
        for field in gen:
            finfo = self.pf._get_field_info(*field)
            try:
                finfo.check_available(self)
            except NeedsOriginalGrid:
                fill.append(field)
        gen = [f for f in gen if f not in fill]
        return fill, gen

    def _fill_fields(self, fields):
        output_fields = [np.zeros(self.ActiveDimensions, dtype="float64")
                         for field in fields]
        for chunk in self._data_source.chunks(fields, "io"):
            input_fields = [chunk[field] for field in fields]
            fill_region(input_fields, output_fields, self.level,
                        self.global_startindex, chunk.icoords, chunk.ires)
        for name, v in zip(fields, output_fields):
            self[name] = v

class LevelState(object):
    current_dx = None
    current_dims = None
    current_level = None
    global_startindex = None
    old_global_startindex = None
    domain_iwidth = None
    fields = None
    data_source = None

class YTSmoothedCoveringGridBase(YTCoveringGridBase):
    _type_name = "smoothed_covering_grid"
    filename = None
    @wraps(YTCoveringGridBase.__init__)
    def __init__(self, *args, **kwargs):
        """A 3D region with all data extracted and interpolated to a
        single, specified resolution.  (Identical to covering_grid,
        except that it interpolates.)

        Smoothed covering grids start at level 0, interpolating to
        fill the region to level 1, replacing any cells actually
        covered by level 1 data, and then recursively repeating this
        process until it reaches the specified `level`.

        Parameters
        ----------
        level : int
            The resolution level data is uniformly gridded at
        left_edge : array_like
            The left edge of the region to be extracted
        dims : array_like
            Number of cells along each axis of resulting covering_grid
        fields : array_like, optional
            A list of fields that you'd like pre-generated for your object

        Example
        -------
        cube = pf.h.smoothed_covering_grid(2, left_edge=[0.0, 0.0, 0.0], \
                                  dims=[128, 128, 128])
        """
        self._base_dx = (
              (self.pf.domain_right_edge - self.pf.domain_left_edge) /
               self.pf.domain_dimensions.astype("float64"))
        self.global_endindex = None
        YTCoveringGridBase.__init__(self, *args, **kwargs)
        self._final_start_index = self.global_startindex

    def _setup_data_source(self, level_state = None):
        if level_state is None: return
        # We need a buffer region to allow for zones that contribute to the
        # interpolation but are not directly inside our bounds
        level_state.data_source = self.pf.h.region(
            self.center,
            self.left_edge - level_state.current_dx,
            self.right_edge + level_state.current_dx)
        level_state.data_source.min_level = level_state.current_level
        level_state.data_source.max_level = level_state.current_level

    def _fill_fields(self, fields):
        ls = self._initialize_level_state(fields)
        for level in range(self.level + 1):
            tot = 0
            for chunk in ls.data_source.chunks(fields, "io"):
                chunk[fields[0]]
                input_fields = [chunk[field] for field in fields]
                tot += fill_region(input_fields, ls.fields, ls.current_level,
                            ls.global_startindex, chunk.icoords,
                            chunk.ires)
            self._update_level_state(ls)
        for name, v in zip(fields, ls.fields):
            if self.level > 0: v = v[1:-1,1:-1,1:-1]
            self[name] = v

    def _initialize_level_state(self, fields):
        ls = LevelState()
        ls.current_dx = self._base_dx
        ls.current_level = 0
        LL = self.left_edge - self.pf.domain_left_edge
        ls.global_startindex = np.rint(LL / ls.current_dx).astype('int64') - 1
        ls.domain_iwidth = np.rint((self.pf.domain_right_edge -
                    self.pf.domain_left_edge)/ls.current_dx).astype('int64')
        if self.level > 0:
            # We use one grid cell at LEAST, plus one buffer on all sides
            width = self.right_edge-self.left_edge
            idims = np.rint(width/ls.current_dx).astype('int64') + 2
        elif self.level == 0:
            ls.global_startindex = np.array(np.floor(LL/ls.current_dx), dtype='int64')
            idims = np.rint((self.ActiveDimensions*self.dds)/ls.current_dx).astype('int64')
        ls.current_dims = idims.astype("int32")
        ls.fields = [np.zeros(idims, dtype="float64")-999 for field in fields]
        self._setup_data_source(ls)
        return ls

    def _update_level_state(self, level_state):
        ls = level_state
        if ls.current_level >= self.level: return
        ls.current_level += 1
        ls.current_dx = self._base_dx / self.pf.refine_by**ls.current_level
        self._setup_data_source(ls)
        LL = self.left_edge - self.pf.domain_left_edge
        ls.old_global_startindex = ls.global_startindex
        ls.global_startindex = np.rint(LL / ls.current_dx).astype('int64') - 1
        ls.domain_iwidth = np.rint(self.pf.domain_width/ls.current_dx).astype('int64') 
        rf = float(self.pf.refine_by)
        input_left = (level_state.old_global_startindex + 0.5) * rf 
        width = (self.ActiveDimensions*self.dds)
        output_dims = np.rint(width/level_state.current_dx+0.5).astype("int32") + 2
        level_state.current_dims = output_dims
        new_fields = []
        for input_field in level_state.fields:
            output_field = np.zeros(output_dims, dtype="float64")
            output_left = self.global_startindex + 0.5
            ghost_zone_interpolate(rf, input_field, input_left,
                                   output_field, output_left)
            new_fields.append(output_field)
        level_state.fields = new_fields

class YTSurfaceBase(YTSelectionContainer3D, ParallelAnalysisInterface):
    _type_name = "surface"
    _con_args = ("data_source", "surface_field", "field_value")
    vertices = None
    def __init__(self, data_source, surface_field, field_value):
        r"""This surface object identifies isocontours on a cell-by-cell basis,
        with no consideration of global connectedness, and returns the vertices
        of the Triangles in that isocontour.

        This object simply returns the vertices of all the triangles
        calculated by the marching cubes algorithm; for more complex
        operations, such as identifying connected sets of cells above a given
        threshold, see the extract_connected_sets function.  This is more
        useful for calculating, for instance, total isocontour area, or
        visualizing in an external program (such as `MeshLab
        <http://meshlab.sf.net>`_.)  The object has the properties .vertices
        and will sample values if a field is requested.  The values are
        interpolated to the center of a given face.
        
        Parameters
        ----------
        data_source : AMR3DDataObject
            This is the object which will used as a source
        surface_field : string
            Any field that can be obtained in a data object.  This is the field
            which will be isocontoured.
        field_value : float
            The value at which the isocontour should be calculated.

        References
        ----------

        .. [1] Marching Cubes: http://en.wikipedia.org/wiki/Marching_cubes

        Examples
        --------
        This will create a data object, find a nice value in the center, and
        output the vertices to "triangles.obj" after rescaling them.

        >>> sp = pf.h.sphere("max", (10, "kpc")
        >>> surf = pf.h.surface(sp, "Density", 5e-27)
        >>> print surf["Temperature"]
        >>> print surf.vertices
        >>> bounds = [(sp.center[i] - 5.0/pf['kpc'],
        ...            sp.center[i] + 5.0/pf['kpc']) for i in range(3)]
        >>> surf.export_ply("my_galaxy.ply", bounds = bounds)
        """
        ParallelAnalysisInterface.__init__(self)
        self.data_source = data_source
        self.surface_field = surface_field
        self.field_value = field_value
        self.vertex_samples = YTFieldData()
        center = data_source.get_field_parameter("center")
        AMRData.__init__(self, center = center, fields = None, pf =
                         data_source.pf)
        self._grids = self.data_source._grids.copy()

    def get_data(self, fields = None, sample_type = "face"):
        if isinstance(fields, list) and len(fields) > 1:
            for field in fields: self.get_data(field)
            return
        elif isinstance(fields, list):
            fields = fields[0]
        # Now we have a "fields" value that is either a string or None
        pb = get_pbar("Extracting (sampling: %s)" % fields,
                      len(list(self._get_grid_objs())))
        verts = []
        samples = []
        for i,g in enumerate(self._get_grid_objs()):
            pb.update(i)
            my_verts = self._extract_isocontours_from_grid(
                            g, self.surface_field, self.field_value,
                            fields, sample_type)
            if fields is not None:
                my_verts, svals = my_verts
                samples.append(svals)
            verts.append(my_verts)
        pb.finish()
        verts = np.concatenate(verts).transpose()
        verts = self.comm.par_combine_object(verts, op='cat', datatype='array')
        self.vertices = verts
        if fields is not None:
            samples = np.concatenate(samples)
            samples = self.comm.par_combine_object(samples, op='cat',
                                datatype='array')
            if sample_type == "face":
                self[fields] = samples
            elif sample_type == "vertex":
                self.vertex_samples[fields] = samples
        

    def _extract_isocontours_from_grid(self, grid, field, value,
                                       sample_values = None,
                                       sample_type = "face"):
        mask = self.data_source._get_cut_mask(grid) * grid.child_mask
        vals = grid.get_vertex_centered_data(field, no_ghost = False)
        if sample_values is not None:
            svals = grid.get_vertex_centered_data(sample_values)
        else:
            svals = None
        sample_type = {"face":1, "vertex":2}[sample_type]
        my_verts = march_cubes_grid(value, vals, mask, grid.LeftEdge,
                                    grid.dds, svals, sample_type)
        return my_verts

    def calculate_flux(self, field_x, field_y, field_z, fluxing_field = None):
        r"""This calculates the flux over the surface.

        This function will conduct marching cubes on all the cells in a given
        data container (grid-by-grid), and then for each identified triangular
        segment of an isocontour in a given cell, calculate the gradient (i.e.,
        normal) in the isocontoured field, interpolate the local value of the
        "fluxing" field, the area of the triangle, and then return:

        area * local_flux_value * (n dot v)

        Where area, local_value, and the vector v are interpolated at the barycenter
        (weighted by the vertex values) of the triangle.  Note that this
        specifically allows for the field fluxing across the surface to be
        *different* from the field being contoured.  If the fluxing_field is
        not specified, it is assumed to be 1.0 everywhere, and the raw flux
        with no local-weighting is returned.

        Additionally, the returned flux is defined as flux *into* the surface,
        not flux *out of* the surface.
        
        Parameters
        ----------
        field_x : string
            The x-component field
        field_y : string
            The y-component field
        field_z : string
            The z-component field
        fluxing_field : string, optional
            The field whose passage over the surface is of interest.  If not
            specified, assumed to be 1.0 everywhere.

        Returns
        -------
        flux : float
            The summed flux.  Note that it is not currently scaled; this is
            simply the code-unit area times the fields.

        References
        ----------

        .. [1] Marching Cubes: http://en.wikipedia.org/wiki/Marching_cubes

        Examples
        --------

        This will create a data object, find a nice value in the center, and
        calculate the metal flux over it.

        >>> sp = pf.h.sphere("max", (10, "kpc")
        >>> surf = pf.h.surface(sp, "Density", 5e-27)
        >>> flux = surf.calculate_flux(
        ...     "x-velocity", "y-velocity", "z-velocity", "Metal_Density")
        """
        flux = 0.0
        pb = get_pbar("Fluxing %s" % fluxing_field,
                len(list(self._get_grid_objs())))
        for i, g in enumerate(self._get_grid_objs()):
            pb.update(i)
            flux += self._calculate_flux_in_grid(g,
                    field_x, field_y, field_z, fluxing_field)
        pb.finish()
        flux = self.comm.mpi_allreduce(flux, op="sum")
        return flux

    def _calculate_flux_in_grid(self, grid, 
                    field_x, field_y, field_z, fluxing_field = None):
        mask = self.data_source._get_cut_mask(grid) * grid.child_mask
        vals = grid.get_vertex_centered_data(self.surface_field)
        if fluxing_field is None:
            ff = np.ones(vals.shape, dtype="float64")
        else:
            ff = grid.get_vertex_centered_data(fluxing_field)
        xv, yv, zv = [grid.get_vertex_centered_data(f) for f in 
                     [field_x, field_y, field_z]]
        return march_cubes_grid_flux(self.field_value, vals, xv, yv, zv,
                    ff, mask, grid.LeftEdge, grid.dds)

    def export_ply(self, filename, bounds = None, color_field = None,
                   color_map = "algae", color_log = True, sample_type = "face"):
        r"""This exports the surface to the PLY format, suitable for visualization
        in many different programs (e.g., MeshLab).

        Parameters
        ----------
        filename : string
            The file this will be exported to.  This cannot be a file-like object.
        bounds : list of tuples
            The bounds the vertices will be normalized to.  This is of the format:
            [(xmin, xmax), (ymin, ymax), (zmin, zmax)].  Defaults to the full
            domain.
        color_field : string
            Should a field be sample and colormapped?
        color_map : string
            Which color map should be applied?
        color_log : bool
            Should the color field be logged before being mapped?

        Examples
        --------

        >>> sp = pf.h.sphere("max", (10, "kpc")
        >>> surf = pf.h.surface(sp, "Density", 5e-27)
        >>> print surf["Temperature"]
        >>> print surf.vertices
        >>> bounds = [(sp.center[i] - 5.0/pf['kpc'],
        ...            sp.center[i] + 5.0/pf['kpc']) for i in range(3)]
        >>> surf.export_ply("my_galaxy.ply", bounds = bounds)
        """
        if self.vertices is None:
            self.get_data(color_field, sample_type)
        elif color_field is not None:
            if sample_type == "face" and \
                color_field not in self.field_data:
                self[color_field]
            elif sample_type == "vertex" and \
                color_field not in self.vertex_data:
                self.get_data(color_field, sample_type)
        self._export_ply(filename, bounds, color_field, color_map, color_log,
                         sample_type)

    def _color_samples(self, cs, color_log, color_map, arr):
            if color_log: cs = np.log10(cs)
            mi, ma = cs.min(), cs.max()
            cs = (cs - mi) / (ma - mi)
            from yt.visualization.image_writer import map_to_colors
            cs = map_to_colors(cs, color_map)
            arr["red"][:] = cs[0,:,0]
            arr["green"][:] = cs[0,:,1]
            arr["blue"][:] = cs[0,:,2]

    @parallel_root_only
    def _export_ply(self, filename, bounds = None, color_field = None,
                   color_map = "algae", color_log = True, sample_type = "face"):
        if isinstance(filename, file):
            f = filename
        else:
            f = open(filename, "wb")
        if bounds is None:
            DLE = self.pf.domain_left_edge
            DRE = self.pf.domain_right_edge
            bounds = [(DLE[i], DRE[i]) for i in range(3)]
        nv = self.vertices.shape[1]
        vs = [("x", "<f"), ("y", "<f"), ("z", "<f"),
              ("red", "uint8"), ("green", "uint8"), ("blue", "uint8") ]
        fs = [("ni", "uint8"), ("v1", "<i4"), ("v2", "<i4"), ("v3", "<i4"),
              ("red", "uint8"), ("green", "uint8"), ("blue", "uint8") ]
        f.write("ply\n")
        f.write("format binary_little_endian 1.0\n")
        f.write("element vertex %s\n" % (nv))
        f.write("property float x\n")
        f.write("property float y\n")
        f.write("property float z\n")
        if color_field is not None and sample_type == "vertex":
            f.write("property uchar red\n")
            f.write("property uchar green\n")
            f.write("property uchar blue\n")
            v = np.empty(self.vertices.shape[1], dtype=vs)
            cs = self.vertex_samples[color_field]
            self._color_samples(cs, color_log, color_map, v)
        else:
            v = np.empty(self.vertices.shape[1], dtype=vs[:3])
        f.write("element face %s\n" % (nv/3))
        f.write("property list uchar int vertex_indices\n")
        if color_field is not None and sample_type == "face":
            f.write("property uchar red\n")
            f.write("property uchar green\n")
            f.write("property uchar blue\n")
            # Now we get our samples
            cs = self[color_field]
            arr = np.empty(cs.shape[0], dtype=np.dtype(fs))
            self._color_samples(cs, color_log, color_map, arr)
        else:
            arr = np.empty(nv/3, np.dtype(fs[:-3]))
        for i, ax in enumerate("xyz"):
            # Do the bounds first since we cast to f32
            tmp = self.vertices[i,:]
            np.subtract(tmp, bounds[i][0], tmp)
            w = bounds[i][1] - bounds[i][0]
            np.divide(tmp, w, tmp)
            np.subtract(tmp, 0.5, tmp) # Center at origin.
            v[ax][:] = tmp 
        f.write("end_header\n")
        v.tofile(f)
        arr["ni"][:] = 3
        vi = np.arange(nv, dtype="<i")
        vi.shape = (nv/3, 3)
        arr["v1"][:] = vi[:,0]
        arr["v2"][:] = vi[:,1]
        arr["v3"][:] = vi[:,2]
        arr.tofile(f)
        if filename is not f:
            f.close()

    def export_sketchfab(self, title, description, api_key = None,
                            color_field = None, color_map = "algae",
                            color_log = True, bounds = None):
        r"""This exports Surfaces to SketchFab.com, where they can be viewed
        interactively in a web browser.

        SketchFab.com is a proprietary web service that provides WebGL
        rendering of models.  This routine will use temporary files to
        construct a compressed binary representation (in .PLY format) of the
        Surface and any optional fields you specify and upload it to
        SketchFab.com.  It requires an API key, which can be found on your
        SketchFab.com dashboard.  You can either supply the API key to this
        routine directly or you can place it in the variable
        "sketchfab_api_key" in your ~/.yt/config file.  This function is
        parallel-safe.

        Parameters
        ----------
        title : string
            The title for the model on the website
        description : string
            How you want the model to be described on the website
        api_key : string
            Optional; defaults to using the one in the config file
        color_field : string
            If specified, the field by which the surface will be colored
        color_map : string
            The name of the color map to use to map the color field
        color_log : bool
            Should the field be logged before being mapped to RGB?
        bounds : list of tuples
            [ (xmin, xmax), (ymin, ymax), (zmin, zmax) ] within which the model
            will be scaled and centered.  Defaults to the full domain.

        Returns
        -------
        URL : string
            The URL at which your model can be viewed.

        Examples
        --------

        from yt.mods import *
        pf = load("redshift0058")
        dd = pf.h.sphere("max", (200, "kpc"))
        rho = 5e-27

        bounds = [(dd.center[i] - 100.0/pf['kpc'],
                   dd.center[i] + 100.0/pf['kpc']) for i in range(3)]

        surf = pf.h.surface(dd, "Density", rho)

        rv = surf.export_sketchfab(
            title = "Testing Upload",
            description = "A simple test of the uploader",
            color_field = "Temperature",
            color_map = "hot",
            color_log = True,
            bounds = bounds
        )
        """
        api_key = api_key or ytcfg.get("yt","sketchfab_api_key")
        if api_key in (None, "None"):
            raise YTNoAPIKey("SketchFab.com", "sketchfab_api_key")
        import zipfile, json
        from tempfile import TemporaryFile

        ply_file = TemporaryFile()
        self.export_ply(ply_file, bounds, color_field, color_map, color_log,
                        sample_type = "vertex")
        ply_file.seek(0)
        # Greater than ten million vertices and we throw an error but dump
        # to a file.
        if self.vertices.shape[1] > 1e7:
            tfi = 0
            fn = "temp_model_%03i.ply" % tfi
            while os.path.exists(fn):
                fn = "temp_model_%03i.ply" % tfi
                tfi += 1
            open(fn, "wb").write(ply_file.read())
            raise YTTooManyVertices(self.vertices.shape[1], fn)

        zfs = TemporaryFile()
        with zipfile.ZipFile(zfs, "w", zipfile.ZIP_DEFLATED) as zf:
            zf.writestr("yt_export.ply", ply_file.read())
        zfs.seek(0)

        zfs.seek(0)
        data = {
            'title': title,
            'token': api_key,
            'description': description,
            'fileModel': zfs,
            'filenameModel': "yt_export.zip",
        }
        upload_id = self._upload_to_sketchfab(data)
        upload_id = self.comm.mpi_bcast(upload_id, root = 0)
        return upload_id

    @parallel_root_only
    def _upload_to_sketchfab(self, data):
        import urllib2, json
        from yt.utilities.poster.encode import multipart_encode
        from yt.utilities.poster.streaminghttp import register_openers
        register_openers()
        datamulti, headers = multipart_encode(data)
        request = urllib2.Request("https://api.sketchfab.com/v1/models",
                        datamulti, headers)
        rv = urllib2.urlopen(request).read()
        rv = json.loads(rv)
        upload_id = rv.get("result", {}).get("id", None)
        if upload_id:
            mylog.info("Model uploaded to: https://sketchfab.com/show/%s",
                       upload_id)
        else:
            mylog.error("Problem uploading.")
        return upload_id


