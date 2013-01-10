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
    parallel_objects

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
        xd = self.pf.domain_dimensions[x_dict[self.axis]]
        yd = self.pf.domain_dimensions[y_dict[self.axis]]
        return QuadTree(np.array([xd,yd], dtype='int64'), nvals,
                        style = self.proj_style)

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
        YTSelectionContainer3D.__init__(self,
            field_parameters.get("center", None),
            pf, field_parameters)
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
