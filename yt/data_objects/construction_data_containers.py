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
    QuadTree, ghost_zone_interpolate
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

class YTOverlapProjBase(YTSelectionContainer2D):
    _top_node = "/Projections"
    _key_fields = YTSelectionContainer2D._key_fields + ['weight_field']
    _type_name = "overlap_proj"
    _con_args = ('axis', 'weight_field')
    def __init__(self, field, axis, weight_field = None,
                 max_level = None, center = None, pf = None,
                 source=None, node_name = None, field_cuts = None,
                 preload_style=None, serialize=True, **kwargs):
        """
        This is a data object corresponding to a line integral through the
        simulation domain.

        This object is typically accessed through the `proj` object that
        hangs off of hierarchy objects.  AMRProj is a projection of a `field`
        along an `axis`.  The field can have an associated `weight_field`, in
        which case the values are multiplied by a weight before being summed,
        and then divided by the sum of that weight; the two fundamental modes
        of operating are direct line integral (no weighting) and average along
        a line of sight (weighting.)  Note also that lines of sight are
        integrated at every projected finest-level cell

        Parameters
        ----------
        field : string
            This is the field which will be "projected" along the axis.  If
            multiple are specified (in a list) they will all be projected in
            the first pass.
        axis : int or axis
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
            Either 'level', 'all', or None (default).  Defines how grids are
            loaded -- either level by level, or all at once.  Only applicable
            during parallel runs.
        serialize : bool, optional
            Whether we should store this projection in the .yt file or not.
        kwargs : dict of items
            Any additional values are passed as field parameters that can be
            accessed by generated fields.

        Examples
        --------

        >>> pf = load("RedshiftOutput0005")
        >>> proj = pf.h.proj("Density", "x")
        >>> print proj["Density"]
        """
        YTSelectionContainer2D.__init__(self, axis, field, pf, node_name = None, **kwargs)
        self.weight_field = weight_field
        self._field_cuts = field_cuts
        self.serialize = serialize
        self._set_center(center)
        if center is not None: self.set_field_parameter('center',center)
        self._node_name = node_name
        self._initialize_source(source)
        self._grids = self.source._grids
        if max_level == None:
            max_level = self.hierarchy.max_level
        if self.source is not None:
            max_level = min(max_level, self.source.grid_levels.max())
        self._max_level = max_level
        self._weight = weight_field
        self.preload_style = preload_style
        self.func = np.sum # for the future
        self.__retval_coords = {}
        self.__retval_fields = {}
        self.__retval_coarse = {}
        self.__overlap_masks = {}
        self._deserialize(node_name)
        self._refresh_data()
        if self._okay_to_serialize and self.serialize: self._serialize(node_name=self._node_name)

    def _convert_field_name(self, field):
        if field == "weight_field": return "weight_field_%s" % self._weight
        if field in self._key_fields: return field
        return "%s_%s" % (field, self._weight)

    def _initialize_source(self, source = None):
        if source is None:
            check, source = self.partition_hierarchy_2d(self.axis)
            self._check_region = check
            #self._okay_to_serialize = (not check)
        else:
            self._distributed = False
            self._okay_to_serialize = False
            self._check_region = True
        self.source = source
        if self._field_cuts is not None:
            # Override if field cuts are around; we don't want to serialize!
            self._check_region = True
            self._okay_to_serialize = False
        if self._node_name is not None:
            self._node_name = "%s/%s" % (self._top_node,self._node_name)
            self._okay_to_serialize = True

    def __calculate_overlap(self, level):
        s = self.source
        mylog.info("Generating overlap masks for level %s", level)
        i = 0
        pbar = get_pbar("Reading and masking grids ", len(s._grids))
        mylog.debug("Examining level %s", level)
        grids = s.select_grid_indices(level)
        RE = s.grid_right_edge[grids]
        LE = s.grid_left_edge[grids]
        for grid in s._grids[grids]:
            pbar.update(i)
            self.__overlap_masks[grid.id] = \
                grid._generate_overlap_masks(self.axis, LE, RE)
            i += 1
        pbar.finish()
        mylog.info("Finished calculating overlap.")

    def _get_dls(self, grid, fields):
        # Place holder for a time when maybe we will not be doing just
        # a single dx for every field.
        dls = []
        convs = []
        for field in fields + [self._weight]:
            if field is None: continue
            dls.append(just_one(grid['d%s' % axis_names[self.axis]]))
            convs.append(self.pf.units[self.pf.field_info[field].projection_conversion])
        return np.array(dls), np.array(convs)

    def __project_level(self, level, fields):
        grids_to_project = self.source.select_grids(level)
        dls, convs = self._get_dls(grids_to_project[0], fields)
        zero_out = (level != self._max_level)
        pbar = get_pbar('Projecting  level % 2i / % 2i ' \
                          % (level, self._max_level), len(grids_to_project))
        for pi, grid in enumerate(grids_to_project):
            g_coords, g_fields = self._project_grid(grid, fields, zero_out)
            self.__retval_coords[grid.id] = g_coords
            self.__retval_fields[grid.id] = g_fields
            for fi in range(len(fields)): g_fields[fi] *= dls[fi]
            if self._weight is not None: g_coords[3] *= dls[-1]
            pbar.update(pi)
            grid.clear_data()
        pbar.finish()
        self.__combine_grids_on_level(level) # In-place
        if level > 0 and level <= self._max_level:
            self.__refine_to_level(level) # In-place
        coord_data = []
        field_data = []
        for grid in grids_to_project:
            coarse = self.__retval_coords[grid.id][2]==0 # Where childmask = 0
            fine = ~coarse
            coord_data.append([pi[fine] for pi in self.__retval_coords[grid.id]])
            field_data.append([pi[fine] for pi in self.__retval_fields[grid.id]])
            self.__retval_coords[grid.id] = [pi[coarse] for pi in self.__retval_coords[grid.id]]
            self.__retval_fields[grid.id] = [pi[coarse] for pi in self.__retval_fields[grid.id]]
        coord_data = np.concatenate(coord_data, axis=1)
        field_data = np.concatenate(field_data, axis=1)
        if self._weight is not None:
            field_data = field_data / coord_data[3,:].reshape((1,coord_data.shape[1]))
        else:
            field_data *= convs[...,np.newaxis]
        mylog.info("Level %s done: %s final", \
                   level, coord_data.shape[1])
        pdx = grids_to_project[0].dds[x_dict[self.axis]] # this is our dl
        pdy = grids_to_project[0].dds[y_dict[self.axis]] # this is our dl
        return coord_data, pdx, pdy, field_data

    def __combine_grids_on_level(self, level):
        grids = self.source.select_grids(level)
        grids_i = self.source.select_grid_indices(level)
        pbar = get_pbar('Combining   level % 2i / % 2i ' \
                          % (level, self._max_level), len(grids))
        # We have an N^2 check, so we try to be as quick as possible
        # and to skip as many as possible
        for pi, grid1 in enumerate(grids):
            pbar.update(pi)
            if self.__retval_coords[grid1.id][0].shape[0] == 0: continue
            for grid2 in self.source._grids[grids_i][self.__overlap_masks[grid1.id]]:
                if self.__retval_coords[grid2.id][0].shape[0] == 0 \
                  or grid1.id == grid2.id:
                    continue
                args = [] # First is source, then destination
                args += self.__retval_coords[grid2.id] + [self.__retval_fields[grid2.id]]
                args += self.__retval_coords[grid1.id] + [self.__retval_fields[grid1.id]]
                args.append(1) # Refinement factor
                args.append(np.ones(args[0].shape, dtype='int64'))
                kk = CombineGrids(*args)
                goodI = args[-1].astype('bool')
                self.__retval_coords[grid2.id] = \
                    [coords[goodI] for coords in self.__retval_coords[grid2.id]]
                self.__retval_fields[grid2.id] = \
                    [fields[goodI] for fields in self.__retval_fields[grid2.id]]
        pbar.finish()

    def __refine_to_level(self, level):
        grids = self.source.select_grids(level)
        grids_up = self.source.select_grid_indices(level - 1)
        pbar = get_pbar('Refining to level % 2i / % 2i ' \
                          % (level, self._max_level), len(grids))
        for pi, grid1 in enumerate(grids):
            pbar.update(pi)
            for parent in ensure_list(grid1.Parent):
                if parent.id not in self.__overlap_masks: continue
                for grid2 in self.source._grids[grids_up][self.__overlap_masks[parent.id]]:
                    if self.__retval_coords[grid2.id][0].shape[0] == 0: continue
                    args = []
                    args += self.__retval_coords[grid2.id] + [self.__retval_fields[grid2.id]]
                    args += self.__retval_coords[grid1.id] + [self.__retval_fields[grid1.id]]
                    # Refinement factor, which is same in all directions.  Note
                    # that this complicated rounding is because sometimes
                    # epsilon differences in dds between the grids causes this
                    # to round to up or down from the expected value.
                    args.append(int(np.rint(grid2.dds / grid1.dds)[0]))
                    args.append(np.ones(args[0].shape, dtype='int64'))
                    kk = CombineGrids(*args)
                    goodI = args[-1].astype('bool')
                    self.__retval_coords[grid2.id] = \
                        [coords[goodI] for coords in self.__retval_coords[grid2.id]]
                    self.__retval_fields[grid2.id] = \
                        [fields[goodI] for fields in self.__retval_fields[grid2.id]]
        for grid1 in self.source.select_grids(level-1):
            if not self._check_region and self.__retval_coords[grid1.id][0].size != 0:
                mylog.error("Something messed up, and %s still has %s points of data",
                            grid1, self.__retval_coords[grid1.id][0].size)
                mylog.error("Please contact the yt-users mailing list.")
                raise ValueError(grid1, self.__retval_coords[grid1.id])
        pbar.finish()

    def get_data(self, fields):
        fields = ensure_list(fields)
        self._obtain_fields(fields, self._node_name)
        fields = [f for f in fields if f not in self.field_data]
        if len(fields) == 0: return
        coord_data = []
        field_data = []
        pdxs = []
        pdys = []
        # We do this here, but I am not convinced it should be done here
        # It is probably faster, as it consolidates IO, but if we did it in
        # _project_level, then it would be more memory conservative
        if self.preload_style == 'all':
            print "Preloading %s grids and getting %s" % (
                    len(self.source._grids), self.get_dependencies(fields))
            self.comm.preload(self.source._grids,
                          self.get_dependencies(fields), self.hierarchy.io)
        for level in range(0, self._max_level+1):
            if self.preload_style == 'level':
                self.comm.preload(self.source.select_grids(level),
                              self.get_dependencies(fields), self.hierarchy.io)
            self.__calculate_overlap(level)
            my_coords, my_pdx, my_pdy, my_fields = \
                self.__project_level(level, fields)
            coord_data.append(my_coords)
            field_data.append(my_fields)
            pdxs.append(my_pdx * np.ones(my_coords.shape[1], dtype='float64'))
            pdys.append(my_pdx * np.ones(my_coords.shape[1], dtype='float64'))
            if self._check_region and False:
                check=self.__cleanup_level(level - 1)
                if len(check) > 0: all_data.append(check)
            # Now, we should clean up after ourselves...
            for grid in self.source.select_grids(level - 1):
                del self.__retval_coords[grid.id]
                del self.__retval_fields[grid.id]
                del self.__overlap_masks[grid.id]
            mylog.debug("End of projecting level level %s, memory usage %0.3e",
                        level, get_memory_usage()/1024.)
        coord_data = np.concatenate(coord_data, axis=1)
        field_data = np.concatenate(field_data, axis=1)
        pdxs = np.concatenate(pdxs, axis=1)
        pdys = np.concatenate(pdys, axis=1)
        # We now convert to half-widths and center-points
        data = {}
        data['pdx'] = pdxs; del pdxs
        data['pdy'] = pdys; del pdys
        ox = self.pf.domain_left_edge[x_dict[self.axis]]
        oy = self.pf.domain_left_edge[y_dict[self.axis]]
        data['px'] = (coord_data[0,:]+0.5) * data['pdx'] + ox
        data['py'] = (coord_data[1,:]+0.5) * data['pdx'] + oy
        data['weight_field'] = coord_data[3,:].copy()
        del coord_data
        data['pdx'] *= 0.5
        data['pdy'] *= 0.5
        data['fields'] = field_data
        # Now we run the finalizer, which is ignored if we don't need it
        data = self.comm.par_combine_object(data, datatype='dict', op='cat')
        field_data = np.vsplit(data.pop('fields'), len(fields))
        for fi, field in enumerate(fields):
            self[field] = field_data[fi].ravel()
            if self.serialize: self._store_fields(field, self._node_name)
        for i in data.keys(): self[i] = data.pop(i)
        mylog.info("Projection completed")

    def add_fields(self, fields, weight = "CellMassMsun"):
        pass

    def _project_grid(self, grid, fields, zero_out):
        # We split this next bit into two sections to try to limit the IO load
        # on the system.  This way, we perserve grid state (@restore_grid_state
        # in _get_data_from_grid *and* we attempt not to load weight data
        # independently of the standard field data.
        if self._weight is None:
            weight_data = np.ones(grid.ActiveDimensions, dtype='float64')
            if zero_out: weight_data[grid.child_indices] = 0
            masked_data = [fd.astype('float64') * weight_data
                           for fd in self._get_data_from_grid(grid, fields)]
        else:
            fields_to_get = list(set(fields + [self._weight]))
            field_data = dict(zip(
                fields_to_get, self._get_data_from_grid(grid, fields_to_get)))
            weight_data = field_data[self._weight].copy().astype('float64')
            if zero_out: weight_data[grid.child_indices] = 0
            masked_data  = [field_data[field].copy().astype('float64') * weight_data
                                for field in fields]
            del field_data
        # if we zero it out here, then we only have to zero out the weight!
        full_proj = [self.func(field, axis=self.axis) for field in masked_data]
        weight_proj = self.func(weight_data, axis=self.axis)
        if (self._check_region and not self.source._is_fully_enclosed(grid)) or self._field_cuts is not None:
            used_data = self._get_points_in_region(grid).astype('bool')
            used_points = np.where(np.logical_or.reduce(used_data, self.axis))
        else:
            used_data = np.array([1.0], dtype='bool')
            used_points = slice(None)
        if zero_out:
            subgrid_mask = np.logical_and.reduce(
                                np.logical_or(grid.child_mask,
                                             ~used_data),
                                self.axis).astype('int64')
        else:
            subgrid_mask = np.ones(full_proj[0].shape, dtype='int64')
        xind, yind = [arr[used_points].ravel() for arr in np.indices(full_proj[0].shape)]
        start_index = grid.get_global_startindex()
        xpoints = (xind + (start_index[x_dict[self.axis]])).astype('int64')
        ypoints = (yind + (start_index[y_dict[self.axis]])).astype('int64')
        return ([xpoints, ypoints,
                subgrid_mask[used_points].ravel(),
                weight_proj[used_points].ravel()],
                [data[used_points].ravel() for data in full_proj])

    def _get_points_in_region(self, grid):
        pointI = self.source._get_point_indices(grid, use_child_mask=False)
        point_mask = np.zeros(grid.ActiveDimensions)
        point_mask[pointI] = 1.0
        if self._field_cuts is not None:
            for cut in self._field_cuts:
                point_mask *= eval(cut)
        return point_mask

    def _get_data_from_grid(self, grid, fields):
        fields = ensure_list(fields)
        if self._check_region:
            bad_points = self._get_points_in_region(grid)
        else:
            bad_points = 1.0
        return [grid[field] * bad_points for field in fields]

    def _gen_node_name(self):
        return  "%s/%s" % \
            (self._top_node, self.axis)


class YTCoveringGridBase(YTSelectionContainer3D):
    _spatial = True
    _type_name = "covering_grid"
    _con_args = ('level', 'left_edge', 'ActiveDimensions')
    def __init__(self, level, left_edge, dims, fields = None,
                 pf = None, num_ghost_zones = 0, use_pbar = True, **kwargs):
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
        YTSelectionContainer3D.__init__(self, center=kwargs.pop("center", None),
                           pf=pf, **kwargs)
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
        self.get_data(fields)

    def _get_list_of_grids(self, buffer = 0.0):
        if self._grids is not None: return
        if np.any(self.left_edge - buffer < self.pf.domain_left_edge) or \
           np.any(self.right_edge + buffer > self.pf.domain_right_edge):
            grids,ind = self.pf.hierarchy.get_periodic_box_grids_below_level(
                            self.left_edge - buffer,
                            self.right_edge + buffer, self.level)
        else:
            grids,ind = self.pf.hierarchy.get_box_grids_below_level(
                self.left_edge - buffer,
                self.right_edge + buffer, self.level)
        sort_ind = np.argsort(self.pf.h.grid_levels.ravel()[ind])
        self._grids = self.pf.hierarchy.grids[ind][(sort_ind,)][::-1]

    def _refresh_data(self):
        YTSelectionContainer3D._refresh_data(self)
        self['dx'] = self.dds[0] * np.ones(self.ActiveDimensions, dtype='float64')
        self['dy'] = self.dds[1] * np.ones(self.ActiveDimensions, dtype='float64')
        self['dz'] = self.dds[2] * np.ones(self.ActiveDimensions, dtype='float64')

    def get_data(self, fields):
        if self._grids is None:
            self._get_list_of_grids()
        fields = ensure_list(fields)
        obtain_fields = []
        for field in fields:
            if self.field_data.has_key(field): continue
            if field not in self.hierarchy.field_list:
                try:
                    #print "Generating", field
                    self._generate_field(field)
                    continue
                except NeedsOriginalGrid, ngt_exception:
                    pass
            obtain_fields.append(field)
            self[field] = np.zeros(self.ActiveDimensions, dtype='float64') -999
        if len(obtain_fields) == 0: return
        mylog.debug("Getting fields %s from %s possible grids",
                   obtain_fields, len(self._grids))
        if self._use_pbar: pbar = \
                get_pbar('Searching grids for values ', len(self._grids))
        count = self.ActiveDimensions.prod()
        for i, grid in enumerate(self._grids):
            if self._use_pbar: pbar.update(i)
            count -= self._get_data_from_grid(grid, obtain_fields)
            if count <= 0: break
        if self._use_pbar: pbar.finish()
        if count > 0 or np.any(self[obtain_fields[0]] == -999):
            # and self.dx < self.hierarchy.grids[0].dx:
            n_bad = np.where(self[obtain_fields[0]]==-999)[0].size
            mylog.error("Covering problem: %s cells are uncovered", n_bad)
            raise KeyError(n_bad)

    def _generate_field(self, field):
        if self.pf.field_info.has_key(field):
            # First we check the validator; this might even raise!
            self.pf.field_info[field].check_available(self)
            self[field] = self.pf.field_info[field](self)
        else: # Can't find the field, try as it might
            raise KeyError(field)

    def flush_data(self, fields=None):
        """
        Any modifications made to the data in this object are pushed back
        to the originating grids, except the cells where those grids are both
        below the current level `and` have child cells.
        """
        self._get_list_of_grids()
        # We don't generate coordinates here.
        for grid in self._grids:
            self._flush_data_to_grid(grid, ensure_list(fields))

    def _get_data_from_grid(self, grid, fields):
        ll = int(grid.Level == self.level)
        ref_ratio = self.pf.refine_by**(self.level - grid.Level)
        g_fields = [grid[field].astype("float64") for field in fields]
        c_fields = [self[field] for field in fields]
        count = FillRegion(ref_ratio,
            grid.get_global_startindex(), self.global_startindex,
            c_fields, g_fields,
            self.ActiveDimensions, grid.ActiveDimensions,
            grid.child_mask, self.domain_width, ll, 0)
        return count

    def _flush_data_to_grid(self, grid, fields):
        ll = int(grid.Level == self.level)
        ref_ratio = self.pf.refine_by**(self.level - grid.Level)
        g_fields = []
        for field in fields:
            if not grid.has_key(field): grid[field] = \
               np.zeros(grid.ActiveDimensions, dtype=self[field].dtype)
            g_fields.append(grid[field])
        c_fields = [self[field] for field in fields]
        FillRegion(ref_ratio,
            grid.get_global_startindex(), self.global_startindex,
            c_fields, g_fields,
            self.ActiveDimensions, grid.ActiveDimensions,
            grid.child_mask, self.domain_width, ll, 1)

    @property
    def LeftEdge(self):
        return self.left_edge

    @property
    def RightEdge(self):
        return self.right_edge


class YTSmoothedCoveringGridBase(YTCoveringGridBase):
    _type_name = "smoothed_covering_grid"
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
        # We need a buffer region to allow for zones that contribute to the
        # interpolation but are not directly inside our bounds
        buffer = ((self.pf.domain_right_edge - self.pf.domain_left_edge)
                 / self.pf.domain_dimensions).max()
        self._base_region = self.pf.geometry.region(
            self.center,
            self.left_edge - buffer,
            self.right_edge + buffer)

    def get_data(self, field):
        self._get_list_of_grids()
        # We don't generate coordinates here.
        if field == None:
            fields = self.fields[:]
        else:
            fields = ensure_list(field)
        fields_to_get = []
        for field in fields:
            if self.field_data.has_key(field): continue
            if field not in self.hierarchy.field_list:
                try:
                    #print "Generating", field
                    self._generate_field(field)
                    continue
                except NeedsOriginalGrid, ngt_exception:
                    pass
            fields_to_get.append(field)
        if len(fields_to_get) == 0: return
        # Note that, thanks to some trickery, we have different dimensions
        # on the field than one might think from looking at the dx and the
        # L/R edges.
        # We jump-start our task here
        mylog.debug("Getting fields %s from %s possible grids",
                   fields_to_get, len(self._grids))
        self._update_level_state(0, fields_to_get)
        if self._use_pbar: pbar = \
                get_pbar('Searching grids for values ', len(self._grids))
        # The grids are assumed to be pre-sorted
        last_level = 0
        for gi, grid in enumerate(self._grids):
            if self._use_pbar: pbar.update(gi)
            if grid.Level > last_level and grid.Level <= self.level:
                mylog.debug("Updating level state to %s", last_level + 1)
                self._update_level_state(last_level + 1)
                self._refine(1, fields_to_get)
                last_level = grid.Level
            self._get_data_from_grid(grid, fields_to_get)
        while last_level < self.level:
            mylog.debug("Grid-free refinement %s to %s", last_level, last_level + 1)
            self._update_level_state(last_level + 1)
            self._refine(1, fields_to_get)
            last_level += 1
        if self.level > 0:
            for field in fields_to_get:
                self[field] = self[field][1:-1,1:-1,1:-1]
                if np.any(self[field] == -999):
                    # and self.dx < self.hierarchy.grids[0].dx:
                    n_bad = (self[field]==-999).sum()
                    mylog.error("Covering problem: %s cells are uncovered", n_bad)
                    raise KeyError(n_bad)
        if self._use_pbar: pbar.finish()

    def _update_level_state(self, level, fields = None):
        dx = self._base_dx / self.pf.refine_by**level
        self.field_data['cdx'] = dx[0]
        self.field_data['cdy'] = dx[1]
        self.field_data['cdz'] = dx[2]
        LL = self.left_edge - self.pf.domain_left_edge
        self._old_global_startindex = self.global_startindex
        self.global_startindex = np.rint(LL / dx).astype('int64') - 1
        self.domain_width = np.rint((self.pf.domain_right_edge -
                    self.pf.domain_left_edge)/dx).astype('int64')
        if level == 0 and self.level > 0:
            # We use one grid cell at LEAST, plus one buffer on all sides
            idims = np.rint((self.right_edge-self.left_edge)/dx).astype('int64') + 2
            fields = ensure_list(fields)
            for field in fields:
                self.field_data[field] = np.zeros(idims,dtype='float64')-999
            self._cur_dims = idims.astype("int32")
        elif level == 0 and self.level == 0:
            DLE = self.pf.domain_left_edge
            self.global_startindex = np.array(np.floor(LL/ dx), dtype='int64')
            idims = np.rint((self.ActiveDimensions*self.dds)/dx).astype('int64')
            fields = ensure_list(fields)
            for field in fields:
                self.field_data[field] = np.zeros(idims,dtype='float64')-999
            self._cur_dims = idims.astype("int32")

    def _refine(self, dlevel, fields):
        rf = float(self.pf.refine_by**dlevel)

        input_left = (self._old_global_startindex + 0.5) * rf 
        dx = np.fromiter((self['cd%s' % ax] for ax in 'xyz'), count=3, dtype='float64')
        output_dims = np.rint((self.ActiveDimensions*self.dds)/dx+0.5).astype('int32') + 2
        self._cur_dims = output_dims

        for field in fields:
            output_field = np.zeros(output_dims, dtype="float64")
            output_left = self.global_startindex + 0.5
            ghost_zone_interpolate(rf, self[field], input_left,
                                   output_field, output_left)
            self.field_data[field] = output_field

    @restore_field_information_state
    def _get_data_from_grid(self, grid, fields):
        fields = ensure_list(fields)
        g_fields = [grid[field].astype("float64") for field in fields]
        c_fields = [self.field_data[field] for field in fields]
        count = FillRegion(1,
            grid.get_global_startindex(), self.global_startindex,
            c_fields, g_fields, 
            self._cur_dims, grid.ActiveDimensions,
            grid.child_mask, self.domain_width, 1, 0)
        return count

    def flush_data(self, *args, **kwargs):
        raise KeyError("Can't do this")

class YTFixedResProjectionBase(YTSelectionContainer2D):
    _top_node = "/Projections"
    _type_name = "fixed_res_proj"
    _con_args = ('axis', 'field', 'weight_field')
    def __init__(self, axis, level, left_edge, dims, pf=None,
                 field_parameters = None):
        """
        This is a data structure that projects grids, but only to fixed (rather
        than variable) resolution.

        This object is typically accessed through the `fixed_res_proj` object
        that hangs off of hierarchy objects.  This projection mechanism is much
        simpler than the standard, variable-resolution projection.  Rather than
        attempt to identify the highest-resolution element along every possible
        line of sight, this data structure simply deposits every cell into one
        of a fixed number of bins.  It is suitable for inline analysis, and it
        should scale nicely.

        Parameters
        ----------
        axis : int
            The axis along which to project.  Can be 0, 1, or 2 for x, y, z.
        level : int
            This is the level to which values will be projected.  Note that
            the pixel size in the projection will be identical to a cell at
            this level of refinement in the simulation.
        left_edge : array of ints
            The left edge, in level-local integer coordinates, of the
            projection
        dims : array of ints
            The dimensions of the projection (which, in concert with the
            left_edge, serves to define its right edge.)
        field_parameters : dict of items
            Any additional values are passed as field parameters that can be
            accessed by generated fields.

        Examples
        --------

        >>> pf = load("RedshiftOutput0005")
        >>> fproj = pf.h.fixed_res_proj(1, [0, 0, 0], [64, 64, 64], ["Density"])
        >>> print fproj["Density"]
        """
        YTSelectionContainer2D.__init__(self, axis, pf, field_parameters)
        self.left_edge = np.array(left_edge)
        self.level = level
        self.dds = self.pf.h.select_grids(self.level)[0].dds.copy()
        self.dims = np.array([dims]*2)
        self.ActiveDimensions = np.array([dims]*3, dtype='int32')
        self.right_edge = self.left_edge + self.ActiveDimensions*self.dds
        self.global_startindex = np.rint((self.left_edge - self.pf.domain_left_edge)
                                         /self.dds).astype('int64')
        self._dls = {}
        self.domain_width = np.rint((self.pf.domain_right_edge -
                    self.pf.domain_left_edge)/self.dds).astype('int64')
        self._refresh_data()

    def _get_list_of_grids(self):
        if self._grids is not None: return
        if np.any(self.left_edge < self.pf.domain_left_edge) or \
           np.any(self.right_edge > self.pf.domain_right_edge):
            grids,ind = self.pf.hierarchy.get_periodic_box_grids(
                            self.left_edge, self.right_edge)
        else:
            grids,ind = self.pf.hierarchy.get_box_grids(
                            self.left_edge, self.right_edge)
        level_ind = (self.pf.hierarchy.grid_levels.ravel()[ind] <= self.level)
        sort_ind = np.argsort(self.pf.h.grid_levels.ravel()[ind][level_ind])
        self._grids = self.pf.hierarchy.grids[ind][level_ind][(sort_ind,)][::-1]

    def _generate_coords(self):
        xax = x_dict[self.axis]
        yax = y_dict[self.axis]
        ci = self.left_edge + self.dds*0.5
        cf = self.left_edge + self.dds*(self.ActiveDimensions-0.5)
        cx = np.mgrid[ci[xax]:cf[xax]:self.ActiveDimensions[xax]*1j]
        cy = np.mgrid[ci[yax]:cf[yax]:self.ActiveDimensions[yax]*1j]
        blank = np.ones( (self.ActiveDimensions[xax],
                          self.ActiveDimensions[yax]), dtype='float64')
        self['px'] = cx[None,:] * blank
        self['py'] = cx[:,None] * blank
        self['pdx'] = self.dds[xax]
        self['pdy'] = self.dds[yax]

    #@time_execution
    def get_data(self, fields):
        """
        Iterates over the list of fields and generates/reads them all.
        """
        self._get_list_of_grids()
        if not self.has_key('pdx'):
            self._generate_coords()
        fields_to_get = ensure_list(fields)
        if len(fields_to_get) == 0: return
        temp_data = {}
        for field in fields_to_get:
            self[field] = np.zeros(self.dims, dtype='float64')
        dls = self.__setup_dls(fields_to_get)
        for i,grid in enumerate(self._get_grids()):
            mylog.debug("Getting fields from %s", i)
            self._get_data_from_grid(grid, fields_to_get, dls)
        mylog.info("IO completed; summing")
        for field in fields_to_get:
            self[field] = self.comm.mpi_allreduce(self[field], op='sum')
            conv = self.pf.units[self.pf.field_info[field].projection_conversion]
            self[field] *= conv

    def __setup_dls(self, fields):
        dls = {}
        for level in range(self.level+1):
            dls[level] = []
            grid = self.select_grids(level)[0]
            for field in fields:
                if field is None: continue
                dls[level].append(float(just_one(grid['d%s' % axis_names[self.axis]])))
        return dls

    #@restore_grid_state
    def _get_data_from_grid(self, grid, fields, dls):
        g_fields = [grid[field].astype("float64") for field in fields]
        c_fields = [self[field] for field in fields]
        ref_ratio = self.pf.refine_by**(self.level - grid.Level)
        FillBuffer(ref_ratio,
            grid.get_global_startindex(), self.global_startindex,
            c_fields, g_fields,
            self.ActiveDimensions, grid.ActiveDimensions,
            grid.child_mask, self.domain_width, dls[grid.Level],
            self.axis)
