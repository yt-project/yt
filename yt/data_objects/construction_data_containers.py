"""
Data containers that require processing before they can be utilized.




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
import math
import weakref
import itertools
import shelve
from functools import wraps
import fileinput
from re import finditer
import os.path

from yt.config import ytcfg
from yt.funcs import *
from yt.utilities.logger import ytLogger
from .data_containers import \
    YTSelectionContainer1D, YTSelectionContainer2D, YTSelectionContainer3D, \
    restore_field_information_state, YTFieldData
from yt.utilities.lib.QuadTree import \
    QuadTree
from yt.utilities.lib.Interpolators import \
    ghost_zone_interpolate
from yt.utilities.lib.misc_utilities import \
    fill_region
from yt.utilities.lib.marching_cubes import \
    march_cubes_grid, march_cubes_grid_flux
from yt.utilities.data_point_utilities import CombineGrids,\
    DataCubeRefine, DataCubeReplace, FillRegion, FillBuffer
from yt.utilities.minimal_representation import \
    MinimalProjectionData
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    parallel_objects, parallel_root_only, communication_system
from yt.units.unit_object import Unit
import yt.geometry.particle_deposit as particle_deposit
from yt.utilities.grid_data_format.writer import write_to_gdf
from yt.frontends.stream.api import load_uniform_grid

from yt.fields.field_exceptions import \
    NeedsGridType,\
    NeedsOriginalGrid,\
    NeedsDataField,\
    NeedsProperty,\
    NeedsParameter
from yt.fields.derived_field import \
    TranslationFunc

class YTStreamlineBase(YTSelectionContainer1D):
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
    ds : dataset object
        Passed in to access the index
    kwargs : dict of items
        Any additional values are passed as field parameters that can be
        accessed by generated fields.

    Examples
    --------

    >>> from yt.visualization.api import Streamlines
    >>> streamlines = Streamlines(ds, [0.5]*3)
    >>> streamlines.integrate_through_volume()
    >>> stream = streamlines.path(0)
    >>> matplotlib.pylab.semilogy(stream['t'], stream['Density'], '-x')

    """
    _type_name = "streamline"
    _con_args = ('positions')
    sort_by = 't'
    def __init__(self, positions, length = 1.0, fields=None, ds=None, **kwargs):
        YTSelectionContainer1D.__init__(self, ds, fields, **kwargs)
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
        LE = self.ds.grid_left_edge
        RE = self.ds.grid_right_edge
        # Check left faces first
        min_streampoint = np.min(self.positions, axis=0)
        max_streampoint = np.max(self.positions, axis=0)
        p = np.all((min_streampoint <= RE) & (max_streampoint > LE), axis=1)
        self._grids = self.index.grids[p]

    def _get_data_from_grid(self, grid, field):
        # No child masking here; it happens inside the mask cut
        mask = self._get_cut_mask(grid)
        if field == 'dts': return self._dts[grid.id]
        if field == 't': return self._ts[grid.id]
        return grid[field].flat[mask]


    def _get_cut_mask(self, grid):
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
    """
    This is a data object corresponding to a line integral through the
    simulation domain.

    This object is typically accessed through the `proj` object that
    hangs off of index objects.  YTQuadTreeProj is a projection of a
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
    center : array_like, optional
        The 'center' supplied to fields that use it.  Note that this does
        not have to have `coord` as one value.  Strictly optional.
    data_source : `yt.data_objects.data_containers.YTSelectionContainer`, optional
        If specified, this will be the data source used for selecting
        regions to project.
    method : string, optional
        The method of projection to be performed.
        "integrate" : integration along the axis
        "mip" : maximum intensity projection
        "sum" : same as "integrate", except that we don't multiply by the path length
        WARNING: The "sum" option should only be used for uniform resolution grid
        datasets, as other datasets may result in unphysical images.
    style : string, optional
        The same as the method keyword.  Deprecated as of version 3.0.2.
        Please use method keyword instead.
    field_parameters : dict of items
        Values to be passed as field parameters that can be
        accessed by generated fields.

    Examples
    --------

    >>> ds = load("RedshiftOutput0005")
    >>> prj = ds.proj("density", 0)
    >>> print proj["density"]
    """
    _key_fields = YTSelectionContainer2D._key_fields + ['weight_field']
    _type_name = "proj"
    _con_args = ('axis', 'field', 'weight_field')
    _container_fields = ('px', 'py', 'pdx', 'pdy', 'weight_field')
    def __init__(self, field, axis, weight_field = None,
                 center = None, ds = None, data_source = None,
                 style = None, method = "integrate",
                 field_parameters = None):
        YTSelectionContainer2D.__init__(self, axis, ds, field_parameters)
        # Style is deprecated, but if it is set, then it trumps method
        # keyword.  TODO: Remove this keyword and this check at some point in
        # the future.
        if style is not None:
            method = style
        if method == "sum":
            self.method = "integrate"
            self._sum_only = True
        else:
            self.method = method
            self._sum_only = False
        if self.method == "mip":
            self.func = np.max
        elif self.method == "integrate":
            self.func = np.sum # for the future
        else:
            raise NotImplementedError(self.method)
        self._set_center(center)
        self._projected_units = {}
        if data_source is None: data_source = self.ds.all_data()
        for k, v in data_source.field_parameters.items():
            if k not in self.field_parameters or \
              self._is_default_field_parameter(k):
                self.set_field_parameter(k, v)
        self.data_source = data_source
        if weight_field is None:
            self.weight_field = weight_field
        else:
            self.weight_field = self._determine_fields(weight_field)[0]

        field = field or []
        field = self._determine_fields(ensure_list(field))

        if not self.deserialize(field):
            self.get_data(field)
            self.serialize()

    @property
    def blocks(self):
        return self.data_source.blocks

    @property
    def field(self):
        return [k for k in self.field_data.keys() if k not in self._container_fields]

    @property
    def _mrep(self):
        return MinimalProjectionData(self)

    def hub_upload(self):
        self._mrep.upload()

    def deserialize(self, fields):
        if not ytcfg.get("yt", "serialize"):
            return False

        deserialized_successfully = False
        store_file = self.ds.parameter_filename + '.yt'
        if os.path.isfile(store_file):
            deserialized_successfully = self._mrep.restore(store_file, self.ds)

            if deserialized_successfully:
                mylog.info("Using previous projection data from %s" % store_file)
                for field, field_data in self._mrep.field_data.items():
                    self[field] = field_data
        return deserialized_successfully

    def serialize(self):
        if not ytcfg.get("yt", "serialize"):
            return
        self._mrep.store(self.ds.parameter_filename + '.yt')

    def _get_tree(self, nvals):
        xax = self.ds.coordinates.x_axis[self.axis]
        yax = self.ds.coordinates.y_axis[self.axis]
        xd = self.ds.domain_dimensions[xax]
        yd = self.ds.domain_dimensions[yax]
        bounds = (self.ds.domain_left_edge[xax],
                  self.ds.domain_right_edge[yax],
                  self.ds.domain_left_edge[xax],
                  self.ds.domain_right_edge[yax])
        return QuadTree(np.array([xd,yd], dtype='int64'), nvals,
                        bounds, method = self.method)

    def get_data(self, fields = None):
        fields = fields or []
        fields = self._determine_fields(ensure_list(fields))
        # We need a new tree for every single set of fields we add
        if len(fields) == 0: return
        chunk_fields = fields[:]
        if self.weight_field is not None:
            chunk_fields.append(self.weight_field)
        tree = self._get_tree(len(fields))
        # This only needs to be done if we are in parallel; otherwise, we can
        # safely build the mesh as we go.
        if communication_system.communicators[-1].size > 1:
            for chunk in self.data_source.chunks([], "io", local_only = False):
                self._initialize_chunk(chunk, tree)
        self._initialize_projected_units(fields)
        with self.data_source._field_parameter_state(self.field_parameters):
            for chunk in parallel_objects(self.data_source.chunks(
                                          chunk_fields, "io", local_only = True)):
                mylog.debug("Adding chunk (%s) to tree (%0.3e GB RAM)", chunk.ires.size,
                    get_memory_usage()/1024.)
                self._handle_chunk(chunk, fields, tree)
        # Note that this will briefly double RAM usage
        if self.method == "mip":
            merge_style = -1
            op = "max"
        elif self.method == "integrate":
            merge_style = 1
            op = "sum"
        else:
            raise NotImplementedError
        # TODO: Add the combine operation
        xax = self.ds.coordinates.x_axis[self.axis]
        yax = self.ds.coordinates.y_axis[self.axis]
        ox = self.ds.domain_left_edge[xax]
        oy = self.ds.domain_left_edge[yax]
        px, py, pdx, pdy, nvals, nwvals = tree.get_all(False, merge_style)
        nvals = self.comm.mpi_allreduce(nvals, op=op)
        nwvals = self.comm.mpi_allreduce(nwvals, op=op)
        np.multiply(px, self.ds.domain_width[xax], px)
        np.add(px, ox, px)
        np.multiply(pdx, self.ds.domain_width[xax], pdx)

        np.multiply(py, self.ds.domain_width[yax], py)
        np.add(py, oy, py)
        np.multiply(pdy, self.ds.domain_width[yax], pdy)
        if self.weight_field is not None:
            # If there are 0s remaining in the weight vals
            # this will not throw an error, but silently
            # return nans for vals where dividing by 0
            # Leave as NaNs to be auto-masked by Matplotlib
            with np.errstate(invalid='ignore'):
                np.divide(nvals, nwvals[:,None], nvals)
        # We now convert to half-widths and center-points
        data = {}
        code_length = self.ds.domain_width.units
        data['px'] = self.ds.arr(px, code_length)
        data['py'] = self.ds.arr(py, code_length)
        data['weight_field'] = nwvals
        data['pdx'] = self.ds.arr(pdx, code_length)
        data['pdy'] = self.ds.arr(pdy, code_length)
        data['fields'] = nvals
        # Now we run the finalizer, which is ignored if we don't need it
        fd = data['fields']
        field_data = np.hsplit(data.pop('fields'), len(fields))
        for fi, field in enumerate(fields):
            finfo = self.ds._get_field_info(*field)
            mylog.debug("Setting field %s", field)
            input_units = self._projected_units[field]
            self[field] = self.ds.arr(field_data[fi].ravel(), input_units)
        for i in data.keys(): self[i] = data.pop(i)
        mylog.info("Projection completed")

    def _initialize_chunk(self, chunk, tree):
        icoords = chunk.icoords
        xax = self.ds.coordinates.x_axis[self.axis]
        yax = self.ds.coordinates.y_axis[self.axis]
        i1 = icoords[:,xax]
        i2 = icoords[:,yax]
        ilevel = chunk.ires * self.ds.ires_factor
        tree.initialize_chunk(i1, i2, ilevel)

    def _initialize_projected_units(self, fields):
        for field in fields:
            field_unit = Unit(self.ds.field_info[field].units,
                              registry=self.ds.unit_registry)
            if self.method == "mip" or self._sum_only:
                path_length_unit = Unit(registry=self.ds.unit_registry)
            else:
                ax_name = self.ds.coordinates.axis_name[self.axis]
                path_element_name = ("index", "path_element_%s" % (ax_name))
                path_length_unit = self.ds.field_info[path_element_name].units
                path_length_unit = Unit(path_length_unit,
                                        registry=self.ds.unit_registry)
                # Only convert to CGS for path elements that aren't angles
                if not path_length_unit.is_dimensionless:
                    path_length_unit = path_length_unit.get_cgs_equivalent()
            if self.weight_field is None:
                self._projected_units[field] = field_unit*path_length_unit
            else:
                self._projected_units[field] = field_unit

    def _handle_chunk(self, chunk, fields, tree):
        if self.method == "mip" or self._sum_only:
            dl = self.ds.quan(1.0, "")
        else:
            # This gets explicitly converted to cm
            ax_name = self.ds.coordinates.axis_name[self.axis]
            dl = chunk["index", "path_element_%s" % (ax_name)]
            # This is done for cases where our path element does not have a CGS
            # equivalent.  Once "preferred units" have been implemented, this
            # will not be necessary at all, as the final conversion will occur
            # at the display layer.
            if not dl.units.is_dimensionless:
                dl.convert_to_units("cm")
        v = np.empty((chunk.ires.size, len(fields)), dtype="float64")
        for i, field in enumerate(fields):
            d = chunk[field] * dl
            v[:,i] = d
        if self.weight_field is not None:
            w = chunk[self.weight_field]
            np.multiply(v, w[:,None], v)
            np.multiply(w, dl, w)
        else:
            w = np.ones(chunk.ires.size, dtype="float64")
        icoords = chunk.icoords
        xax = self.ds.coordinates.x_axis[self.axis]
        yax = self.ds.coordinates.y_axis[self.axis]
        i1 = icoords[:,xax]
        i2 = icoords[:,yax]
        ilevel = chunk.ires * self.ds.ires_factor
        tree.add_chunk_to_tree(i1, i2, ilevel, v, w)

    def to_pw(self, fields=None, center='c', width=None, origin='center-window'):
        r"""Create a :class:`~yt.visualization.plot_window.PWViewerMPL` from this
        object.

        This is a bare-bones mechanism of creating a plot window from this
        object, which can then be moved around, zoomed, and on and on.  All
        behavior of the plot window is relegated to that routine.
        """
        pw = self._get_pw(fields, center, width, origin, 'Projection')
        return pw

class YTCoveringGridBase(YTSelectionContainer3D):
    """A 3D region with all data extracted to a single, specified
    resolution.  Left edge should align with a cell boundary, but
    defaults to the closest cell boundary.

    Parameters
    ----------
    level : int
        The resolution level data to which data will be gridded. Level
        0 is the root grid dx for that dataset.
    left_edge : array_like
        The left edge of the region to be extracted.  Specify units by supplying
        a YTArray, otherwise code length units are assumed.
    dims : array_like
        Number of cells along each axis of resulting covering_grid
    fields : array_like, optional
        A list of fields that you'd like pre-generated for your object

    Examples
    --------
    >>> cube = ds.covering_grid(2, left_edge=[0.0, 0.0, 0.0], \
    ...                          dims=[128, 128, 128])
    """
    _spatial = True
    _type_name = "covering_grid"
    _con_args = ('level', 'left_edge', 'ActiveDimensions')
    _container_fields = (("index", "dx"),
                         ("index", "dy"),
                         ("index", "dz"),
                         ("index", "x"),
                         ("index", "y"),
                         ("index", "z"))
    _base_grid = None
    def __init__(self, level, left_edge, dims, fields = None,
                 ds = None, num_ghost_zones = 0, use_pbar = True,
                 field_parameters = None):
        if field_parameters is None:
            center = None
        else:
            center = field_parameters.get("center", None)
        YTSelectionContainer3D.__init__(self,
            center, ds, field_parameters)

        self.level = level

        if not iterable(left_edge):
            left_edge = [left_edge]*self.ds.dimensionality
        if len(left_edge) != self.ds.dimensionality:
            raise RuntimeError(
                "Length of left_edge must match the dimensionality of the "
                "dataset")
        if hasattr(left_edge, 'units'):
            le_units = left_edge.units
        else:
            le_units = 'code_length'
        self.left_edge = self.ds.arr(left_edge, le_units)

        if not iterable(dims):
            dims = [dims]*self.ds.dimensionality
        if len(dims) != self.ds.dimensionality:
            raise RuntimeError(
                "Length of dims must match the dimensionality of the dataset")
        self.ActiveDimensions = np.array(dims, dtype='int32')

        rdx = self.ds.domain_dimensions*self.ds.relative_refinement(0, level)
        rdx[np.where(np.array(dims) - 2 * num_ghost_zones <= 1)] = 1   # issue 602
        self.base_dds = self.ds.domain_width / self.ds.domain_dimensions
        self.dds = self.ds.domain_width / rdx.astype("float64")
        self.right_edge = self.left_edge + self.ActiveDimensions*self.dds
        self._num_ghost_zones = num_ghost_zones
        self._use_pbar = use_pbar
        self.global_startindex = np.rint(
            (self.left_edge-self.ds.domain_left_edge)/self.dds).astype('int64')
        self.domain_width = np.rint((self.ds.domain_right_edge -
                    self.ds.domain_left_edge)/self.dds).astype('int64')
        self._setup_data_source()
        self.get_data(fields)

    @property
    def icoords(self):
        ic = np.indices(self.ActiveDimensions).astype("int64")
        return np.column_stack([i.ravel() + gi for i, gi in
            zip(ic, self.get_global_startindex())])

    @property
    def fwidth(self):
        fw = np.ones((self.ActiveDimensions.prod(), 3), dtype="float64")
        fw *= self.dds
        return fw

    @property
    def fcoords(self):
        LE = self.LeftEdge + self.dds/2.0
        RE = self.RightEdge - self.dds/2.0
        N = self.ActiveDimensions
        fc = np.mgrid[LE[0]:RE[0]:N[0]*1j,
                      LE[1]:RE[1]:N[1]*1j,
                      LE[2]:RE[2]:N[2]*1j]
        return np.column_stack([f.ravel() for f in fc])

    @property
    def ires(self):
        tr = np.ones(self.ActiveDimensions.prod(), dtype="int64")
        tr *= self.level
        return tr

    def _reshape_vals(self, arr):
        if len(arr.shape) == 3: return arr
        return arr.reshape(self.ActiveDimensions, order="C")

    @property
    def shape(self):
        return tuple(self.ActiveDimensions.tolist())

    def _setup_data_source(self):
        LE = self.left_edge - self.base_dds
        RE = self.right_edge + self.base_dds
        if not all(self.ds.periodicity):
            for i in range(3):
                if self.ds.periodicity[i]: continue
                LE[i] = max(LE[i], self.ds.domain_left_edge[i])
                RE[i] = min(RE[i], self.ds.domain_right_edge[i])
        self._data_source = self.ds.region(self.center, LE, RE)
        self._data_source.min_level = 0
        self._data_source.max_level = self.level
        self._pdata_source = self.ds.region(self.center,
            self.left_edge, self.right_edge)
        self._pdata_source.min_level = 0
        self._pdata_source.max_level = self.level

    def get_data(self, fields = None):
        if fields is None: return
        fields = self._determine_fields(ensure_list(fields))
        fields_to_get = [f for f in fields if f not in self.field_data]
        fields_to_get = self._identify_dependencies(fields_to_get)
        if len(fields_to_get) == 0: return
        fill, gen, part, alias = self._split_fields(fields_to_get)
        if len(part) > 0: self._fill_particles(part)
        if len(fill) > 0: self._fill_fields(fill)
        for a, f in sorted(alias.items()):
            self[a] = f(self)
            self.field_data[a].convert_to_units(f.output_units)
        if len(gen) > 0: self._generate_fields(gen)

    def _split_fields(self, fields_to_get):
        fill, gen = self.index._split_fields(fields_to_get)
        particles = []
        alias = {}
        for field in gen:
            if field[0] == 'deposit':
                fill.append(field)
                continue
            finfo = self.ds._get_field_info(*field)
            if finfo._function.func_name == "_TranslationFunc":
                alias[field] = finfo
                continue
            try:
                finfo.check_available(self)
            except NeedsOriginalGrid:
                fill.append(field)
        for field in fill:
            finfo = self.ds._get_field_info(*field)
            if finfo.particle_type:
                particles.append(field)
        gen = [f for f in gen if f not in fill and f not in alias]
        fill = [f for f in fill if f not in particles]
        return fill, gen, particles, alias

    def _fill_particles(self, part):
        for p in part:
            self[p] = self._pdata_source[p]

    def _fill_fields(self, fields):
        fields = [f for f in fields if f not in self.field_data]
        if len(fields) == 0: return
        output_fields = [np.zeros(self.ActiveDimensions, dtype="float64")
                         for field in fields]
        domain_dims = self.ds.domain_dimensions.astype("int64") \
                    * self.ds.relative_refinement(0, self.level)
        for chunk in self._data_source.chunks(fields, "io"):
            input_fields = [chunk[field] for field in fields]
            # NOTE: This usage of "refine_by" is actually *okay*, because it's
            # being used with respect to iref, which is *already* scaled!
            fill_region(input_fields, output_fields, self.level,
                        self.global_startindex, chunk.icoords, chunk.ires,
                        domain_dims, self.ds.refine_by)
        for name, v in zip(fields, output_fields):
            fi = self.ds._get_field_info(*name)
            self[name] = self.ds.arr(v, fi.units)

    def _generate_container_field(self, field):
        rv = self.ds.arr(np.ones(self.ActiveDimensions, dtype="float64"),
                             "")
        if field == ("index", "dx"):
            np.multiply(rv, self.dds[0], rv)
        elif field == ("index", "dy"):
            np.multiply(rv, self.dds[1], rv)
        elif field == ("index", "dz"):
            np.multiply(rv, self.dds[2], rv)
        elif field == ("index", "x"):
            x = np.mgrid[self.left_edge[0] + 0.5*self.dds[0]:
                         self.right_edge[0] - 0.5*self.dds[0]:
                         self.ActiveDimensions[0] * 1j]
            np.multiply(rv, x[:,None,None], rv)
        elif field == ("index", "y"):
            y = np.mgrid[self.left_edge[1] + 0.5*self.dds[1]:
                         self.right_edge[1] - 0.5*self.dds[1]:
                         self.ActiveDimensions[1] * 1j]
            np.multiply(rv, y[None,:,None], rv)
        elif field == ("index", "z"):
            z = np.mgrid[self.left_edge[2] + 0.5*self.dds[2]:
                         self.right_edge[2] - 0.5*self.dds[2]:
                         self.ActiveDimensions[2] * 1j]
            np.multiply(rv, z[None,None,:], rv)
        else:
            raise KeyError(field)
        return rv

    @property
    def LeftEdge(self):
        return self.left_edge

    @property
    def RightEdge(self):
        return self.right_edge

    def deposit(self, positions, fields = None, method = None):
        cls = getattr(particle_deposit, "deposit_%s" % method, None)
        if cls is None:
            raise YTParticleDepositionNotImplemented(method)
        op = cls(self.ActiveDimensions.prod()) # We allocate number of zones, not number of octs
        op.initialize()
        op.process_grid(self, positions, fields)
        vals = op.finalize()
        return vals.reshape(self.ActiveDimensions, order="C")

    def write_to_gdf(self, gdf_path, fields, nprocs=1, field_units=None,
                     **kwargs):
        r"""
        Write the covering grid data to a GDF file.

        Parameters
        ----------
        gdf_path : string
            Pathname of the GDF file to write.
        fields : list of strings
            Fields to write to the GDF file.
        nprocs : integer, optional
            Split the covering grid into *nprocs* subgrids before
            writing to the GDF file. Default: 1
        field_units : dictionary, optional
            Dictionary of units to convert fields to. If not set, fields are
            in their default units.
        All remaining keyword arguments are passed to
        yt.utilities.grid_data_format.writer.write_to_gdf.

        Examples
        --------
        >>> cube.write_to_gdf("clumps.h5", ["density","temperature"], nprocs=16,
        ...                   clobber=True)
        """
        data = {}
        for field in fields:
            if field in field_units:
                units = field_units[field]
            else:
                units = str(self[field].units)
            data[field] = (self[field].in_units(units).v, units)
        le = self.left_edge.v
        re = self.right_edge.v
        bbox = np.array([[l,r] for l,r in zip(le, re)])
        ds = load_uniform_grid(data, self.ActiveDimensions, bbox=bbox,
                               length_unit=self.ds.length_unit,
                               time_unit=self.ds.time_unit,
                               mass_unit=self.ds.mass_unit, nprocs=nprocs,
                               sim_time=self.ds.current_time.v)
        write_to_gdf(ds, gdf_path, **kwargs)

class YTArbitraryGridBase(YTCoveringGridBase):
    """A 3D region with arbitrary bounds and dimensions.

    In contrast to the Covering Grid, this object accepts a left edge, a right
    edge, and dimensions.  This allows it to be used for creating 3D particle
    deposition fields that are independent of the underlying mesh, whether that
    is yt-generated or from the simulation data.  For example, arbitrary boxes
    around particles can be drawn and particle deposition fields can be
    created.  This object will refuse to generate any fluid fields.

    Parameters
    ----------
    left_edge : array_like
        The left edge of the region to be extracted
    right_edge : array_like
        The left edge of the region to be extracted
    dims : array_like
        Number of cells along each axis of resulting grid.

    Examples
    --------
    >>> obj = ds.arbitrary_grid([0.0, 0.0, 0.0], [0.99, 0.99, 0.99],
    ...                          dims=[128, 128, 128])
    """
    _spatial = True
    _type_name = "arbitrary_grid"
    _con_args = ('left_edge', 'right_edge', 'ActiveDimensions')
    _container_fields = (("index", "dx"),
                         ("index", "dy"),
                         ("index", "dz"),
                         ("index", "x"),
                         ("index", "y"),
                         ("index", "z"))
    def __init__(self, left_edge, right_edge, dims,
                 ds = None, field_parameters = None):
        if field_parameters is None:
            center = None
        else:
            center = field_parameters.get("center", None)
        YTSelectionContainer3D.__init__(self, center, ds, field_parameters)
        self.left_edge = np.array(left_edge)
        self.right_edge = np.array(right_edge)
        self.ActiveDimensions = np.array(dims, dtype='int32')
        if self.ActiveDimensions.size == 1:
            self.ActiveDimensions = np.array([dims, dims, dims], dtype="int32")
        self.dds = self.base_dds = (self.right_edge - self.left_edge)/self.ActiveDimensions
        self.level = 99
        self._setup_data_source()

    def _fill_fields(self, fields):
        raise NotImplementedError

class LevelState(object):
    current_dx = None
    current_dims = None
    current_level = None
    global_startindex = None
    old_global_startindex = None
    domain_iwidth = None
    fields = None
    data_source = None

    # These are all cached here as numpy arrays, without units, in
    # code_lengths.
    domain_width = None
    domain_left_edge = None
    domain_right_edge = None
    left_edge = None
    right_edge = None
    base_dx = None
    dds = None

class YTSmoothedCoveringGridBase(YTCoveringGridBase):
    """A 3D region with all data extracted and interpolated to a
    single, specified resolution. (Identical to covering_grid,
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
        Number of cells along each axis of resulting covering_grid.
    fields : array_like, optional
        A list of fields that you'd like pre-generated for your object

    Example
    -------
    cube = ds.smoothed_covering_grid(2, left_edge=[0.0, 0.0, 0.0], \
                              dims=[128, 128, 128])
    """
    _type_name = "smoothed_covering_grid"
    filename = None
    @wraps(YTCoveringGridBase.__init__)
    def __init__(self, *args, **kwargs):
        self._base_dx = (
              (self.ds.domain_right_edge - self.ds.domain_left_edge) /
               self.ds.domain_dimensions.astype("float64"))
        self.global_endindex = None
        YTCoveringGridBase.__init__(self, *args, **kwargs)
        self._final_start_index = self.global_startindex

    def _setup_data_source(self, level_state = None):
        if level_state is None: return
        # We need a buffer region to allow for zones that contribute to the
        # interpolation but are not directly inside our bounds
        level_state.data_source = self.ds.region(
            self.center,
            level_state.left_edge - level_state.current_dx,
            level_state.right_edge + level_state.current_dx)
        level_state.data_source.min_level = level_state.current_level
        level_state.data_source.max_level = level_state.current_level
        self._pdata_source = self.ds.region(
            self.center,
            level_state.left_edge - level_state.current_dx,
            level_state.right_edge + level_state.current_dx)
        self._pdata_source.min_level = level_state.current_level
        self._pdata_source.max_level = level_state.current_level


    def _fill_fields(self, fields):
        fields = [f for f in fields if f not in self.field_data]
        if len(fields) == 0: return
        ls = self._initialize_level_state(fields)
        for level in range(self.level + 1):
            domain_dims = self.ds.domain_dimensions.astype("int64") \
                        * self.ds.relative_refinement(0, self.level)
            for chunk in ls.data_source.chunks(fields, "io"):
                chunk[fields[0]]
                input_fields = [chunk[field] for field in fields]
                # NOTE: This usage of "refine_by" is actually *okay*, because it's
                # being used with respect to iref, which is *already* scaled!
                fill_region(input_fields, ls.fields, ls.current_level,
                            ls.global_startindex, chunk.icoords,
                            chunk.ires, domain_dims, self.ds.refine_by)
            self._update_level_state(ls)
        for name, v in zip(fields, ls.fields):
            if self.level > 0: v = v[1:-1,1:-1,1:-1]
            fi = self.ds._get_field_info(*name)
            self[name] = self.ds.arr(v, fi.units)

    def _initialize_level_state(self, fields):
        ls = LevelState()
        ls.domain_width = self.ds.domain_width
        ls.domain_left_edge = self.ds.domain_left_edge
        ls.domain_right_edge = self.ds.domain_right_edge
        ls.left_edge = self.left_edge
        ls.right_edge = self.right_edge
        ls.base_dx = self._base_dx
        ls.dds = self.dds
        for att in ("domain_width", "domain_left_edge", "domain_right_edge",
                    "left_edge", "right_edge", "base_dx", "dds"):
            setattr(ls, att, getattr(ls, att).in_units("code_length").d)
        ls.current_dx = ls.base_dx
        ls.current_level = 0
        LL = self.left_edge - self.ds.domain_left_edge
        ls.global_startindex = np.rint(LL / ls.current_dx).astype('int64') - 1
        ls.domain_iwidth = np.rint((self.ds.domain_right_edge -
                    self.ds.domain_left_edge)/ls.current_dx).astype('int64')
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
        rf = float(self.ds.relative_refinement(
                    ls.current_level, ls.current_level + 1))
        ls.current_level += 1
        ls.current_dx = ls.base_dx / \
            self.ds.relative_refinement(0, ls.current_level)
        self._setup_data_source(ls)
        LL = ls.left_edge - ls.domain_left_edge
        ls.old_global_startindex = ls.global_startindex
        ls.global_startindex = np.rint(LL / ls.current_dx).astype('int64') - 1
        ls.domain_iwidth = np.rint(ls.domain_width/ls.current_dx).astype('int64')
        input_left = (level_state.old_global_startindex + 0.5) * rf
        width = (self.ActiveDimensions*ls.dds)
        output_dims = np.rint(width/level_state.current_dx+0.5).astype("int32") + 2
        level_state.current_dims = output_dims
        new_fields = []
        for input_field in level_state.fields:
            output_field = np.zeros(output_dims, dtype="float64")
            output_left = level_state.global_startindex + 0.5
            ghost_zone_interpolate(rf, input_field, input_left,
                                   output_field, output_left)
            new_fields.append(output_field)
        level_state.fields = new_fields

class YTSurfaceBase(YTSelectionContainer3D):
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
    data_source : YTSelectionContainer
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

    >>> sp = ds.sphere("max", (10, "kpc")
    >>> surf = ds.surface(sp, "Density", 5e-27)
    >>> print surf["Temperature"]
    >>> print surf.vertices
    >>> bounds = [(sp.center[i] - 5.0/ds['kpc'],
    ...            sp.center[i] + 5.0/ds['kpc']) for i in range(3)]
    >>> surf.export_ply("my_galaxy.ply", bounds = bounds)
    """
    _type_name = "surface"
    _con_args = ("data_source", "surface_field", "field_value")
    _container_fields = (("index", "dx"),
                         ("index", "dy"),
                         ("index", "dz"),
                         ("index", "x"),
                         ("index", "y"),
                         ("index", "z"))
    vertices = None
    def __init__(self, data_source, surface_field, field_value):
        self.data_source = data_source
        self.surface_field = surface_field
        self.field_value = field_value
        self.vertex_samples = YTFieldData()
        center = data_source.get_field_parameter("center")
        super(YTSurfaceBase, self).__init__(center = center, ds =
                    data_source.ds )

    def _generate_container_field(self, field):
        self.get_data(field)
        return self[field]

    def get_data(self, fields = None, sample_type = "face"):
        if isinstance(fields, list) and len(fields) > 1:
            for field in fields: self.get_data(field)
            return
        elif isinstance(fields, list):
            fields = fields[0]
        # Now we have a "fields" value that is either a string or None
        mylog.info("Extracting (sampling: %s)" % (fields,))
        verts = []
        samples = []
        deps = self._determine_fields(self.surface_field)
        deps = self._identify_dependencies(deps, spatial=True)
        for io_chunk in parallel_objects(self.data_source.chunks(deps, "io",
                                         preload_fields = deps)):
            for block, mask in self.data_source.blocks:
                my_verts = self._extract_isocontours_from_grid(
                                block, self.surface_field, self.field_value,
                                mask, fields, sample_type)
                if fields is not None:
                    my_verts, svals = my_verts
                    samples.append(svals)
                verts.append(my_verts)
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
                                       mask, sample_values = None,
                                       sample_type = "face"):
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

        >>> sp = ds.sphere("max", (10, "kpc")
        >>> surf = ds.surface(sp, "Density", 5e-27)
        >>> flux = surf.calculate_flux(
        ...     "velocity_x", "velocity_y", "velocity_z", "Metal_Density")
        """
        flux = 0.0
        mylog.info("Fluxing %s", fluxing_field)
        deps = [field_x, field_y, field_z]
        if fluxing_field is not None: deps.append(fluxing_field)
        deps = self._determine_fields(deps)
        deps = self._identify_dependencies(deps)
        for io_chunk in parallel_objects(self.data_source.chunks(deps, "io",
                                preload_fields = deps)):
            for block, mask in self.data_source.blocks:
                flux += self._calculate_flux_in_grid(block, mask,
                        field_x, field_y, field_z, fluxing_field)
        flux = self.comm.mpi_allreduce(flux, op="sum")
        return flux

    def _calculate_flux_in_grid(self, grid, mask,
            field_x, field_y, field_z, fluxing_field = None):
        vals = grid.get_vertex_centered_data(self.surface_field)
        if fluxing_field is None:
            ff = np.ones(vals.shape, dtype="float64")
        else:
            ff = grid.get_vertex_centered_data(fluxing_field)
        xv, yv, zv = [grid.get_vertex_centered_data(f)
                      for f in [field_x, field_y, field_z]]
        return march_cubes_grid_flux(self.field_value, vals, xv, yv, zv,
                    ff, mask, grid.LeftEdge, grid.dds)

    @property
    def triangles(self):
        if self.vertices is None:
            self.get_data()
        vv = np.empty((self.vertices.shape[1]/3, 3, 3), dtype="float64")
        for i in range(3):
            for j in range(3):
                vv[:,i,j] = self.vertices[j,i::3]
        return vv

    def export_obj(self, filename, transparency = 1.0, dist_fac = None,
                   color_field = None, emit_field = None, color_map = "algae",
                   color_log = True, emit_log = True, plot_index = None,
                   color_field_max = None, color_field_min = None,
                   emit_field_max = None, emit_field_min = None):
        r"""This exports the surface to the OBJ format, suitable for visualization
        in many different programs (e.g., Blender).  NOTE: this exports an .obj file
        and an .mtl file, both with the general 'filename' as a prefix.
        The .obj file points to the .mtl file in its header, so if you move the 2
        files, make sure you change the .obj header to account for this. ALSO NOTE:
        the emit_field needs to be a combination of the other 2 fields used to
        have the emissivity track with the color.

        Parameters
         ----------
        filename : string
            The file this will be exported to.  This cannot be a file-like object.
            Note - there are no file extentions included - both obj & mtl files
            are created.
        transparency : float
            This gives the transparency of the output surface plot.  Values
            from 0.0 (invisible) to 1.0 (opaque).
        dist_fac : float
            Divide the axes distances by this amount.
        color_field : string
            Should a field be sample and colormapped?
        emit_field : string
            Should we track the emissivity of a field?
              NOTE: this should be a combination of the other 2 fields being used.
        color_map : string
            Which color map should be applied?
        color_log : bool
            Should the color field be logged before being mapped?
        emit_log : bool
            Should the emitting field be logged before being mapped?
        plot_index : integer
            Index of plot for multiple plots.  If none, then only 1 plot.
        color_field_max : float
            Maximum value of the color field across all surfaces.
        color_field_min : float
            Minimum value of the color field across all surfaces.
        emit_field_max : float
            Maximum value of the emitting field across all surfaces.
        emit_field_min : float
            Minimum value of the emitting field across all surfaces.

        Examples
        --------

        >>> sp = ds.sphere("max", (10, "kpc"))
        >>> trans = 1.0
        >>> distf = 3.1e18*1e3 # distances into kpc
        >>> surf = ds.surface(sp, "Density", 5e-27)
        >>> surf.export_obj("my_galaxy", transparency=trans, dist_fac = distf)

        >>> sp = ds.sphere("max", (10, "kpc"))
        >>> mi, ma = sp.quantities['Extrema']('Temperature')[0]
        >>> rhos = [1e-24, 1e-25]
        >>> trans = [0.5, 1.0]
        >>> distf = 3.1e18*1e3 # distances into kpc
        >>> for i, r in enumerate(rhos):
        ...     surf = ds.surface(sp,'Density',r)
        ...     surf.export_obj("my_galaxy", transparency=trans[i],
        ...                      color_field='Temperature', dist_fac = distf,
        ...                      plot_index = i, color_field_max = ma,
        ...                      color_field_min = mi)

        >>> sp = ds.sphere("max", (10, "kpc"))
        >>> rhos = [1e-24, 1e-25]
        >>> trans = [0.5, 1.0]
        >>> distf = 3.1e18*1e3 # distances into kpc
        >>> def _Emissivity(field, data):
        ...     return (data['Density']*data['Density']*np.sqrt(data['Temperature']))
        >>> add_field("Emissivity", function=_Emissivity, units=r"\rm{g K}/\rm{cm}^{6}")
        >>> for i, r in enumerate(rhos):
        ...     surf = ds.surface(sp,'Density',r)
        ...     surf.export_obj("my_galaxy", transparency=trans[i],
        ...                      color_field='Temperature', emit_field = 'Emissivity',
        ...                      dist_fac = distf, plot_index = i)

        """
        if self.vertices is None:
            if color_field is not None:
                self.get_data(color_field,"face")
        elif color_field is not None:
            if color_field not in self.field_data:
                self[color_field]
        if color_field is None:
            self.get_data(self.surface_field,'face')
        if emit_field is not None:
            if color_field not in self.field_data:
                self[emit_field]
        only_on_root(self._export_obj, filename, transparency, dist_fac, color_field, emit_field,
                             color_map, color_log, emit_log, plot_index, color_field_max,
                             color_field_min, emit_field_max, emit_field_min)

    def _color_samples_obj(self, cs, em, color_log, emit_log, color_map, arr,
                           color_field_max, color_field_min, color_field,
                           emit_field_max, emit_field_min, emit_field): # this now holds for obj files
        from sys import version
        if color_field is not None:
            if color_log: cs = np.log10(cs)
        if emit_field is not None:
            if emit_log: em = np.log10(em)
        if color_field is not None:
            if color_field_min is None:
                if version >= '3':
                    cs = [float(field) for field in cs]
                    cs = np.array(cs)
                mi = cs.min()
            else:
                mi = color_field_min
                if color_log: mi = np.log10(mi)
            if color_field_max is None:
                if version >= '3':
                    cs = [float(field) for field in cs]
                    cs = np.array(cs)
                ma = cs.max()
            else:
                ma = color_field_max
                if color_log: ma = np.log10(ma)
            cs = (cs - mi) / (ma - mi)
        else:
            cs[:] = 1.0
        # to get color indicies for OBJ formatting
        from yt.visualization._colormap_data import color_map_luts
        lut = color_map_luts[color_map]
        x = np.mgrid[0.0:1.0:lut[0].shape[0]*1j]
        arr["cind"][:] = (np.interp(cs,x,x)*(lut[0].shape[0]-1)).astype("uint8")
        # now, get emission
        if emit_field is not None:
            if emit_field_min is None:
                if version >= '3':
                    em = [float(field) for field in em]
                    em = np.array(em)
                emi = em.min()
            else:
                emi = emit_field_min
                if emit_log: emi = np.log10(emi)
            if emit_field_max is None:
                if version >= '3':
                    em = [float(field) for field in em]
                    em = np.array(em)
                ema = em.max()
            else:
                ema = emit_field_max
                if emit_log: ema = np.log10(ema)
            em = (em - emi)/(ema - emi)
            x = np.mgrid[0.0:255.0:2j] # assume 1 emissivity per color
            arr["emit"][:] = (np.interp(em,x,x))*2.0 # for some reason, max emiss = 2
        else:
            arr["emit"][:] = 0.0


    @parallel_root_only
    def _export_obj(self, filename, transparency, dist_fac = None,
                    color_field = None, emit_field = None, color_map = "algae",
                    color_log = True, emit_log = True, plot_index = None,
                    color_field_max = None, color_field_min = None,
                    emit_field_max = None, emit_field_min = None):
        from sys import version
        from io import IOBase
        if plot_index is None:
            plot_index = 0
        if version < '3':
            checker = file
        else:
            checker = IOBase
        if isinstance(filename, checker):
            fobj = filename + '.obj'
            fmtl = filename + '.mtl'
        else:
            if plot_index == 0:
                fobj = open(filename + '.obj', "w")
                fmtl = open(filename + '.mtl', 'w')
                cc = 1
            else:
                # read in last vertex
                linesave = ''
                for line in fileinput.input(filename + '.obj'):
                    if line[0] == 'f':
                        linesave = line
                p = [m.start() for m in finditer(' ', linesave)]
                cc = int(linesave[p[len(p)-1]:])+1
                fobj = open(filename + '.obj', "a")
                fmtl = open(filename + '.mtl', 'a')
        ftype = [("cind", "uint8"), ("emit", "float")]
        vtype = [("x","float"),("y","float"), ("z","float")]
        if plot_index == 0:
            fobj.write("# yt OBJ file\n")
            fobj.write("# www.yt-project.com\n")
            fobj.write("mtllib " + filename + '.mtl\n\n')  # use this material file for the faces
            fmtl.write("# yt MLT file\n")
            fmtl.write("# www.yt-project.com\n\n")
        #(0) formulate vertices
        nv = self.vertices.shape[1] # number of groups of vertices
        f = np.empty(nv/self.vertices.shape[0], dtype=ftype) # store sets of face colors
        v = np.empty(nv, dtype=vtype) # stores vertices
        if color_field is not None:
            cs = self[color_field]
        else:
            cs = np.empty(self.vertices.shape[1]/self.vertices.shape[0])
        if emit_field is not None:
            em = self[emit_field]
        else:
            em = np.empty(self.vertices.shape[1]/self.vertices.shape[0])
        self._color_samples_obj(cs, em, color_log, emit_log, color_map, f,
                                color_field_max, color_field_min,  color_field,
                                emit_field_max, emit_field_min, emit_field) # map color values to color scheme
        from yt.visualization._colormap_data import color_map_luts # import colors for mtl file
        lut = color_map_luts[color_map] # enumerate colors
        # interpolate emissivity to enumerated colors
        emiss = np.interp(np.mgrid[0:lut[0].shape[0]],np.mgrid[0:len(cs)],f["emit"][:])
        if dist_fac is None: # then normalize by bounds
            DLE = self.pf.domain_left_edge
            DRE = self.pf.domain_right_edge
            bounds = [(DLE[i], DRE[i]) for i in range(3)]
            for i, ax in enumerate("xyz"):
                # Do the bounds first since we cast to f32
                tmp = self.vertices[i,:]
                np.subtract(tmp, bounds[i][0], tmp)
                w = bounds[i][1] - bounds[i][0]
                np.divide(tmp, w, tmp)
                np.subtract(tmp, 0.5, tmp) # Center at origin.
                v[ax][:] = tmp
        else:
            for i, ax in enumerate("xyz"):
                tmp = self.vertices[i,:]
                np.divide(tmp, dist_fac, tmp)
                v[ax][:] = tmp
        #(1) write all colors per surface to mtl file
        for i in range(0,lut[0].shape[0]):
            omname = "material_" + str(i) + '_' + str(plot_index)  # name of the material
            fmtl.write("newmtl " + omname +'\n') # the specific material (color) for this face
            fmtl.write("Ka %.6f %.6f %.6f\n" %(0.0, 0.0, 0.0)) # ambient color, keep off
            fmtl.write("Kd %.6f %.6f %.6f\n" %(lut[0][i], lut[1][i], lut[2][i])) # color of face
            fmtl.write("Ks %.6f %.6f %.6f\n" %(0.0, 0.0, 0.0)) # specular color, keep off
            fmtl.write("d %.6f\n" %(transparency))  # transparency
            fmtl.write("em %.6f\n" %(emiss[i])) # emissivity per color
            fmtl.write("illum 2\n") # not relevant, 2 means highlights on?
            fmtl.write("Ns %.6f\n\n" %(0.0)) #keep off, some other specular thing
        #(2) write vertices
        for i in range(0,self.vertices.shape[1]):
            fobj.write("v %.6f %.6f %.6f\n" %(v["x"][i], v["y"][i], v["z"][i]))
        fobj.write("#done defining vertices\n\n")
        #(3) define faces and materials for each face
        for i in range(0,self.triangles.shape[0]):
            omname = 'material_' + str(f["cind"][i]) + '_' + str(plot_index) # which color to use
            fobj.write("usemtl " + omname + '\n') # which material to use for this face (color)
            fobj.write("f " + str(cc) + ' ' + str(cc+1) + ' ' + str(cc+2) + '\n\n') # vertices to color
            cc = cc+3
        fmtl.close()
        fobj.close()


    def export_blender(self,  transparency = 1.0, dist_fac = None,
                   color_field = None, emit_field = None, color_map = "algae",
                   color_log = True, emit_log = True, plot_index = None,
                   color_field_max = None, color_field_min = None,
                   emit_field_max = None, emit_field_min = None):
        r"""This exports the surface to the OBJ format, suitable for visualization
        in many different programs (e.g., Blender).  NOTE: this exports an .obj file
        and an .mtl file, both with the general 'filename' as a prefix.
        The .obj file points to the .mtl file in its header, so if you move the 2
        files, make sure you change the .obj header to account for this. ALSO NOTE:
        the emit_field needs to be a combination of the other 2 fields used to
        have the emissivity track with the color.

        Parameters
        ----------
        filename : string
            The file this will be exported to.  This cannot be a file-like object.
            Note - there are no file extentions included - both obj & mtl files
            are created.
        transparency : float
            This gives the transparency of the output surface plot.  Values
            from 0.0 (invisible) to 1.0 (opaque).
        dist_fac : float
            Divide the axes distances by this amount.
        color_field : string
            Should a field be sample and colormapped?
        emit_field : string
            Should we track the emissivity of a field?
              NOTE: this should be a combination of the other 2 fields being used.
        color_map : string
            Which color map should be applied?
        color_log : bool
            Should the color field be logged before being mapped?
        emit_log : bool
            Should the emitting field be logged before being mapped?
        plot_index : integer
            Index of plot for multiple plots.  If none, then only 1 plot.
        color_field_max : float
            Maximum value of the color field across all surfaces.
        color_field_min : float
            Minimum value of the color field across all surfaces.
        emit_field_max : float
            Maximum value of the emitting field across all surfaces.
        emit_field_min : float
            Minimum value of the emitting field across all surfaces.

        Examples
        --------

        >>> sp = ds.sphere("max", (10, "kpc"))
        >>> trans = 1.0
        >>> distf = 3.1e18*1e3 # distances into kpc
        >>> surf = ds.surface(sp, "Density", 5e-27)
        >>> surf.export_obj("my_galaxy", transparency=trans, dist_fac = distf)

        >>> sp = ds.sphere("max", (10, "kpc"))
        >>> mi, ma = sp.quantities['Extrema']('Temperature')[0]
        >>> rhos = [1e-24, 1e-25]
        >>> trans = [0.5, 1.0]
        >>> distf = 3.1e18*1e3 # distances into kpc
        >>> for i, r in enumerate(rhos):
        ...     surf = ds.surface(sp,'Density',r)
        ...     surf.export_obj("my_galaxy", transparency=trans[i],
        ...                      color_field='Temperature', dist_fac = distf,
        ...                      plot_index = i, color_field_max = ma,
        ...                      color_field_min = mi)

        >>> sp = ds.sphere("max", (10, "kpc"))
        >>> rhos = [1e-24, 1e-25]
        >>> trans = [0.5, 1.0]
        >>> distf = 3.1e18*1e3 # distances into kpc
        >>> def _Emissivity(field, data):
        ...     return (data['Density']*data['Density']*np.sqrt(data['Temperature']))
        >>> add_field("Emissivity", function=_Emissivity, units=r"\rm{g K}/\rm{cm}^{6}")
        >>> for i, r in enumerate(rhos):
        ...     surf = ds.surface(sp,'Density',r)
        ...     surf.export_obj("my_galaxy", transparency=trans[i],
        ...                      color_field='Temperature', emit_field = 'Emissivity',
        ...                      dist_fac = distf, plot_index = i)

        """
        if self.vertices is None:
            if color_field is not None:
                self.get_data(color_field,"face")
        elif color_field is not None:
            if color_field not in self.field_data:
                self[color_field]
        if color_field is None:
            self.get_data(self.surface_field,'face')
        if emit_field is not None:
            if color_field not in self.field_data:
                self[emit_field]
        fullverts, colors, alpha, emisses, colorindex = only_on_root(self._export_blender,
                                                                transparency, dist_fac, color_field, emit_field,
                                                                color_map, color_log, emit_log, plot_index,
                                                                color_field_max,
                                                                color_field_min, emit_field_max, emit_field_min)
        return fullverts, colors, alpha, emisses, colorindex

    def _export_blender(self, transparency, dist_fac = None,
                    color_field = None, emit_field = None, color_map = "algae",
                    color_log = True, emit_log = True, plot_index = None,
                    color_field_max = None, color_field_min = None,
                    emit_field_max = None, emit_field_min = None):
        import io
        from sys import version
        if plot_index is None:
            plot_index = 0
            vmax=0
        ftype = [("cind", "uint8"), ("emit", "float")]
        vtype = [("x","float"),("y","float"), ("z","float")]
        #(0) formulate vertices
        nv = self.vertices.shape[1] # number of groups of vertices
        f = np.empty(nv/self.vertices.shape[0], dtype=ftype) # store sets of face colors
        v = np.empty(nv, dtype=vtype) # stores vertices
        if color_field is not None:
            cs = self[color_field]
        else:
            cs = np.empty(self.vertices.shape[1]/self.vertices.shape[0])
        if emit_field is not None:
            em = self[emit_field]
        else:
            em = np.empty(self.vertices.shape[1]/self.vertices.shape[0])
        self._color_samples_obj(cs, em, color_log, emit_log, color_map, f,
                                color_field_max, color_field_min, color_field,
                                emit_field_max, emit_field_min, emit_field) # map color values to color scheme
        from yt.visualization._colormap_data import color_map_luts # import colors for mtl file
        lut = color_map_luts[color_map] # enumerate colors
        # interpolate emissivity to enumerated colors
        emiss = np.interp(np.mgrid[0:lut[0].shape[0]],np.mgrid[0:len(cs)],f["emit"][:])
        if dist_fac is None: # then normalize by bounds
            DLE = self.ds.domain_left_edge
            DRE = self.ds.domain_right_edge
            bounds = [(DLE[i], DRE[i]) for i in range(3)]
            for i, ax in enumerate("xyz"):
                # Do the bounds first since we cast to f32
                tmp = self.vertices[i,:]
                np.subtract(tmp, bounds[i][0], tmp)
                w = bounds[i][1] - bounds[i][0]
                np.divide(tmp, w, tmp)
                np.subtract(tmp, 0.5, tmp) # Center at origin.
                v[ax][:] = tmp
        else:
            for i, ax in enumerate("xyz"):
                tmp = self.vertices[i,:]
                np.divide(tmp, dist_fac, tmp)
                v[ax][:] = tmp
        return  v, lut, transparency, emiss, f['cind']


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

        >>> sp = ds.sphere("max", (10, "kpc")
        >>> surf = ds.surface(sp, "Density", 5e-27)
        >>> print surf["Temperature"]
        >>> print surf.vertices
        >>> bounds = [(sp.center[i] - 5.0/ds['kpc'],
        ...            sp.center[i] + 5.0/ds['kpc']) for i in range(3)]
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
            DLE = self.ds.domain_left_edge
            DRE = self.ds.domain_right_edge
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

        >>> from yt.mods import *
        >>> ds = load("redshift0058")
        >>> dd = ds.sphere("max", (200, "kpc"))
        >>> rho = 5e-27
        >>> bounds = [(dd.center[i] - 100.0/ds['kpc'],
        ...            dd.center[i] + 100.0/ds['kpc']) for i in range(3)]
        ...
        >>> surf = ds.surface(dd, "Density", rho)
        >>> rv = surf.export_sketchfab(
        ...     title = "Testing Upload",
        ...     description = "A simple test of the uploader",
        ...     color_field = "Temperature",
        ...     color_map = "hot",
        ...     color_log = True,
        ...     bounds = bounds)
        ...
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


