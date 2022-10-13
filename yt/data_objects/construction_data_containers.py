import fileinput
import io
import os
import warnings
import zipfile
from functools import partial, wraps
from re import finditer
from tempfile import NamedTemporaryFile, TemporaryFile
from typing import Tuple

import numpy as np
from more_itertools import always_iterable
from tqdm import tqdm

from yt._maintenance.deprecation import issue_deprecation_warning
from yt.config import ytcfg
from yt.data_objects.field_data import YTFieldData
from yt.data_objects.selection_objects.data_selection_objects import (
    YTSelectionContainer1D,
    YTSelectionContainer2D,
    YTSelectionContainer3D,
)
from yt.fields.field_exceptions import NeedsGridType, NeedsOriginalGrid
from yt.frontends.sph.data_structures import ParticleDataset
from yt.funcs import (
    get_memory_usage,
    is_sequence,
    iter_fields,
    mylog,
    only_on_root,
    validate_moment,
)
from yt.geometry import particle_deposit as particle_deposit
from yt.geometry.coordinates.cartesian_coordinates import all_data
from yt.loaders import load_uniform_grid
from yt.units._numpy_wrapper_functions import uconcatenate
from yt.units.unit_object import Unit  # type: ignore
from yt.units.yt_array import YTArray
from yt.utilities.exceptions import (
    YTNoAPIKey,
    YTNotInsideNotebook,
    YTParticleDepositionNotImplemented,
    YTTooManyVertices,
)
from yt.utilities.grid_data_format.writer import write_to_gdf
from yt.utilities.lib.cyoctree import CyOctree
from yt.utilities.lib.interpolators import ghost_zone_interpolate
from yt.utilities.lib.marching_cubes import march_cubes_grid, march_cubes_grid_flux
from yt.utilities.lib.misc_utilities import fill_region, fill_region_float
from yt.utilities.lib.pixelization_routines import (
    interpolate_sph_grid_gather,
    interpolate_sph_positions_gather,
    normalization_1d_utility,
    normalization_3d_utility,
    pixelize_sph_kernel_arbitrary_grid,
)
from yt.utilities.lib.quad_tree import QuadTree
from yt.utilities.math_utils import compute_stddev_image
from yt.utilities.minimal_representation import MinimalProjectionData
from yt.utilities.parallel_tools.parallel_analysis_interface import (
    communication_system,
    parallel_objects,
    parallel_root_only,
)
from yt.visualization.color_maps import get_colormap_lut


class YTStreamline(YTSelectionContainer1D):
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
    >>> streamlines = Streamlines(ds, [0.5] * 3)
    >>> streamlines.integrate_through_volume()
    >>> stream = streamlines.path(0)
    >>> fig, ax = plt.subplots()
    >>> ax.set_yscale("log")
    >>> ax.plot(stream["t"], stream[("gas", "density")], "-x")
    """

    _type_name = "streamline"
    _con_args = ("positions",)
    sort_by = "t"

    def __init__(self, positions, length=1.0, fields=None, ds=None, **kwargs):
        YTSelectionContainer1D.__init__(self, ds, fields, **kwargs)
        self.positions = positions
        self.dts = np.empty_like(positions[:, 0])
        self.dts[:-1] = np.sqrt(
            np.sum((self.positions[1:] - self.positions[:-1]) ** 2, axis=1)
        )
        self.dts[-1] = self.dts[-1]
        self.length = length
        self.dts /= length
        self.ts = np.add.accumulate(self.dts)
        self._set_center(self.positions[0])
        self.set_field_parameter("center", self.positions[0])
        self._dts, self._ts = {}, {}
        # self._refresh_data()

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
        if field == "dts":
            return self._dts[grid.id]
        if field == "t":
            return self._ts[grid.id]
        return grid[field].flat[mask]

    def _get_cut_mask(self, grid):
        points_in_grid = np.all(self.positions > grid.LeftEdge, axis=1) & np.all(
            self.positions <= grid.RightEdge, axis=1
        )
        pids = np.where(points_in_grid)[0]
        mask = np.zeros(points_in_grid.sum(), dtype="int64")
        dts = np.zeros(points_in_grid.sum(), dtype="float64")
        ts = np.zeros(points_in_grid.sum(), dtype="float64")
        for mi, (i, pos) in enumerate(zip(pids, self.positions[points_in_grid])):
            if not points_in_grid[i]:
                continue
            ci = ((pos - grid.LeftEdge) / grid.dds).astype("int")
            if grid.child_mask[ci[0], ci[1], ci[2]] == 0:
                continue
            for j in range(3):
                ci[j] = min(ci[j], grid.ActiveDimensions[j] - 1)
            mask[mi] = np.ravel_multi_index(ci, grid.ActiveDimensions)
            dts[mi] = self.dts[i]
            ts[mi] = self.ts[i]
        self._dts[grid.id] = dts
        self._ts[grid.id] = ts
        return mask


class YTProj(YTSelectionContainer2D):
    _key_fields = YTSelectionContainer2D._key_fields + ["weight_field"]
    _con_args = ("axis", "field", "weight_field")
    _container_fields = ("px", "py", "pdx", "pdy", "weight_field")

    def __init__(
        self,
        field,
        axis,
        weight_field=None,
        center=None,
        ds=None,
        data_source=None,
        method="integrate",
        field_parameters=None,
        max_level=None,
        *,
        moment=1,
        style=None,
    ):
        super().__init__(axis, ds, field_parameters)
        if style is not None:
            issue_deprecation_warning(
                "The 'style' keyword argument is a deprecated alias for 'method'. "
                "Please use method directly.",
                since="3.2",
                removal="4.2",
                stacklevel=4,
            )
            method = style
        if method == "mip":
            issue_deprecation_warning(
                "The 'mip' method value is a deprecated alias for 'max'. "
                "Please use max directly.",
                since="4.1",
                stacklevel=4,
            )
            method = "max"
        if method == "sum":
            self.method = "integrate"
            self._sum_only = True
        else:
            self.method = method
            self._sum_only = False
        if self.method in ["max", "mip"]:
            self.func = np.max
        elif self.method == "min":
            self.func = np.min
        elif self.method == "integrate":
            self.func = np.sum  # for the future
        else:
            raise NotImplementedError(self.method)
        validate_moment(moment, weight_field)
        self.moment = moment
        self._set_center(center)
        self._projected_units = {}
        if data_source is None:
            data_source = self.ds.all_data()
        if max_level is not None:
            data_source.max_level = max_level
        for k, v in data_source.field_parameters.items():
            if k not in self.field_parameters or self._is_default_field_parameter(k):
                self.set_field_parameter(k, v)
        self.data_source = data_source
        if weight_field is None:
            self.weight_field = weight_field
        else:
            self.weight_field = self._determine_fields(weight_field)[0]

        for f in self._determine_fields(field):
            nodal_flag = self.ds._get_field_info(f).nodal_flag
            if any(nodal_flag):
                raise RuntimeError(
                    "Nodal fields are currently not supported for projections."
                )

    @property
    def blocks(self):
        return self.data_source.blocks

    @property
    def field(self):
        return [k for k in self.field_data.keys() if k not in self._container_fields]

    def get_data(self, fields=None):
        fields = self._determine_fields(fields)
        sfields = []
        if self.moment == 2:

            def _sq_field(field, data, fname: Tuple[str, str]):
                return data[fname] ** 2

            for fname in fields:
                fd = self.ds._get_field_info(*fname)
                self.ds.add_field(
                    (fname[0], f"tmp_{fname[1]}_squared"),
                    partial(_sq_field, fname=fname),
                    sampling_type=fd.sampling_type,
                    units=f"({fd.units})*({fd.units})",
                )
                sfields.append((fname[0], f"tmp_{fname[1]}_squared"))
        nfields = len(fields)
        nsfields = len(sfields)
        # We need a new tree for every single set of fields we add
        if nfields == 0:
            return
        if isinstance(self.ds, ParticleDataset):
            return
        tree = self._get_tree(nfields + nsfields)
        # This only needs to be done if we are in parallel; otherwise, we can
        # safely build the mesh as we go.
        if communication_system.communicators[-1].size > 1:
            for chunk in self.data_source.chunks([], "io", local_only=False):
                self._initialize_chunk(chunk, tree)
        _units_initialized = False
        with self.data_source._field_parameter_state(self.field_parameters):
            for chunk in parallel_objects(
                self.data_source.chunks([], "io", local_only=True)
            ):
                if not _units_initialized:
                    self._initialize_projected_units(fields, chunk)
                    _units_initialized = True
                self._handle_chunk(chunk, fields + sfields, tree)
        # if there's less than nprocs chunks, units won't be initialized
        # on all processors, so sync with _projected_units on rank 0
        projected_units = self.comm.mpi_bcast(self._projected_units)
        self._projected_units = projected_units
        # Note that this will briefly double RAM usage
        if self.method in ["max", "mip"]:
            merge_style = -1
            op = "max"
        elif self.method == "min":
            merge_style = -2
            op = "min"
        elif self.method == "integrate":
            merge_style = 1
            op = "sum"
        else:
            raise NotImplementedError
        # TODO: Add the combine operation
        xax = self.ds.coordinates.x_axis[self.axis]
        yax = self.ds.coordinates.y_axis[self.axis]

        ix, iy, ires, nvals, nwvals = tree.get_all(False, merge_style)
        px, pdx = self.ds.index._icoords_to_fcoords(
            ix[:, None], ires // self.ds.ires_factor, axes=(xax,)
        )
        py, pdy = self.ds.index._icoords_to_fcoords(
            iy[:, None], ires // self.ds.ires_factor, axes=(yax,)
        )
        px = px.ravel()
        py = py.ravel()
        pdx = pdx.ravel()
        pdy = pdy.ravel()
        np.multiply(pdx, 0.5, pdx)
        np.multiply(pdy, 0.5, pdy)

        nvals = self.comm.mpi_allreduce(nvals, op=op)
        nwvals = self.comm.mpi_allreduce(nwvals, op=op)
        if self.weight_field is not None:
            # If there are 0s remaining in the weight vals
            # this will not throw an error, but silently
            # return nans for vals where dividing by 0
            # Leave as NaNs to be auto-masked by Matplotlib
            with np.errstate(invalid="ignore"):
                np.divide(nvals, nwvals[:, None], nvals)
        # We now convert to half-widths and center-points
        data = {}
        code_length = self.ds.domain_width.units
        data["px"] = self.ds.arr(px, code_length)
        data["py"] = self.ds.arr(py, code_length)
        data["weight_field"] = nwvals
        data["pdx"] = self.ds.arr(pdx, code_length)
        data["pdy"] = self.ds.arr(pdy, code_length)
        data["fields"] = nvals
        # Now we run the finalizer, which is ignored if we don't need it
        field_data = np.hsplit(data.pop("fields"), nfields + nsfields)
        for fi, field in enumerate(fields):
            mylog.debug("Setting field %s", field)
            input_units = self._projected_units[field]
            fvals = field_data[fi].ravel()
            if self.moment == 2:
                fvals = compute_stddev_image(field_data[fi + nfields].ravel(), fvals)
            self[field] = self.ds.arr(fvals, input_units)
        for i in list(data.keys()):
            self[i] = data.pop(i)
        mylog.info("Projection completed")
        if self.moment == 2:
            for field in sfields:
                self.ds.field_info.pop(field)
        self.tree = tree

    def to_pw(self, fields=None, center="c", width=None, origin="center-window"):
        r"""Create a :class:`~yt.visualization.plot_window.PWViewerMPL` from this
        object.

        This is a bare-bones mechanism of creating a plot window from this
        object, which can then be moved around, zoomed, and on and on.  All
        behavior of the plot window is relegated to that routine.
        """
        pw = self._get_pw(fields, center, width, origin, "Projection")
        return pw

    def plot(self, fields=None):
        if hasattr(self.data_source, "left_edge") and hasattr(
            self.data_source, "right_edge"
        ):
            left_edge = self.data_source.left_edge
            right_edge = self.data_source.right_edge
            center = (left_edge + right_edge) / 2.0
            width = right_edge - left_edge
            xax = self.ds.coordinates.x_axis[self.axis]
            yax = self.ds.coordinates.y_axis[self.axis]
            lx, rx = left_edge[xax], right_edge[xax]
            ly, ry = left_edge[yax], right_edge[yax]
            width = (rx - lx), (ry - ly)
        else:
            width = self.ds.domain_width
            center = self.ds.domain_center
        pw = self._get_pw(fields, center, width, "native", "Projection")
        try:
            pw.show()
        except YTNotInsideNotebook:
            pass
        return pw

    def _initialize_projected_units(self, fields, chunk):
        for field in self.data_source._determine_fields(fields):
            if field in self._projected_units:
                continue
            finfo = self.ds._get_field_info(*field)
            if finfo.units is None:
                # First time calling a units="auto" field, infer units and cache
                # for future field accesses.
                finfo.units = str(chunk[field].units)
            field_unit = Unit(finfo.output_units, registry=self.ds.unit_registry)
            if self.method in ("min", "max") or self._sum_only:
                path_length_unit = Unit(registry=self.ds.unit_registry)
            else:
                ax_name = self.ds.coordinates.axis_name[self.axis]
                path_element_name = ("index", f"path_element_{ax_name}")
                path_length_unit = self.ds.field_info[path_element_name].units
                path_length_unit = Unit(
                    path_length_unit, registry=self.ds.unit_registry
                )
                # Only convert to appropriate unit system for path
                # elements that aren't angles
                if not path_length_unit.is_dimensionless:
                    path_length_unit = path_length_unit.get_base_equivalent(
                        unit_system=self.ds.unit_system
                    )
            if self.weight_field is None:
                self._projected_units[field] = field_unit * path_length_unit
            else:
                self._projected_units[field] = field_unit


class YTParticleProj(YTProj):
    """
    A projection operation optimized for SPH particle data.
    """

    _type_name = "particle_proj"

    def __init__(
        self,
        field,
        axis,
        weight_field=None,
        center=None,
        ds=None,
        data_source=None,
        method="integrate",
        field_parameters=None,
        max_level=None,
        *,
        moment=1,
        style=None,
    ):
        super().__init__(
            field,
            axis,
            weight_field,
            center,
            ds,
            data_source,
            method,
            field_parameters,
            max_level,
            moment=moment,
            style=style,
        )

    def _handle_chunk(self, chunk, fields, tree):
        raise NotImplementedError("Particle projections have not yet been implemented")


class YTQuadTreeProj(YTProj):
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
        "mip" : maximum intensity projection (deprecated)
        "max" : maximum intensity projection
        "min" : minimum intensity projection
        "sum" : same as "integrate", except that we don't multiply by the path length
        WARNING: The "sum" option should only be used for uniform resolution grid
        datasets, as other datasets may result in unphysical images.
    style : string, optional
        The same as the method keyword.  Deprecated as of version 3.0.2.
        Please use method keyword instead.
    field_parameters : dict of items
        Values to be passed as field parameters that can be
        accessed by generated fields.
    moment : integer, optional
        for a weighted projection, moment = 1 (the default) corresponds to a
        weighted average. moment = 2 corresponds to a weighted standard
        deviation.

    Examples
    --------

    >>> ds = load("RedshiftOutput0005")
    >>> prj = ds.proj(("gas", "density"), 0)
    >>> print(proj[("gas", "density")])
    """

    _type_name = "quad_proj"

    def __init__(
        self,
        field,
        axis,
        weight_field=None,
        center=None,
        ds=None,
        data_source=None,
        method="integrate",
        field_parameters=None,
        max_level=None,
        *,
        moment=1,
        style=None,
    ):
        super().__init__(
            field,
            axis,
            weight_field,
            center,
            ds,
            data_source,
            method,
            field_parameters,
            max_level,
            moment=moment,
            style=style,
        )

        if not self.deserialize(field):
            self.get_data(field)
            self.serialize()

    @property
    def _mrep(self):
        return MinimalProjectionData(self)

    def deserialize(self, fields):
        if not ytcfg.get("yt", "serialize"):
            return False
        for field in fields:
            self[field] = None
        deserialized_successfully = False
        store_file = self.ds.parameter_filename + ".yt"
        if os.path.isfile(store_file):
            deserialized_successfully = self._mrep.restore(store_file, self.ds)

            if deserialized_successfully:
                mylog.info("Using previous projection data from %s", store_file)
                for field, field_data in self._mrep.field_data.items():
                    self[field] = field_data
        if not deserialized_successfully:
            for field in fields:
                del self[field]
        return deserialized_successfully

    def serialize(self):
        if not ytcfg.get("yt", "serialize"):
            return
        self._mrep.store(self.ds.parameter_filename + ".yt")

    def _get_tree(self, nvals):
        xax = self.ds.coordinates.x_axis[self.axis]
        yax = self.ds.coordinates.y_axis[self.axis]
        xd = self.ds.domain_dimensions[xax]
        yd = self.ds.domain_dimensions[yax]
        bounds = (
            self.ds.domain_left_edge[xax],
            self.ds.domain_right_edge[xax],
            self.ds.domain_left_edge[yax],
            self.ds.domain_right_edge[yax],
        )
        return QuadTree(
            np.array([xd, yd], dtype="int64"), nvals, bounds, method=self.method
        )

    def _initialize_chunk(self, chunk, tree):
        icoords = chunk.icoords
        xax = self.ds.coordinates.x_axis[self.axis]
        yax = self.ds.coordinates.y_axis[self.axis]
        i1 = icoords[:, xax]
        i2 = icoords[:, yax]
        ilevel = chunk.ires * self.ds.ires_factor
        tree.initialize_chunk(i1, i2, ilevel)

    def _handle_chunk(self, chunk, fields, tree):
        mylog.debug(
            "Adding chunk (%s) to tree (%0.3e GB RAM)",
            chunk.ires.size,
            get_memory_usage() / 1024.0,
        )
        if self.method in ("min", "max") or self._sum_only:
            dl = self.ds.quan(1.0, "")
        else:
            ax_name = self.ds.coordinates.axis_name[self.axis]
            dl = chunk["index", f"path_element_{ax_name}"]
            # This is done for cases where our path element does not have a CGS
            # equivalent.  Once "preferred units" have been implemented, this
            # will not be necessary at all, as the final conversion will occur
            # at the display layer.
            if not dl.units.is_dimensionless:
                dl.convert_to_units(self.ds.unit_system["length"])
        v = np.empty((chunk.ires.size, len(fields)), dtype="float64")
        for i, field in enumerate(fields):
            d = chunk[field] * dl
            v[:, i] = d
        if self.weight_field is not None:
            w = chunk[self.weight_field]
            np.multiply(v, w[:, None], v)
            np.multiply(w, dl, w)
        else:
            w = np.ones(chunk.ires.size, dtype="float64")
        icoords = chunk.icoords
        xax = self.ds.coordinates.x_axis[self.axis]
        yax = self.ds.coordinates.y_axis[self.axis]
        i1 = icoords[:, xax]
        i2 = icoords[:, yax]
        ilevel = chunk.ires * self.ds.ires_factor
        tree.add_chunk_to_tree(i1, i2, ilevel, v, w)


class YTCoveringGrid(YTSelectionContainer3D):
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
    num_ghost_zones : integer, optional
        The number of padding ghost zones used when accessing fields.
    data_source :
        An existing data object to intersect with the covering grid. Grid points
        outside the data_source will exist as empty values.

    Examples
    --------
    >>> cube = ds.covering_grid(2, left_edge=[0.0, 0.0, 0.0], dims=[128, 128, 128])
    """

    _spatial = True
    _type_name = "covering_grid"
    _con_args = ("level", "left_edge", "ActiveDimensions")
    _container_fields = (
        ("index", "dx"),
        ("index", "dy"),
        ("index", "dz"),
        ("index", "x"),
        ("index", "y"),
        ("index", "z"),
    )
    _base_grid = None

    def __init__(
        self,
        level,
        left_edge,
        dims,
        fields=None,
        ds=None,
        num_ghost_zones=0,
        use_pbar=True,
        field_parameters=None,
        *,
        data_source=None,
    ):
        if field_parameters is None:
            center = None
        else:
            center = field_parameters.get("center", None)
        super().__init__(center, ds, field_parameters, data_source=data_source)

        self.level = level
        self.left_edge = self._sanitize_edge(left_edge)
        self.ActiveDimensions = self._sanitize_dims(dims)

        rdx = self.ds.domain_dimensions * self.ds.relative_refinement(0, level)

        # normalize dims as a non-zero dim array
        dims = np.array(list(always_iterable(dims)))
        rdx[np.where(dims - 2 * num_ghost_zones <= 1)] = 1  # issue 602
        self.base_dds = self.ds.domain_width / self.ds.domain_dimensions
        self.dds = self.ds.domain_width / rdx.astype("float64")
        self.right_edge = self.left_edge + self.ActiveDimensions * self.dds
        self._num_ghost_zones = num_ghost_zones
        self._use_pbar = use_pbar
        self.global_startindex = (
            np.rint((self.left_edge - self.ds.domain_left_edge) / self.dds).astype(
                "int64"
            )
            + self.ds.domain_offset
        )
        self._setup_data_source()
        self.get_data(fields)

    def get_global_startindex(self):
        r"""Get the global start index of the covering grid."""
        return self.global_startindex

    def to_xarray(self, fields=None):
        r"""Export this fixed-resolution object to an xarray Dataset

        This function will take a regularized grid and optionally a list of
        fields and return an xarray Dataset object.  If xarray is not
        importable, this will raise ImportError.

        Parameters
        ----------
        fields : list of strings or tuple field names, default None
            If this is supplied, it is the list of fields to be exported into
            the data frame.  If not supplied, whatever fields presently exist
            will be used.

        Returns
        -------
        arr : Dataset
            The data contained in the object.

        Examples
        --------

        >>> dd = ds.r[::256j, ::256j, ::256j]
        >>> xf1 = dd.to_xarray([("gas", "density"), ("gas", "temperature")])
        >>> dd[("gas", "velocity_magnitude")]
        >>> xf2 = dd.to_xarray()
        """
        import xarray as xr

        data = {}
        coords = {}
        for f in fields or self.field_data.keys():
            data[f] = {
                "dims": (
                    "x",
                    "y",
                    "z",
                ),
                "data": self[f],
                "attrs": {"units": str(self[f].uq)},
            }
        # We have our data, so now we generate both our coordinates and our metadata.
        LE = self.LeftEdge + self.dds / 2.0
        RE = self.RightEdge - self.dds / 2.0
        N = self.ActiveDimensions
        u = str(LE.uq)
        for i, ax in enumerate("xyz"):
            coords[ax] = {
                "dims": (ax,),
                "data": np.mgrid[LE[i] : RE[i] : N[i] * 1j],
                "attrs": {"units": u},
            }
        return xr.Dataset.from_dict({"data_vars": data, "coords": coords})

    @property
    def icoords(self):
        ic = np.indices(self.ActiveDimensions).astype("int64")
        return np.column_stack(
            [i.ravel() + gi for i, gi in zip(ic, self.get_global_startindex())]
        )

    @property
    def fwidth(self):
        fw = np.ones((self.ActiveDimensions.prod(), 3), dtype="float64")
        fw *= self.dds
        return fw

    @property
    def fcoords(self):
        LE = self.LeftEdge + self.dds / 2.0
        RE = self.RightEdge - self.dds / 2.0
        N = self.ActiveDimensions
        fc = np.mgrid[
            LE[0] : RE[0] : N[0] * 1j,
            LE[1] : RE[1] : N[1] * 1j,
            LE[2] : RE[2] : N[2] * 1j,
        ]
        return np.column_stack([f.ravel() for f in fc])

    @property
    def ires(self):
        tr = np.ones(self.ActiveDimensions.prod(), dtype="int64")
        tr *= self.level
        return tr

    def set_field_parameter(self, name, val):
        super().set_field_parameter(name, val)
        if self._data_source is not None:
            self._data_source.set_field_parameter(name, val)

    def _sanitize_dims(self, dims):
        if not is_sequence(dims):
            dims = [dims] * len(self.ds.domain_left_edge)
        if len(dims) != len(self.ds.domain_left_edge):
            raise RuntimeError(
                "Length of dims must match the dimensionality of the dataset"
            )
        return np.array(dims, dtype="int32")

    def _sanitize_edge(self, edge):
        if not is_sequence(edge):
            edge = [edge] * len(self.ds.domain_left_edge)
        if len(edge) != len(self.ds.domain_left_edge):
            raise RuntimeError(
                "Length of edges must match the dimensionality of the dataset"
            )
        if hasattr(edge, "units"):
            if edge.units.registry is self.ds.unit_registry:
                return edge
            edge_units = edge.units.copy()
            edge_units.registry = self.ds.unit_registry
        else:
            edge_units = "code_length"
        return self.ds.arr(edge, edge_units, dtype="float64")

    def _reshape_vals(self, arr):
        if len(arr.shape) == 3:
            return arr
        return arr.reshape(self.ActiveDimensions, order="C")

    @property
    def shape(self):
        return tuple(self.ActiveDimensions.tolist())

    def _setup_data_source(self):

        reg = self.ds.region(self.center, self.left_edge, self.right_edge)
        if self._data_source is None:
            # note: https://github.com/yt-project/yt/pull/4063 implemented
            # a data_source kwarg for YTCoveringGrid, but not YTArbitraryGrid
            # so as of 4063, this will always be True for YTArbitraryGrid
            # instances.
            self._data_source = reg
        else:
            self._data_source = self.ds.intersection([self._data_source, reg])

        self._data_source.min_level = 0
        self._data_source.max_level = self.level
        # This triggers "special" behavior in the RegionSelector to ensure we
        # select *cells* whose bounding boxes overlap with our region, not just
        # their cell centers.
        self._data_source.loose_selection = True

    def get_data(self, fields=None):
        if fields is None:
            return
        fields = self._determine_fields(fields)
        fields_to_get = [f for f in fields if f not in self.field_data]
        fields_to_get = self._identify_dependencies(fields_to_get)
        if len(fields_to_get) == 0:
            return
        try:
            fill, gen, part, alias = self._split_fields(fields_to_get)
        except NeedsGridType as e:
            if self._num_ghost_zones == 0:
                raise RuntimeError(
                    "Attempting to access a field that needs ghost zones, but "
                    "num_ghost_zones = %s. You should create the covering grid "
                    "with nonzero num_ghost_zones." % self._num_ghost_zones
                ) from e
            else:
                raise

        # checking if we have a sph particles
        if len(part) == 0:
            is_sph_field = False
        else:
            is_sph_field = self.ds.field_info[part[0]].is_sph_field

        if len(part) > 0 and len(alias) == 0:
            if is_sph_field:
                self._fill_sph_particles(fields)
                for field in fields:
                    if field in gen:
                        gen.remove(field)
            else:
                self._fill_particles(part)

        if len(fill) > 0:
            self._fill_fields(fill)
        for a, f in sorted(alias.items()):
            if f.sampling_type == "particle" and not is_sph_field:
                self[a] = self._data_source[f]
            else:
                self[a] = f(self)
            self.field_data[a].convert_to_units(f.output_units)

        if len(gen) > 0:
            part_gen = []
            cell_gen = []
            for field in gen:
                finfo = self.ds.field_info[field]
                if finfo.sampling_type == "particle":
                    part_gen.append(field)
                else:
                    cell_gen.append(field)
            self._generate_fields(cell_gen)
            for p in part_gen:
                self[p] = self._data_source[p]

    def _split_fields(self, fields_to_get):
        fill, gen = self.index._split_fields(fields_to_get)
        particles = []
        alias = {}
        for field in gen:
            finfo = self.ds._get_field_info(*field)
            if finfo.is_alias:
                alias[field] = finfo
                continue
            try:
                finfo.check_available(self)
            except NeedsOriginalGrid:
                fill.append(field)
        for field in fill:
            finfo = self.ds._get_field_info(*field)
            if finfo.sampling_type == "particle":
                particles.append(field)
        gen = [f for f in gen if f not in fill and f not in alias]
        fill = [f for f in fill if f not in particles]
        return fill, gen, particles, alias

    def _fill_particles(self, part):
        for p in part:
            self[p] = self._data_source[p]

    def _fill_sph_particles(self, fields):
        # checks that we have the field and gets information
        fields = [f for f in fields if f not in self.field_data]
        if len(fields) == 0:
            return

        smoothing_style = getattr(self.ds, "sph_smoothing_style", "scatter")
        normalize = getattr(self.ds, "use_sph_normalization", True)

        bounds, size = self._get_grid_bounds_size()

        period = self.ds.coordinates.period.copy()
        if hasattr(period, "in_units"):
            period = period.in_units("code_length").d
        # TODO maybe there is a better way of handling this
        is_periodic = int(any(self.ds.periodicity))

        if smoothing_style == "scatter":
            for field in fields:
                fi = self.ds._get_field_info(field)
                ptype = fi.name[0]
                if ptype not in self.ds._sph_ptypes:
                    raise KeyError(f"{ptype} is not a SPH particle type!")
                buff = np.zeros(size, dtype="float64")
                if normalize:
                    buff_den = np.zeros(size, dtype="float64")

                pbar = tqdm(desc=f"Interpolating SPH field {field}")
                for chunk in self._data_source.chunks([field], "io"):
                    px = chunk[(ptype, "particle_position_x")].in_base("code").d
                    py = chunk[(ptype, "particle_position_y")].in_base("code").d
                    pz = chunk[(ptype, "particle_position_z")].in_base("code").d
                    hsml = chunk[(ptype, "smoothing_length")].in_base("code").d
                    mass = chunk[(ptype, "particle_mass")].in_base("code").d
                    dens = chunk[(ptype, "density")].in_base("code").d
                    field_quantity = chunk[field].d

                    pixelize_sph_kernel_arbitrary_grid(
                        buff,
                        px,
                        py,
                        pz,
                        hsml,
                        mass,
                        dens,
                        field_quantity,
                        bounds,
                        pbar=pbar,
                        check_period=is_periodic,
                        period=period,
                    )
                    if normalize:
                        pixelize_sph_kernel_arbitrary_grid(
                            buff_den,
                            px,
                            py,
                            pz,
                            hsml,
                            mass,
                            dens,
                            np.ones(dens.shape[0]),
                            bounds,
                            pbar=pbar,
                            check_period=is_periodic,
                            period=period,
                        )

                if normalize:
                    normalization_3d_utility(buff, buff_den)

                self[field] = self.ds.arr(buff, fi.units)
                pbar.close()

        if smoothing_style == "gather":
            num_neighbors = getattr(self.ds, "num_neighbors", 32)
            for field in fields:
                buff = np.zeros(size, dtype="float64")

                fields_to_get = [
                    "particle_position",
                    "density",
                    "particle_mass",
                    "smoothing_length",
                    field[1],
                ]
                all_fields = all_data(self.ds, field[0], fields_to_get, kdtree=True)

                fi = self.ds._get_field_info(field)
                interpolate_sph_grid_gather(
                    buff,
                    all_fields["particle_position"],
                    bounds,
                    all_fields["smoothing_length"],
                    all_fields["particle_mass"],
                    all_fields["density"],
                    all_fields[field[1]].in_units(fi.units),
                    self.ds.index.kdtree,
                    use_normalization=normalize,
                    num_neigh=num_neighbors,
                )

                self[field] = self.ds.arr(buff, fi.units)

    def _fill_fields(self, fields):
        fields = [f for f in fields if f not in self.field_data]
        if len(fields) == 0:
            return
        output_fields = [
            np.zeros(self.ActiveDimensions, dtype="float64") for field in fields
        ]
        domain_dims = self.ds.domain_dimensions.astype(
            "int64"
        ) * self.ds.relative_refinement(0, self.level)
        refine_by = self.ds.refine_by
        if not is_sequence(self.ds.refine_by):
            refine_by = [refine_by, refine_by, refine_by]
        refine_by = np.array(refine_by, dtype="i8")
        for chunk in parallel_objects(self._data_source.chunks(fields, "io")):
            input_fields = [chunk[field] for field in fields]
            # NOTE: This usage of "refine_by" is actually *okay*, because it's
            # being used with respect to iref, which is *already* scaled!
            fill_region(
                input_fields,
                output_fields,
                self.level,
                self.global_startindex,
                chunk.icoords,
                chunk.ires,
                domain_dims,
                refine_by,
            )
        if self.comm.size > 1:
            for i in range(len(fields)):
                output_fields[i] = self.comm.mpi_allreduce(output_fields[i], op="sum")
        for name, v in zip(fields, output_fields):
            fi = self.ds._get_field_info(*name)
            self[name] = self.ds.arr(v, fi.units)

    def _generate_container_field(self, field):
        rv = self.ds.arr(np.ones(self.ActiveDimensions, dtype="float64"), "")
        axis_name = self.ds.coordinates.axis_name
        if field == ("index", f"d{axis_name[0]}"):
            np.multiply(rv, self.dds[0], rv)
        elif field == ("index", f"d{axis_name[1]}"):
            np.multiply(rv, self.dds[1], rv)
        elif field == ("index", f"d{axis_name[2]}"):
            np.multiply(rv, self.dds[2], rv)
        elif field == ("index", axis_name[0]):
            x = np.mgrid[
                self.left_edge[0]
                + 0.5 * self.dds[0] : self.right_edge[0]
                - 0.5 * self.dds[0] : self.ActiveDimensions[0] * 1j
            ]
            np.multiply(rv, x[:, None, None], rv)
        elif field == ("index", axis_name[1]):
            y = np.mgrid[
                self.left_edge[1]
                + 0.5 * self.dds[1] : self.right_edge[1]
                - 0.5 * self.dds[1] : self.ActiveDimensions[1] * 1j
            ]
            np.multiply(rv, y[None, :, None], rv)
        elif field == ("index", axis_name[2]):
            z = np.mgrid[
                self.left_edge[2]
                + 0.5 * self.dds[2] : self.right_edge[2]
                - 0.5 * self.dds[2] : self.ActiveDimensions[2] * 1j
            ]
            np.multiply(rv, z[None, None, :], rv)
        else:
            raise KeyError(field)
        return rv

    @property
    def LeftEdge(self):
        return self.left_edge

    @property
    def RightEdge(self):
        return self.right_edge

    def deposit(self, positions, fields=None, method=None, kernel_name="cubic"):
        cls = getattr(particle_deposit, f"deposit_{method}", None)
        if cls is None:
            raise YTParticleDepositionNotImplemented(method)
        # We allocate number of zones, not number of octs. Everything
        # inside this is Fortran ordered because of the ordering in the
        # octree deposit routines, so we reverse it here to match the
        # convention there
        nvals = tuple(self.ActiveDimensions[::-1])
        # append a dummy dimension because we are only depositing onto
        # one grid
        op = cls(nvals + (1,), kernel_name)
        op.initialize()
        op.process_grid(self, positions, fields)
        # Fortran-ordered, so transpose.
        vals = op.finalize().transpose()
        # squeeze dummy dimension we appended above
        return np.squeeze(vals, axis=0)

    def write_to_gdf(self, gdf_path, fields, nprocs=1, field_units=None, **kwargs):
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
        >>> cube.write_to_gdf(
        ...     "clumps.h5",
        ...     [("gas", "density"), ("gas", "temperature")],
        ...     nprocs=16,
        ...     overwrite=True,
        ... )
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
        bbox = np.array([[l, r] for l, r in zip(le, re)])
        ds = load_uniform_grid(
            data,
            self.ActiveDimensions,
            bbox=bbox,
            length_unit=self.ds.length_unit,
            time_unit=self.ds.time_unit,
            mass_unit=self.ds.mass_unit,
            nprocs=nprocs,
            sim_time=self.ds.current_time.v,
        )
        write_to_gdf(ds, gdf_path, **kwargs)

    def _get_grid_bounds_size(self):
        dd = self.ds.domain_width / 2**self.level
        bounds = np.zeros(6, dtype="float64")

        bounds[0] = self.left_edge[0].in_base("code")
        bounds[1] = bounds[0] + dd[0].d * self.ActiveDimensions[0]
        bounds[2] = self.left_edge[1].in_base("code")
        bounds[3] = bounds[2] + dd[1].d * self.ActiveDimensions[1]
        bounds[4] = self.left_edge[2].in_base("code")
        bounds[5] = bounds[4] + dd[2].d * self.ActiveDimensions[2]
        size = np.ones(3, dtype="int64") * 2**self.level

        return bounds, size

    def to_fits_data(self, fields, length_unit=None):
        r"""Export a set of gridded fields to a FITS file.

        This will export a set of FITS images of either the fields specified
        or all the fields already in the object.

        Parameters
        ----------
        fields : list of strings
            These fields will be pixelized and output. If "None", the keys of the
            FRB will be used.
        length_unit : string, optional
            the length units that the coordinates are written in. The default
            is to use the default length unit of the dataset.
        """
        from yt.visualization.fits_image import FITSImageData

        if length_unit is None:
            length_unit = self.ds.length_unit
        fields = list(iter_fields(fields))
        fid = FITSImageData(self, fields, length_unit=length_unit)
        return fid


class YTArbitraryGrid(YTCoveringGrid):
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
    >>> obj = ds.arbitrary_grid(
    ...     [0.0, 0.0, 0.0], [0.99, 0.99, 0.99], dims=[128, 128, 128]
    ... )
    """

    _spatial = True
    _type_name = "arbitrary_grid"
    _con_args = ("left_edge", "right_edge", "ActiveDimensions")
    _container_fields = (
        ("index", "dx"),
        ("index", "dy"),
        ("index", "dz"),
        ("index", "x"),
        ("index", "y"),
        ("index", "z"),
    )

    def __init__(self, left_edge, right_edge, dims, ds=None, field_parameters=None):
        if field_parameters is None:
            center = None
        else:
            center = field_parameters.get("center", None)
        YTSelectionContainer3D.__init__(self, center, ds, field_parameters)
        self.left_edge = self._sanitize_edge(left_edge)
        self.right_edge = self._sanitize_edge(right_edge)
        self.ActiveDimensions = self._sanitize_dims(dims)
        self.dds = self.base_dds = (
            self.right_edge - self.left_edge
        ) / self.ActiveDimensions
        self.level = 99
        self._setup_data_source()

    def _fill_fields(self, fields):
        fields = [f for f in fields if f not in self.field_data]
        if len(fields) == 0:
            return
        # It may be faster to adapt fill_region_float to fill multiple fields
        # instead of looping here
        for field in fields:
            dest = np.zeros(self.ActiveDimensions, dtype="float64")
            for chunk in self._data_source.chunks(fields, "io"):
                fill_region_float(
                    chunk.fcoords,
                    chunk.fwidth,
                    chunk[field],
                    self.left_edge,
                    self.right_edge,
                    dest,
                    1,
                    self.ds.domain_width,
                    int(any(self.ds.periodicity)),
                )
            fi = self.ds._get_field_info(field)
            self[field] = self.ds.arr(dest, fi.units)

    def _get_grid_bounds_size(self):
        bounds = np.empty(6, dtype="float64")
        bounds[0] = self.left_edge[0].in_base("code")
        bounds[2] = self.left_edge[1].in_base("code")
        bounds[4] = self.left_edge[2].in_base("code")
        bounds[1] = self.right_edge[0].in_base("code")
        bounds[3] = self.right_edge[1].in_base("code")
        bounds[5] = self.right_edge[2].in_base("code")
        size = self.ActiveDimensions

        return bounds, size


class LevelState:
    current_dx = None
    current_dims = None
    current_level = None
    global_startindex = None
    old_global_startindex = None
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


class YTSmoothedCoveringGrid(YTCoveringGrid):
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
    _min_level = None

    @wraps(YTCoveringGrid.__init__)
    def __init__(self, *args, **kwargs):
        ds = kwargs["ds"]
        self._base_dx = (
            ds.domain_right_edge - ds.domain_left_edge
        ) / ds.domain_dimensions.astype("float64")
        self.global_endindex = None
        YTCoveringGrid.__init__(self, *args, **kwargs)
        self._final_start_index = self.global_startindex

    def _setup_data_source(self, level_state=None):
        if level_state is None:
            super()._setup_data_source()
            return
        # We need a buffer region to allow for zones that contribute to the
        # interpolation but are not directly inside our bounds
        level_state.data_source = self.ds.region(
            self.center,
            level_state.left_edge - level_state.current_dx,
            level_state.right_edge + level_state.current_dx,
        )
        level_state.data_source.min_level = level_state.current_level
        level_state.data_source.max_level = level_state.current_level
        self._pdata_source = self.ds.region(
            self.center,
            level_state.left_edge - level_state.current_dx,
            level_state.right_edge + level_state.current_dx,
        )
        self._pdata_source.min_level = level_state.current_level
        self._pdata_source.max_level = level_state.current_level

    def _compute_minimum_level(self):
        # This attempts to determine the minimum level that we should be
        # starting on for this box.  It does this by identifying the minimum
        # level that could contribute to the minimum bounding box at that
        # level; that means that all cells from coarser levels will be replaced.
        if self._min_level is not None:
            return self._min_level
        ils = LevelState()
        min_level = 0
        for l in range(self.level, 0, -1):
            dx = self._base_dx / self.ds.relative_refinement(0, l)
            start_index, end_index, dims = self._minimal_box(dx)
            ils.left_edge = start_index * dx + self.ds.domain_left_edge
            ils.right_edge = ils.left_edge + dx * dims
            ils.current_dx = dx
            ils.current_level = l
            self._setup_data_source(ils)
            # Reset the max_level
            ils.data_source.min_level = 0
            ils.data_source.max_level = l
            ils.data_source.loose_selection = False
            min_level = self.level
            for chunk in ils.data_source.chunks([], "io"):
                # With our odd selection methods, we can sometimes get no-sized ires.
                ir = chunk.ires
                if ir.size == 0:
                    continue
                min_level = min(ir.min(), min_level)
            if min_level >= l:
                break
        self._min_level = min_level
        return min_level

    def _fill_fields(self, fields):
        fields = [f for f in fields if f not in self.field_data]
        if len(fields) == 0:
            return
        ls = self._initialize_level_state(fields)
        min_level = self._compute_minimum_level()
        # NOTE: This usage of "refine_by" is actually *okay*, because it's
        # being used with respect to iref, which is *already* scaled!
        refine_by = self.ds.refine_by
        if not is_sequence(self.ds.refine_by):
            refine_by = [refine_by, refine_by, refine_by]
        refine_by = np.array(refine_by, dtype="i8")

        runtime_errors_count = 0
        for level in range(self.level + 1):
            if level < min_level:
                self._update_level_state(ls)
                continue
            nd = self.ds.dimensionality
            refinement = np.zeros_like(ls.base_dx)
            refinement += self.ds.relative_refinement(0, ls.current_level)
            refinement[nd:] = 1
            domain_dims = self.ds.domain_dimensions * refinement
            domain_dims = domain_dims.astype("int64")
            tot = ls.current_dims.prod()
            for chunk in ls.data_source.chunks(fields, "io"):
                chunk[fields[0]]
                input_fields = [chunk[field] for field in fields]
                tot -= fill_region(
                    input_fields,
                    ls.fields,
                    ls.current_level,
                    ls.global_startindex,
                    chunk.icoords,
                    chunk.ires,
                    domain_dims,
                    refine_by,
                )
            if level == 0 and tot != 0:
                runtime_errors_count += 1
            self._update_level_state(ls)
        if runtime_errors_count:
            warnings.warn(
                "Something went wrong during field computation. "
                "This is likely due to missing ghost-zones support "
                f"in class {type(self.ds)}",
                category=RuntimeWarning,
            )
            mylog.debug("Caught %d runtime errors.", runtime_errors_count)
        for name, v in zip(fields, ls.fields):
            if self.level > 0:
                v = v[1:-1, 1:-1, 1:-1]
            fi = self.ds._get_field_info(*name)
            self[name] = self.ds.arr(v, fi.units)

    def _initialize_level_state(self, fields):
        ls = LevelState()
        ls.domain_width = self.ds.domain_width
        ls.domain_left_edge = self.ds.domain_left_edge
        ls.domain_right_edge = self.ds.domain_right_edge
        ls.base_dx = self._base_dx
        ls.dds = self.dds
        ls.left_edge = self.left_edge
        ls.right_edge = self.right_edge
        for att in (
            "domain_width",
            "domain_left_edge",
            "domain_right_edge",
            "left_edge",
            "right_edge",
            "base_dx",
            "dds",
        ):
            setattr(ls, att, getattr(ls, att).in_units("code_length").d)
        ls.current_dx = ls.base_dx
        ls.current_level = 0
        ls.global_startindex, end_index, idims = self._minimal_box(ls.current_dx)
        ls.current_dims = idims.astype("int32")
        ls.left_edge = ls.global_startindex * ls.current_dx + self.ds.domain_left_edge.d
        ls.right_edge = ls.left_edge + ls.current_dims * ls.current_dx
        ls.fields = [np.zeros(idims, dtype="float64") - 999 for field in fields]
        self._setup_data_source(ls)
        return ls

    def _minimal_box(self, dds):
        LL = self.left_edge.d - self.ds.domain_left_edge.d
        # Nudge in case we're on the edge
        LL += np.finfo(np.float64).eps
        LS = self.right_edge.d - self.ds.domain_left_edge.d
        LS += np.finfo(np.float64).eps
        cell_start = LL / dds  # This is the cell we're inside
        cell_end = LS / dds
        if self.level == 0:
            start_index = np.array(np.floor(cell_start), dtype="int64")
            end_index = np.array(np.ceil(cell_end), dtype="int64")
            dims = np.rint((self.ActiveDimensions * self.dds.d) / dds).astype("int64")
        else:
            # Give us one buffer
            start_index = np.rint(cell_start).astype("int64") - 1
            # How many root cells do we occupy?
            end_index = np.rint(cell_end).astype("int64")
            dims = end_index - start_index + 1
        return start_index, end_index.astype("int64"), dims.astype("int32")

    def _update_level_state(self, level_state):
        ls = level_state
        if ls.current_level >= self.level:
            return
        rf = float(self.ds.relative_refinement(ls.current_level, ls.current_level + 1))
        ls.current_level += 1
        nd = self.ds.dimensionality
        refinement = np.zeros_like(ls.base_dx)
        refinement += self.ds.relative_refinement(0, ls.current_level)
        refinement[nd:] = 1
        ls.current_dx = ls.base_dx / refinement
        ls.old_global_startindex = ls.global_startindex
        ls.global_startindex, end_index, ls.current_dims = self._minimal_box(
            ls.current_dx
        )
        ls.left_edge = ls.global_startindex * ls.current_dx + self.ds.domain_left_edge.d
        ls.right_edge = ls.left_edge + ls.current_dims * ls.current_dx
        input_left = (level_state.old_global_startindex) * rf + 1
        new_fields = []
        for input_field in level_state.fields:
            output_field = np.zeros(ls.current_dims, dtype="float64")
            output_left = level_state.global_startindex + 0.5
            ghost_zone_interpolate(
                rf, input_field, input_left, output_field, output_left
            )
            new_fields.append(output_field)
        level_state.fields = new_fields
        self._setup_data_source(ls)


class YTSurface(YTSelectionContainer3D):
    r"""This surface object identifies isocontours on a cell-by-cell basis,
    with no consideration of global connectedness, and returns the vertices
    of the Triangles in that isocontour.

    This object simply returns the vertices of all the triangles calculated by
    the `marching cubes <https://en.wikipedia.org/wiki/Marching_cubes>`_
    algorithm; for more complex operations, such as identifying connected sets
    of cells above a given threshold, see the extract_connected_sets function.
    This is more useful for calculating, for instance, total isocontour area, or
    visualizing in an external program (such as `MeshLab
    <http://www.meshlab.net>`_.)  The object has the properties .vertices and
    will sample values if a field is requested.  The values are interpolated to
    the center of a given face.

    Parameters
    ----------
    data_source : YTSelectionContainer
        This is the object which will used as a source
    surface_field : string
        Any field that can be obtained in a data object.  This is the field
        which will be isocontoured.
    field_value : float, YTQuantity, or unit tuple
        The value at which the isocontour should be calculated.

    Examples
    --------
    This will create a data object, find a nice value in the center, and
    output the vertices to "triangles.obj" after rescaling them.

    >>> from yt.units import kpc
    >>> sp = ds.sphere("max", (10, "kpc"))
    >>> surf = ds.surface(sp, ("gas", "density"), 5e-27)
    >>> print(surf[("gas", "temperature")])
    >>> print(surf.vertices)
    >>> bounds = [
    ...     (sp.center[i] - 5.0 * kpc, sp.center[i] + 5.0 * kpc) for i in range(3)
    ... ]
    >>> surf.export_ply("my_galaxy.ply", bounds=bounds)
    """
    _type_name = "surface"
    _con_args = ("data_source", "surface_field", "field_value")
    _container_fields = (
        ("index", "dx"),
        ("index", "dy"),
        ("index", "dz"),
        ("index", "x"),
        ("index", "y"),
        ("index", "z"),
    )

    def __init__(self, data_source, surface_field, field_value, ds=None):
        self.data_source = data_source
        self.surface_field = data_source._determine_fields(surface_field)[0]
        finfo = data_source.ds.field_info[self.surface_field]
        try:
            self.field_value = field_value.to(finfo.units)
        except AttributeError:
            if isinstance(field_value, tuple):
                self.field_value = data_source.ds.quan(*field_value)
                self.field_value = self.field_value.to(finfo.units)
            else:
                self.field_value = data_source.ds.quan(field_value, finfo.units)
        self.vertex_samples = YTFieldData()
        center = data_source.get_field_parameter("center")
        super().__init__(center=center, ds=ds)

    def _generate_container_field(self, field):
        self.get_data(field)
        return self[field]

    def get_data(self, fields=None, sample_type="face", no_ghost=False):
        if isinstance(fields, list) and len(fields) > 1:
            for field in fields:
                self.get_data(field)
            return
        elif isinstance(fields, list):
            fields = fields[0]
        # Now we have a "fields" value that is either a string or None
        if fields is not None:
            mylog.info("Extracting (sampling: %s)", fields)
        verts = []
        samples = []
        for _io_chunk in parallel_objects(self.data_source.chunks([], "io")):
            for block, mask in self.data_source.blocks:
                my_verts = self._extract_isocontours_from_grid(
                    block,
                    self.surface_field,
                    self.field_value,
                    mask,
                    fields,
                    sample_type,
                    no_ghost=no_ghost,
                )
                if fields is not None:
                    my_verts, svals = my_verts
                    samples.append(svals)
                verts.append(my_verts)
        verts = np.concatenate(verts).transpose()
        verts = self.comm.par_combine_object(verts, op="cat", datatype="array")
        # verts is an ndarray here and will always be in code units, so we
        # expose it in the public API as a YTArray
        self._vertices = self.ds.arr(verts, "code_length")
        if fields is not None:
            samples = uconcatenate(samples)
            samples = self.comm.par_combine_object(samples, op="cat", datatype="array")
            if sample_type == "face":
                self[fields] = samples
            elif sample_type == "vertex":
                self.vertex_samples[fields] = samples

    def _extract_isocontours_from_grid(
        self,
        grid,
        field,
        value,
        mask,
        sample_values=None,
        sample_type="face",
        no_ghost=False,
    ):
        # TODO: check if multiple fields can be passed here
        vals = grid.get_vertex_centered_data([field], no_ghost=no_ghost)[field]
        if sample_values is not None:
            # TODO: is no_ghost=False correct here?
            svals = grid.get_vertex_centered_data([sample_values])[sample_values]
        else:
            svals = None

        sample_type = {"face": 1, "vertex": 2}[sample_type]
        my_verts = march_cubes_grid(
            value, vals, mask, grid.LeftEdge, grid.dds, svals, sample_type
        )
        return my_verts

    def calculate_flux(self, field_x, field_y, field_z, fluxing_field=None):
        r"""This calculates the flux over the surface.

        This function will conduct `Marching Cubes`_ on all the cells in a
        given data container (grid-by-grid), and then for each identified
        triangular segment of an isocontour in a given cell, calculate the
        gradient (i.e., normal) in the isocontoured field, interpolate the local
        value of the "fluxing" field, the area of the triangle, and then return:

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
        flux : YTQuantity
            The summed flux.

        References
        ----------

        .. _Marching Cubes:
            https://en.wikipedia.org/wiki/Marching_cubes

        Examples
        --------

        This will create a data object, find a nice value in the center, and
        calculate the metal flux over it.

        >>> sp = ds.sphere("max", (10, "kpc"))
        >>> surf = ds.surface(sp, ("gas", "density"), 5e-27)
        >>> flux = surf.calculate_flux(
        ...     ("gas", "velocity_x"),
        ...     ("gas", "velocity_y"),
        ...     ("gas", "velocity_z"),
        ...     ("gas", "metal_density"),
        ... )
        """
        flux = 0.0
        mylog.info("Fluxing %s", fluxing_field)
        for _io_chunk in parallel_objects(self.data_source.chunks([], "io")):
            for block, mask in self.data_source.blocks:
                flux += self._calculate_flux_in_grid(
                    block, mask, field_x, field_y, field_z, fluxing_field
                )
        flux = self.comm.mpi_allreduce(flux, op="sum")
        return flux

    def _calculate_flux_in_grid(
        self, grid, mask, field_x, field_y, field_z, fluxing_field=None
    ):

        vc_fields = [self.surface_field, field_x, field_y, field_z]
        if fluxing_field is not None:
            vc_fields.append(fluxing_field)

        vc_data = grid.get_vertex_centered_data(vc_fields)
        if fluxing_field is None:
            ff = self.ds.arr(
                np.ones_like(vc_data[self.surface_field].d, dtype="float64"),
                "dimensionless",
            )
        else:
            ff = vc_data[fluxing_field]
        surf_vals = vc_data[self.surface_field]
        field_x_vals = vc_data[field_x]
        field_y_vals = vc_data[field_y]
        field_z_vals = vc_data[field_z]
        ret = march_cubes_grid_flux(
            self.field_value,
            surf_vals,
            field_x_vals,
            field_y_vals,
            field_z_vals,
            ff,
            mask,
            grid.LeftEdge,
            grid.dds,
        )
        # assumes all the fluxing fields have the same units
        ret_units = field_x_vals.units * ff.units * grid.dds.units**2
        ret = self.ds.arr(ret, ret_units)
        ret.convert_to_units(self.ds.unit_system[ret_units.dimensions])
        return ret

    _vertices = None

    @property
    def vertices(self):
        if self._vertices is None:
            self.get_data()
        return self._vertices

    @property
    def triangles(self):
        vv = np.empty((self.vertices.shape[1] // 3, 3, 3), dtype="float64")
        vv = self.ds.arr(vv, self.vertices.units)
        for i in range(3):
            for j in range(3):
                vv[:, i, j] = self.vertices[j, i::3]
        return vv

    _surface_area = None

    @property
    def surface_area(self):
        if self._surface_area is not None:
            return self._surface_area
        tris = self.triangles
        x = tris[:, 1, :] - tris[:, 0, :]
        y = tris[:, 2, :] - tris[:, 0, :]
        areas = (x[:, 1] * y[:, 2] - x[:, 2] * y[:, 1]) ** 2
        np.add(areas, (x[:, 2] * y[:, 0] - x[:, 0] * y[:, 2]) ** 2, out=areas)
        np.add(areas, (x[:, 0] * y[:, 1] - x[:, 1] * y[:, 0]) ** 2, out=areas)
        np.sqrt(areas, out=areas)
        self._surface_area = 0.5 * areas.sum()
        return self._surface_area

    def export_obj(
        self,
        filename,
        transparency=1.0,
        dist_fac=None,
        color_field=None,
        emit_field=None,
        color_map=None,
        color_log=True,
        emit_log=True,
        plot_index=None,
        color_field_max=None,
        color_field_min=None,
        emit_field_max=None,
        emit_field_min=None,
    ):
        r"""Export the surface to the OBJ format

        Suitable for visualization in many different programs (e.g., Blender).
        NOTE: this exports an .obj file and an .mtl file, both with the general
        'filename' as a prefix.  The .obj file points to the .mtl file in its
        header, so if you move the 2 files, make sure you change the .obj header
        to account for this. ALSO NOTE: the emit_field needs to be a combination
        of the other 2 fields used to have the emissivity track with the color.

        Parameters
        ----------

        filename : string
            The file this will be exported to.  This cannot be a file-like
            object. If there are no file extensions included - both obj & mtl
            files are created.
        transparency : float
            This gives the transparency of the output surface plot.  Values
            from 0.0 (invisible) to 1.0 (opaque).
        dist_fac : float
            Divide the axes distances by this amount.
        color_field : string
            Should a field be sample and colormapped?
        emit_field : string
            Should we track the emissivity of a field? This should be a
            combination of the other 2 fields being used.
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
        >>> surf = ds.surface(sp, ("gas", "density"), 5e-27)
        >>> surf.export_obj("my_galaxy", transparency=trans)

        >>> sp = ds.sphere("max", (10, "kpc"))
        >>> mi, ma = sp.quantities.extrema("temperature")
        >>> rhos = [1e-24, 1e-25]
        >>> trans = [0.5, 1.0]
        >>> for i, r in enumerate(rhos):
        ...     surf = ds.surface(sp, "density", r)
        ...     surf.export_obj(
        ...         "my_galaxy",
        ...         transparency=trans[i],
        ...         color_field="temperature",
        ...         plot_index=i,
        ...         color_field_max=ma,
        ...         color_field_min=mi,
        ...     )

        >>> sp = ds.sphere("max", (10, "kpc"))
        >>> rhos = [1e-24, 1e-25]
        >>> trans = [0.5, 1.0]
        >>> def _Emissivity(field, data):
        ...     return (
        ...         data[("gas", "density")]
        ...         * data[("gas", "density")]
        ...         * np.sqrt(data[("gas", "temperature")])
        ...     )
        >>> ds.add_field(
        ...     ("gas", "emissivity"),
        ...     function=_Emissivity,
        ...     sampling_type="cell",
        ...     units=r"g**2*sqrt(K)/cm**6",
        ... )
        >>> for i, r in enumerate(rhos):
        ...     surf = ds.surface(sp, "density", r)
        ...     surf.export_obj(
        ...         "my_galaxy",
        ...         transparency=trans[i],
        ...         color_field="temperature",
        ...         emit_field="emissivity",
        ...         plot_index=i,
        ...     )

        """
        if color_map is None:
            color_map = ytcfg.get("yt", "default_colormap")
        if self.vertices is None:
            if color_field is not None:
                self.get_data(color_field, "face")
        elif color_field is not None:
            if color_field not in self.field_data:
                self[color_field]
        if color_field is None:
            self.get_data(self.surface_field, "face")
        if emit_field is not None:
            if color_field not in self.field_data:
                self[emit_field]
        only_on_root(
            self._export_obj,
            filename,
            transparency,
            dist_fac,
            color_field,
            emit_field,
            color_map,
            color_log,
            emit_log,
            plot_index,
            color_field_max,
            color_field_min,
            emit_field_max,
            emit_field_min,
        )

    def _color_samples_obj(
        self,
        cs,
        em,
        color_log,
        emit_log,
        color_map,
        arr,
        color_field_max,
        color_field_min,
        color_field,
        emit_field_max,
        emit_field_min,
        emit_field,
    ):  # this now holds for obj files
        if color_field is not None:
            if color_log:
                cs = np.log10(cs)
        if emit_field is not None:
            if emit_log:
                em = np.log10(em)
        if color_field is not None:
            if color_field_min is None:
                cs = [float(field) for field in cs]
                cs = np.array(cs)
                mi = cs.min()
            else:
                mi = color_field_min
                if color_log:
                    mi = np.log10(mi)
            if color_field_max is None:
                cs = [float(field) for field in cs]
                cs = np.array(cs)
                ma = cs.max()
            else:
                ma = color_field_max
                if color_log:
                    ma = np.log10(ma)
            cs = (cs - mi) / (ma - mi)
        else:
            cs[:] = 1.0
        # to get color indices for OBJ formatting
        lut = get_colormap_lut(color_map)

        x = np.mgrid[0.0 : 1.0 : lut[0].shape[0] * 1j]
        arr["cind"][:] = (np.interp(cs, x, x) * (lut[0].shape[0] - 1)).astype("uint8")
        # now, get emission
        if emit_field is not None:
            if emit_field_min is None:
                em = [float(field) for field in em]
                em = np.array(em)
                emi = em.min()
            else:
                emi = emit_field_min
                if emit_log:
                    emi = np.log10(emi)
            if emit_field_max is None:
                em = [float(field) for field in em]
                em = np.array(em)
                ema = em.max()
            else:
                ema = emit_field_max
                if emit_log:
                    ema = np.log10(ema)
            em = (em - emi) / (ema - emi)
            x = np.mgrid[0.0:255.0:2j]  # assume 1 emissivity per color
            arr["emit"][:] = (
                np.interp(em, x, x)
            ) * 2.0  # for some reason, max emiss = 2
        else:
            arr["emit"][:] = 0.0

    @parallel_root_only
    def _export_obj(
        self,
        filename,
        transparency,
        dist_fac=None,
        color_field=None,
        emit_field=None,
        color_map=None,
        color_log=True,
        emit_log=True,
        plot_index=None,
        color_field_max=None,
        color_field_min=None,
        emit_field_max=None,
        emit_field_min=None,
    ):
        if color_map is None:
            color_map = ytcfg.get("yt", "default_colormap")
        if plot_index is None:
            plot_index = 0
        if isinstance(filename, io.IOBase):
            fobj = filename + ".obj"
            fmtl = filename + ".mtl"
        else:
            if plot_index == 0:
                fobj = open(filename + ".obj", "w")
                fmtl = open(filename + ".mtl", "w")
                cc = 1
            else:
                # read in last vertex
                linesave = ""
                for line in fileinput.input(filename + ".obj"):
                    if line[0] == "f":
                        linesave = line
                p = [m.start() for m in finditer(" ", linesave)]
                cc = int(linesave[p[len(p) - 1] :]) + 1
                fobj = open(filename + ".obj", "a")
                fmtl = open(filename + ".mtl", "a")
        ftype = [("cind", "uint8"), ("emit", "float")]
        vtype = [("x", "float"), ("y", "float"), ("z", "float")]
        if plot_index == 0:
            fobj.write("# yt OBJ file\n")
            fobj.write("# www.yt-project.org\n")
            fobj.write(
                f"mtllib {filename}.mtl\n\n"
            )  # use this material file for the faces
            fmtl.write("# yt MLT file\n")
            fmtl.write("# www.yt-project.org\n\n")
        # (0) formulate vertices
        nv = self.vertices.shape[1]  # number of groups of vertices
        f = np.empty(
            nv // self.vertices.shape[0], dtype=ftype
        )  # store sets of face colors
        v = np.empty(nv, dtype=vtype)  # stores vertices
        if color_field is not None:
            cs = self[color_field]
        else:
            cs = np.empty(self.vertices.shape[1] // self.vertices.shape[0])
        if emit_field is not None:
            em = self[emit_field]
        else:
            em = np.empty(self.vertices.shape[1] // self.vertices.shape[0])
        self._color_samples_obj(
            cs,
            em,
            color_log,
            emit_log,
            color_map,
            f,
            color_field_max,
            color_field_min,
            color_field,
            emit_field_max,
            emit_field_min,
            emit_field,
        )  # map color values to color scheme

        lut = get_colormap_lut(color_map)

        # interpolate emissivity to enumerated colors
        emiss = np.interp(
            np.mgrid[0 : lut[0].shape[0]], np.mgrid[0 : len(cs)], f["emit"][:]
        )
        if dist_fac is None:  # then normalize by bounds
            DLE = self.pf.domain_left_edge
            DRE = self.pf.domain_right_edge
            bounds = [(DLE[i], DRE[i]) for i in range(3)]
            for i, ax in enumerate("xyz"):
                # Do the bounds first since we cast to f32
                tmp = self.vertices[i, :]
                np.subtract(tmp, bounds[i][0], tmp)
                w = bounds[i][1] - bounds[i][0]
                np.divide(tmp, w, tmp)
                np.subtract(tmp, 0.5, tmp)  # Center at origin.
                v[ax][:] = tmp
        else:
            for i, ax in enumerate("xyz"):
                tmp = self.vertices[i, :]
                np.divide(tmp, dist_fac, tmp)
                v[ax][:] = tmp
        # (1) write all colors per surface to mtl file
        for i in range(0, lut[0].shape[0]):
            omname = f"material_{i}_{plot_index}"  # name of the material
            fmtl.write(
                f"newmtl {omname}\n"
            )  # the specific material (color) for this face
            fmtl.write(f"Ka {0.0:.6f} {0.0:.6f} {0.0:.6f}\n")  # ambient color, keep off
            fmtl.write(
                f"Kd {lut[0][i]:.6f} {lut[1][i]:.6f} {lut[2][i]:.6f}\n"
            )  # color of face
            fmtl.write(
                f"Ks {0.0:.6f} {0.0:.6f} {0.0:.6f}\n"
            )  # specular color, keep off
            fmtl.write(f"d {transparency:.6f}\n")  # transparency
            fmtl.write(f"em {emiss[i]:.6f}\n")  # emissivity per color
            fmtl.write("illum 2\n")  # not relevant, 2 means highlights on?
            fmtl.write(f"Ns {0:.6f}\n\n")  # keep off, some other specular thing
        # (2) write vertices
        for i in range(0, self.vertices.shape[1]):
            fobj.write(f"v {v['x'][i]:.6f} {v['y'][i]:.6f} {v['z'][i]:.6f}\n")
        fobj.write("#done defining vertices\n\n")
        # (3) define faces and materials for each face
        for i in range(0, self.triangles.shape[0]):
            omname = f"material_{f['cind'][i]}_{plot_index}"  # which color to use
            fobj.write(
                f"usemtl {omname}\n"
            )  # which material to use for this face (color)
            fobj.write(f"f {cc} {cc+1} {cc+2}\n\n")  # vertices to color
            cc = cc + 3
        fmtl.close()
        fobj.close()

    def export_blender(
        self,
        transparency=1.0,
        dist_fac=None,
        color_field=None,
        emit_field=None,
        color_map=None,
        color_log=True,
        emit_log=True,
        plot_index=None,
        color_field_max=None,
        color_field_min=None,
        emit_field_max=None,
        emit_field_min=None,
    ):
        r"""This exports the surface to the OBJ format, suitable for visualization
        in many different programs (e.g., Blender).  NOTE: this exports an .obj file
        and an .mtl file, both with the general 'filename' as a prefix.
        The .obj file points to the .mtl file in its header, so if you move the 2
        files, make sure you change the .obj header to account for this. ALSO NOTE:
        the emit_field needs to be a combination of the other 2 fields used to
        have the emissivity track with the color.

        Parameters
        ----------
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
        >>> surf = ds.surface(sp, ("gas", "density"), 5e-27)
        >>> surf.export_obj("my_galaxy", transparency=trans)

        >>> sp = ds.sphere("max", (10, "kpc"))
        >>> mi, ma = sp.quantities.extrema("temperature")[0]
        >>> rhos = [1e-24, 1e-25]
        >>> trans = [0.5, 1.0]
        >>> for i, r in enumerate(rhos):
        ...     surf = ds.surface(sp, "density", r)
        ...     surf.export_obj(
        ...         "my_galaxy",
        ...         transparency=trans[i],
        ...         color_field="temperature",
        ...         plot_index=i,
        ...         color_field_max=ma,
        ...         color_field_min=mi,
        ...     )

        >>> sp = ds.sphere("max", (10, "kpc"))
        >>> rhos = [1e-24, 1e-25]
        >>> trans = [0.5, 1.0]
        >>> def _Emissivity(field, data):
        ...     return (
        ...         data[("gas", "density")]
        ...         * data[("gas", "density")]
        ...         * np.sqrt(data[("gas", "temperature")])
        ...     )
        >>> ds.add_field(("gas", "emissivity"), function=_Emissivity, units="g / cm**6")
        >>> for i, r in enumerate(rhos):
        ...     surf = ds.surface(sp, "density", r)
        ...     surf.export_obj(
        ...         "my_galaxy",
        ...         transparency=trans[i],
        ...         color_field="temperature",
        ...         emit_field="emissivity",
        ...         plot_index=i,
        ...     )

        """
        if color_map is None:
            color_map = ytcfg.get("yt", "default_colormap")
        if self.vertices is None:
            if color_field is not None:
                self.get_data(color_field, "face")
        elif color_field is not None:
            if color_field not in self.field_data:
                self[color_field]
        if color_field is None:
            self.get_data(self.surface_field, "face")
        if emit_field is not None:
            if color_field not in self.field_data:
                self[emit_field]
        fullverts, colors, alpha, emisses, colorindex = only_on_root(
            self._export_blender,
            transparency,
            dist_fac,
            color_field,
            emit_field,
            color_map,
            color_log,
            emit_log,
            plot_index,
            color_field_max,
            color_field_min,
            emit_field_max,
            emit_field_min,
        )
        return fullverts, colors, alpha, emisses, colorindex

    def _export_blender(
        self,
        transparency,
        dist_fac=None,
        color_field=None,
        emit_field=None,
        color_map=None,
        color_log=True,
        emit_log=True,
        plot_index=None,
        color_field_max=None,
        color_field_min=None,
        emit_field_max=None,
        emit_field_min=None,
    ):
        if color_map is None:
            color_map = ytcfg.get("yt", "default_colormap")
        if plot_index is None:
            plot_index = 0
        ftype = [("cind", "uint8"), ("emit", "float")]
        vtype = [("x", "float"), ("y", "float"), ("z", "float")]
        # (0) formulate vertices
        nv = self.vertices.shape[1]  # number of groups of vertices
        f = np.empty(
            nv / self.vertices.shape[0], dtype=ftype
        )  # store sets of face colors
        v = np.empty(nv, dtype=vtype)  # stores vertices
        if color_field is not None:
            cs = self[color_field]
        else:
            cs = np.empty(self.vertices.shape[1] / self.vertices.shape[0])
        if emit_field is not None:
            em = self[emit_field]
        else:
            em = np.empty(self.vertices.shape[1] / self.vertices.shape[0])
        self._color_samples_obj(
            cs,
            em,
            color_log,
            emit_log,
            color_map,
            f,
            color_field_max,
            color_field_min,
            color_field,
            emit_field_max,
            emit_field_min,
            emit_field,
        )  # map color values to color scheme

        lut = get_colormap_lut(color_map)

        # interpolate emissivity to enumerated colors
        emiss = np.interp(
            np.mgrid[0 : lut[0].shape[0]], np.mgrid[0 : len(cs)], f["emit"][:]
        )
        if dist_fac is None:  # then normalize by bounds
            DLE = self.ds.domain_left_edge
            DRE = self.ds.domain_right_edge
            bounds = [(DLE[i], DRE[i]) for i in range(3)]
            for i, ax in enumerate("xyz"):
                # Do the bounds first since we cast to f32
                tmp = self.vertices[i, :]
                np.subtract(tmp, bounds[i][0], tmp)
                w = bounds[i][1] - bounds[i][0]
                np.divide(tmp, w, tmp)
                np.subtract(tmp, 0.5, tmp)  # Center at origin.
                v[ax][:] = tmp
        else:
            for i, ax in enumerate("xyz"):
                tmp = self.vertices[i, :]
                np.divide(tmp, dist_fac, tmp)
                v[ax][:] = tmp
        return v, lut, transparency, emiss, f["cind"]

    def export_ply(
        self,
        filename,
        bounds=None,
        color_field=None,
        color_map=None,
        color_log=True,
        sample_type="face",
        no_ghost=False,
    ):
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

        >>> from yt.units import kpc
        >>> sp = ds.sphere("max", (10, "kpc"))
        >>> surf = ds.surface(sp, ("gas", "density"), 5e-27)
        >>> print(surf[("gas", "temperature")])
        >>> print(surf.vertices)
        >>> bounds = [
        ...     (sp.center[i] - 5.0 * kpc, sp.center[i] + 5.0 * kpc) for i in range(3)
        ... ]
        >>> surf.export_ply("my_galaxy.ply", bounds=bounds)
        """
        if color_map is None:
            color_map = ytcfg.get("yt", "default_colormap")
        if self.vertices is None:
            self.get_data(color_field, sample_type, no_ghost=no_ghost)
        elif color_field is not None:
            if sample_type == "face" and color_field not in self.field_data:
                self[color_field]
            elif sample_type == "vertex" and color_field not in self.vertex_samples:
                self.get_data(color_field, sample_type, no_ghost=no_ghost)
        self._export_ply(
            filename, bounds, color_field, color_map, color_log, sample_type
        )

    def _color_samples(self, cs, color_log, color_map, arr):
        if color_log:
            cs = np.log10(cs)
        mi, ma = cs.min(), cs.max()
        cs = (cs - mi) / (ma - mi)
        from yt.visualization.image_writer import map_to_colors

        cs = map_to_colors(cs, color_map)
        arr["red"][:] = cs[0, :, 0]
        arr["green"][:] = cs[0, :, 1]
        arr["blue"][:] = cs[0, :, 2]

    @parallel_root_only
    def _export_ply(
        self,
        filename,
        bounds=None,
        color_field=None,
        color_map=None,
        color_log=True,
        sample_type="face",
    ):
        if color_map is None:
            color_map = ytcfg.get("yt", "default_colormap")
        if hasattr(filename, "read"):
            f = filename
        else:
            f = open(filename, "wb")
        if bounds is None:
            DLE = self.ds.domain_left_edge
            DRE = self.ds.domain_right_edge
            bounds = [(DLE[i], DRE[i]) for i in range(3)]
        elif any([not all([isinstance(be, YTArray) for be in b]) for b in bounds]):
            bounds = [
                tuple(
                    be if isinstance(be, YTArray) else self.ds.quan(be, "code_length")
                    for be in b
                )
                for b in bounds
            ]
        nv = self.vertices.shape[1]
        vs = [
            ("x", "<f"),
            ("y", "<f"),
            ("z", "<f"),
            ("red", "uint8"),
            ("green", "uint8"),
            ("blue", "uint8"),
        ]
        fs = [
            ("ni", "uint8"),
            ("v1", "<i4"),
            ("v2", "<i4"),
            ("v3", "<i4"),
            ("red", "uint8"),
            ("green", "uint8"),
            ("blue", "uint8"),
        ]
        f.write(b"ply\n")
        f.write(b"format binary_little_endian 1.0\n")
        line = "element vertex %i\n" % (nv)
        f.write(line.encode("latin-1"))
        f.write(b"property float x\n")
        f.write(b"property float y\n")
        f.write(b"property float z\n")
        if color_field is not None and sample_type == "vertex":
            f.write(b"property uchar red\n")
            f.write(b"property uchar green\n")
            f.write(b"property uchar blue\n")
            v = np.empty(self.vertices.shape[1], dtype=vs)
            cs = self.vertex_samples[color_field]
            self._color_samples(cs, color_log, color_map, v)
        else:
            v = np.empty(self.vertices.shape[1], dtype=vs[:3])
        line = "element face %i\n" % (nv / 3)
        f.write(line.encode("latin-1"))
        f.write(b"property list uchar int vertex_indices\n")
        if color_field is not None and sample_type == "face":
            f.write(b"property uchar red\n")
            f.write(b"property uchar green\n")
            f.write(b"property uchar blue\n")
            # Now we get our samples
            cs = self[color_field]
            arr = np.empty(cs.shape[0], dtype=np.dtype(fs))
            self._color_samples(cs, color_log, color_map, arr)
        else:
            arr = np.empty(nv // 3, np.dtype(fs[:-3]))
        for i, ax in enumerate("xyz"):
            # Do the bounds first since we cast to f32
            tmp = self.vertices[i, :]
            np.subtract(tmp, bounds[i][0], tmp)
            w = bounds[i][1] - bounds[i][0]
            np.divide(tmp, w, tmp)
            np.subtract(tmp, 0.5, tmp)  # Center at origin.
            v[ax][:] = tmp
        f.write(b"end_header\n")
        v.tofile(f)
        arr["ni"][:] = 3
        vi = np.arange(nv, dtype="<i")
        vi.shape = (nv // 3, 3)
        arr["v1"][:] = vi[:, 0]
        arr["v2"][:] = vi[:, 1]
        arr["v3"][:] = vi[:, 2]
        arr.tofile(f)
        if filename is not f:
            f.close()

    def export_sketchfab(
        self,
        title,
        description,
        api_key=None,
        color_field=None,
        color_map=None,
        color_log=True,
        bounds=None,
        no_ghost=False,
    ):
        r"""This exports Surfaces to SketchFab.com, where they can be viewed
        interactively in a web browser.

        SketchFab.com is a proprietary web service that provides WebGL
        rendering of models.  This routine will use temporary files to
        construct a compressed binary representation (in .PLY format) of the
        Surface and any optional fields you specify and upload it to
        SketchFab.com.  It requires an API key, which can be found on your
        SketchFab.com dashboard.  You can either supply the API key to this
        routine directly or you can place it in the variable
        "sketchfab_api_key" in your ~/.config/yt/yt.toml file.  This function is
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

        >>> from yt.units import kpc
        >>> dd = ds.sphere("max", (200, "kpc"))
        >>> rho = 5e-27
        >>> bounds = [
        ...     (dd.center[i] - 100.0 * kpc, dd.center[i] + 100.0 * kpc)
        ...     for i in range(3)
        ... ]
        ...
        >>> surf = ds.surface(dd, ("gas", "density"), rho)
        >>> rv = surf.export_sketchfab(
        ...     title="Testing Upload",
        ...     description="A simple test of the uploader",
        ...     color_field="temperature",
        ...     color_map="hot",
        ...     color_log=True,
        ...     bounds=bounds,
        ... )
        ...
        """
        if color_map is None:
            color_map = ytcfg.get("yt", "default_colormap")
        api_key = api_key or ytcfg.get("yt", "sketchfab_api_key")
        if api_key in (None, "None"):
            raise YTNoAPIKey("SketchFab.com", "sketchfab_api_key")

        ply_file = TemporaryFile()
        self.export_ply(
            ply_file,
            bounds,
            color_field,
            color_map,
            color_log,
            sample_type="vertex",
            no_ghost=no_ghost,
        )
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

        zfs = NamedTemporaryFile(suffix=".zip")
        with zipfile.ZipFile(zfs, "w", zipfile.ZIP_DEFLATED) as zf:
            zf.writestr("yt_export.ply", ply_file.read())
        zfs.seek(0)

        zfs.seek(0)
        data = {
            "token": api_key,
            "name": title,
            "description": description,
            "tags": "yt",
        }
        files = {"modelFile": zfs}
        upload_id = self._upload_to_sketchfab(data, files)
        upload_id = self.comm.mpi_bcast(upload_id, root=0)
        return upload_id

    @parallel_root_only
    def _upload_to_sketchfab(self, data, files):
        from yt.utilities.on_demand_imports import _requests as requests

        SKETCHFAB_DOMAIN = "sketchfab.com"
        SKETCHFAB_API_URL = f"https://api.{SKETCHFAB_DOMAIN}/v2/models"
        SKETCHFAB_MODEL_URL = f"https://{SKETCHFAB_DOMAIN}/models/"

        try:
            r = requests.post(SKETCHFAB_API_URL, data=data, files=files, verify=False)
        except requests.exceptions.RequestException:
            mylog.exception("An error has occurred")
            return

        result = r.json()

        if r.status_code != requests.codes.created:
            mylog.error("Upload to SketchFab failed with error: %s", result)
            return

        model_uid = result["uid"]
        model_url = SKETCHFAB_MODEL_URL + model_uid
        if model_uid:
            mylog.info("Model uploaded to: %s", model_url)
        else:
            mylog.error("Problem uploading.")

        return model_uid


class YTOctree(YTSelectionContainer3D):
    """A 3D region with all the data filled into an octree. This container
    will deposit particle fields onto octs using a kernel and SPH smoothing.
    The octree is built in a depth-first fashion. Depth-first search (DFS)
    means that tree starts refining at the root node (this is the largest node
    which contains every particles) and refines as far as possible along each
    branch before backtracking.

    Parameters
    ----------
    right_edge : array_like
        The right edge of the region to be extracted.  Specify units by supplying
        a YTArray, otherwise code length units are assumed.
    left_edge : array_like
        The left edge of the region to be extracted.  Specify units by supplying
        a YTArray, otherwise code length units are assumed.
    n_ref: int
        This is the maximum number of particles per leaf in the resulting
        octree.
    ptypes: list
        This is the type of particles to include when building the tree. This
        will default to all particles.

    Examples
    --------

    octree = ds.octree(n_ref=64)
    x_positions_of_cells = octree[('index', 'x')]
    y_positions_of_cells = octree[('index', 'y')]
    z_positions_of_cells = octree[('index', 'z')]
    density_of_gas_in_cells = octree[('gas', 'density')]
    """

    _spatial = True
    _type_name = "octree"

    _con_args = ("left_edge", "right_edge", "n_ref")
    _container_fields = (
        ("index", "dx"),
        ("index", "dy"),
        ("index", "dz"),
        ("index", "x"),
        ("index", "y"),
        ("index", "z"),
        ("index", "depth"),
        ("index", "refined"),
        ("index", "sizes"),
        ("index", "positions"),
    )

    def __init__(
        self,
        left_edge=None,
        right_edge=None,
        n_ref=32,
        ptypes=None,
        ds=None,
        field_parameters=None,
    ):
        if field_parameters is None:
            center = None
        else:
            center = field_parameters.get("center", None)
        YTSelectionContainer3D.__init__(self, center, ds, field_parameters)

        self.left_edge = self._sanitize_edge(left_edge, ds.domain_left_edge)
        self.right_edge = self._sanitize_edge(right_edge, ds.domain_right_edge)
        self.n_ref = n_ref
        self.ptypes = self._sanitize_ptypes(ptypes)

        self._setup_data_source()
        self.tree

    def _generate_tree(self):
        positions = []
        for ptype in self.ptypes:
            positions.append(
                self._data_source[(ptype, "particle_position")]
                .in_units("code_length")
                .d
            )

        positions = (
            np.concatenate(positions) if len(positions) > 0 else np.array(positions)
        )
        if not positions.size:
            mylog.info("No particles found!")
            self._octree = None
            return

        mylog.info("Allocating Octree for %s particles", positions.shape[0])
        self._octree = CyOctree(
            positions,
            left_edge=self.left_edge.to("code_length").d,
            right_edge=self.right_edge.to("code_length").d,
            n_ref=self.n_ref,
        )
        mylog.info("Allocated %s nodes in octree", self._octree.num_nodes)
        mylog.info("Octree bound %s particles", self._octree.bound_particles)

        # Now we store the index data about the octree in the python container
        ds = self.ds
        pos = ds.arr(self._octree.node_positions, "code_length")
        self[("index", "positions")] = pos
        self[("index", "x")] = pos[:, 0]
        self[("index", "y")] = pos[:, 1]
        self[("index", "z")] = pos[:, 2]
        self[("index", "refined")] = self._octree.node_refined

        sizes = ds.arr(self._octree.node_sizes, "code_length")
        self[("index", "sizes")] = sizes
        self[("index", "dx")] = sizes[:, 0]
        self[("index", "dy")] = sizes[:, 1]
        self[("index", "dz")] = sizes[:, 2]
        self[("index", "depth")] = self._octree.node_depth

    @property
    def tree(self):
        """
        The Cython+Python octree instance
        """
        if hasattr(self, "_octree"):
            return self._octree

        self._generate_tree()
        return self._octree

    @property
    def sph_smoothing_style(self):
        smoothing_style = getattr(self.ds, "sph_smoothing_style", "scatter")
        return smoothing_style

    @property
    def sph_normalize(self):
        normalize = getattr(self.ds, "use_sph_normalization", "normalize")
        return normalize

    @property
    def sph_num_neighbors(self):
        num_neighbors = getattr(self.ds, "num_neighbors", 32)
        return num_neighbors

    def _sanitize_ptypes(self, ptypes):
        if ptypes is None:
            return ["all"]

        if not isinstance(ptypes, list):
            ptypes = [ptypes]

        for ptype in ptypes:
            if ptype not in self.ds.particle_types:
                mess = f"{ptype} not found. Particle type must "
                mess += "be in the dataset!"
                raise TypeError(mess)

        return ptypes

    def _setup_data_source(self):
        mylog.info(
            (
                "Allocating octree with spatial range "
                "[%.4e, %.4e, %.4e] code_length to "
                "[%.4e, %.4e, %.4e] code_length"
            ),
            *self.left_edge.to("code_length").d,
            *self.right_edge.to("code_length").d,
        )
        self._data_source = self.ds.region(self.center, self.left_edge, self.right_edge)

    def _sanitize_edge(self, edge, default):
        if edge is None:
            return default.copy()
        if not is_sequence(edge):
            edge = [edge] * len(self.ds.domain_left_edge)
        if len(edge) != len(self.ds.domain_left_edge):
            raise RuntimeError(
                "Length of edges must match the dimensionality of the dataset"
            )
        if hasattr(edge, "units"):
            if edge.units.registry is self.ds.unit_registry:
                return edge
            edge_units = edge.units
        else:
            edge_units = "code_length"
        return self.ds.arr(edge, edge_units)

    def get_data(self, fields=None):
        if fields is None:
            return

        if hasattr(self.ds, "_sph_ptypes"):
            sph_ptypes = self.ds._sph_ptypes
        else:
            sph_ptypes = tuple(
                value
                for value in self.ds.particle_types_raw
                if value in ["PartType0", "Gas", "gas", "io"]
            )

        if len(sph_ptypes) == 0:
            raise RuntimeError
        smoothing_style = self.sph_smoothing_style
        normalize = self.sph_normalize

        if fields[0] in sph_ptypes:
            units = self.ds._get_field_info(fields).units
            if smoothing_style == "scatter":
                self._scatter_smooth(fields, units, normalize)
            else:
                self._gather_smooth(fields, units, normalize)
        elif fields[0] == "index":
            return self[fields]
        else:
            raise NotImplementedError

    def _gather_smooth(self, fields, units, normalize):
        buff = np.zeros(self.tree.num_nodes, dtype="float64")

        # Again, attempt to load num_neighbors from the octree
        num_neighbors = self.sph_num_neighbors

        # For the gather approach we load up all of the data, this like other
        # gather approaches is not memory conservative but with spatial chunking
        # this can be fixed
        fields_to_get = [
            "particle_position",
            "density",
            "particle_mass",
            "smoothing_length",
            fields[1],
        ]
        all_fields = all_data(self.ds, fields[0], fields_to_get, kdtree=True)

        interpolate_sph_positions_gather(
            buff,
            all_fields["particle_position"],
            self[("index", "positions")],
            all_fields["smoothing_length"],
            all_fields["particle_mass"],
            all_fields["density"],
            all_fields[fields[1]].in_units(units),
            self.ds.index.kdtree,
            use_normalization=normalize,
            num_neigh=num_neighbors,
        )

        self[fields] = self.ds.arr(buff[~self[("index", "refined")]], units)

    def _scatter_smooth(self, fields, units, normalize):
        buff = np.zeros(self.tree.num_nodes, dtype="float64")

        if normalize:
            buff_den = np.zeros(buff.shape[0], dtype="float64")
        else:
            buff_den = np.empty(0)

        ptype = fields[0]
        pbar = tqdm(desc=f"Interpolating (scatter) SPH field {fields[0]}")
        for chunk in self._data_source.chunks([fields], "io"):
            px = chunk[(ptype, "particle_position_x")].to("code_length").d
            py = chunk[(ptype, "particle_position_y")].to("code_length").d
            pz = chunk[(ptype, "particle_position_z")].to("code_length").d
            hsml = chunk[(ptype, "smoothing_length")].to("code_length").d
            pmass = chunk[(ptype, "particle_mass")].to("code_mass").d
            pdens = chunk[(ptype, "density")].to("code_mass/code_length**3").d
            field_quantity = chunk[fields].to(units).d

            if px.shape[0] > 0:
                self.tree.interpolate_sph_cells(
                    buff,
                    buff_den,
                    px,
                    py,
                    pz,
                    pmass,
                    pdens,
                    hsml,
                    field_quantity,
                    use_normalization=normalize,
                )

            pbar.update(1)
        pbar.close()

        if normalize:
            normalization_1d_utility(buff, buff_den)

        self[fields] = self.ds.arr(buff[~self[("index", "refined")]], units)
