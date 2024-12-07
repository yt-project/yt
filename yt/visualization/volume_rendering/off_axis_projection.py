import numpy as np

from yt.data_objects.api import ImageArray
from yt.funcs import is_sequence, mylog
from yt.geometry.oct_geometry_handler import OctreeIndex
from yt.units.unit_object import Unit  # type: ignore
from yt.utilities.lib.image_utilities import add_cells_to_image_offaxis
from yt.utilities.lib.partitioned_grid import PartitionedGrid
from yt.utilities.lib.pixelization_routines import (
    normalization_2d_utility,
    off_axis_projection_SPH,
)
from yt.visualization.volume_rendering.lens import PlaneParallelLens

from .render_source import KDTreeVolumeSource
from .scene import Scene
from .transfer_functions import ProjectionTransferFunction
from .utils import data_source_or_all


def off_axis_projection(
    data_source,
    center,
    normal_vector,
    width,
    resolution,
    item,
    weight=None,
    volume=None,
    no_ghost=False,
    interpolated=False,
    north_vector=None,
    depth=None,
    num_threads=1,
    method="integrate",
):
    r"""Project through a dataset, off-axis, and return the image plane.

    This function will accept the necessary items to integrate through a volume
    at an arbitrary angle and return the integrated field of view to the user.
    Note that if a weight is supplied, it will multiply the pre-interpolated
    values together, then create cell-centered values, then interpolate within
    the cell to conduct the integration.

    Parameters
    ----------
    data_source : ~yt.data_objects.static_output.Dataset
                  or ~yt.data_objects.data_containers.YTSelectionDataContainer
        This is the dataset or data object to volume render.
    center : array_like
        The current 'center' of the view port -- the focal point for the
        camera.
    normal_vector : array_like
        The vector between the camera position and the center.
    width : float or list of floats
        The current width of the image.  If a single float, the volume is
        cubical, but if not, it is left/right, top/bottom, front/back
    resolution : int or list of ints
        The number of pixels in each direction.
    item: string
        The field to project through the volume
    weight : optional, default None
        If supplied, the field will be pre-multiplied by this, then divided by
        the integrated value of this field.  This returns an average rather
        than a sum.
    volume : `yt.extensions.volume_rendering.AMRKDTree`, optional
        The volume to ray cast through.  Can be specified for finer-grained
        control, but otherwise will be automatically generated.
    no_ghost: bool, optional
        Optimization option.  If True, homogenized bricks will
        extrapolate out from grid instead of interpolating from
        ghost zones that have to first be calculated.  This can
        lead to large speed improvements, but at a loss of
        accuracy/smoothness in resulting image.  The effects are
        less notable when the transfer function is smooth and
        broad. Default: True
    interpolated : optional, default False
        If True, the data is first interpolated to vertex-centered data,
        then tri-linearly interpolated along the ray. Not suggested for
        quantitative studies.
    north_vector : optional, array_like, default None
        A vector that, if specified, restricts the orientation such that the
        north vector dotted into the image plane points "up". Useful for rotations
    depth: float, tuple[float, str], or unyt_array of size 1.
        specify the depth of the projection region (size along the
        line of sight). If no units are given (unyt_array or second
        tuple element), code units are assumed.
    num_threads: integer, optional, default 1
        Use this many OpenMP threads during projection.
    method : string
        The method of projection.  Valid methods are:

        "integrate" with no weight_field specified : integrate the requested
        field along the line of sight.

        "integrate" with a weight_field specified : weight the requested
        field by the weighting field and integrate along the line of sight.

        "sum" : This method is the same as integrate, except that it does not
        multiply by a path length when performing the integration, and is
        just a straight summation of the field along the given axis. WARNING:
        This should only be used for uniform resolution grid datasets, as other
        datasets may result in unphysical images.
        or camera movements.
    Returns
    -------
    image : array
        An (N,N) array of the final integrated values, in float64 form.

    Examples
    --------

    >>> image = off_axis_projection(
    ...     ds,
    ...     [0.5, 0.5, 0.5],
    ...     [0.2, 0.3, 0.4],
    ...     0.2,
    ...     N,
    ...     ("gas", "temperature"),
    ...     ("gas", "density"),
    ... )
    >>> write_image(np.log10(image), "offaxis.png")

    """
    if method not in ("integrate", "sum"):
        raise NotImplementedError(
            "Only 'integrate' or 'sum' methods are valid for off-axis-projections"
        )

    if interpolated:
        raise NotImplementedError(
            "Only interpolated=False methods are currently implemented "
            "for off-axis-projections"
        )

    data_source = data_source_or_all(data_source)

    item = data_source._determine_fields([item])[0]

    # Assure vectors are numpy arrays as expected by cython code
    normal_vector = np.array(normal_vector, dtype="float64")
    if north_vector is not None:
        north_vector = np.array(north_vector, dtype="float64")
    # Add the normal as a field parameter to the data source
    # so line of sight fields can use it
    data_source.set_field_parameter("axis", normal_vector)

    # Sanitize units
    if not hasattr(center, "units"):
        center = data_source.ds.arr(center, "code_length")
    if not hasattr(width, "units"):
        width = data_source.ds.arr(width, "code_length")
    if depth is not None:
        # handle units (intrinsic or as a tuple),
        # then convert to code length
        # float -> assumed to be in code units
        if isinstance(depth, tuple):
            depth = data_source.ds.arr(np.array([depth[0]]), depth[1])
        if hasattr(depth, "units"):
            depth = depth.to("code_length").d

        # depth = data_source.ds.arr(depth, "code_length")

    if hasattr(data_source.ds, "_sph_ptypes"):
        if method != "integrate":
            raise NotImplementedError("SPH Only allows 'integrate' method")

        sph_ptypes = data_source.ds._sph_ptypes
        fi = data_source.ds.field_info[item]

        raise_error = False

        ptype = sph_ptypes[0]
        ppos = [f"particle_position_{ax}" for ax in "xyz"]
        # Assure that the field we're trying to off-axis project
        # has a field type as the SPH particle type or if the field is an
        # alias to an SPH field or is a 'gas' field
        if item[0] in data_source.ds.known_filters:
            if item[0] not in sph_ptypes:
                raise_error = True
            else:
                ptype = item[0]
                ppos = ["x", "y", "z"]
        elif fi.is_alias:
            if fi.alias_name[0] not in sph_ptypes:
                raise_error = True
            elif item[0] != "gas":
                ptype = item[0]
        else:
            if fi.name[0] not in sph_ptypes and fi.name[0] != "gas":
                raise_error = True

        if raise_error:
            raise RuntimeError(
                "Can only perform off-axis projections for SPH fields, "
                f"Received {item!r}"
            )

        normal = np.array(normal_vector)
        normal = normal / np.linalg.norm(normal)

        # If north_vector is None, we set the default here.
        # This is chosen so that if normal_vector is one of the
        # cartesian coordinate axes, the projection will match
        # the corresponding on-axis projection.
        if north_vector is None:
            vecs = np.identity(3)
            t = np.cross(vecs, normal).sum(axis=1)
            ax = t.argmax()
            east_vector = np.cross(vecs[ax, :], normal).ravel()
            north = np.cross(normal, east_vector).ravel()
        else:
            north = np.array(north_vector)
            north = north / np.linalg.norm(north)
            east_vector = np.cross(north, normal).ravel()

        # if weight is None:
        buf = np.zeros((resolution[0], resolution[1]), dtype="float64")
        mask = np.ones_like(buf, dtype="uint8")

        ## width from fixed_resolution.py is just the size of the domain
        # x_min = center[0] - width[0] / 2
        # x_max = center[0] + width[0] / 2
        # y_min = center[1] - width[1] / 2
        # y_max = center[1] + width[1] / 2
        # z_min = center[2] - width[2] / 2
        # z_max = center[2] + width[2] / 2

        periodic = data_source.ds.periodicity
        le = data_source.ds.domain_left_edge.to("code_length").d
        re = data_source.ds.domain_right_edge.to("code_length").d
        x_min, y_min, z_min = le
        x_max, y_max, z_max = re
        bounds = [x_min, x_max, y_min, y_max, z_min, z_max]
        # only need (rotated) x/y widths
        _width = (width.to("code_length").d)[:2]
        finfo = data_source.ds.field_info[item]
        ounits = finfo.output_units
        kernel_name = None
        if hasattr(data_source.ds, "kernel_name"):
            kernel_name = data_source.ds.kernel_name
        if kernel_name is None:
            kernel_name = "cubic"
        if weight is None:
            for chunk in data_source.chunks([], "io"):
                off_axis_projection_SPH(
                    chunk[ptype, ppos[0]].to("code_length").d,
                    chunk[ptype, ppos[1]].to("code_length").d,
                    chunk[ptype, ppos[2]].to("code_length").d,
                    chunk[ptype, "mass"].to("code_mass").d,
                    chunk[ptype, "density"].to("code_density").d,
                    chunk[ptype, "smoothing_length"].to("code_length").d,
                    bounds,
                    center.to("code_length").d,
                    _width,
                    periodic,
                    chunk[item].in_units(ounits),
                    buf,
                    mask,
                    normal_vector,
                    north,
                    depth=depth,
                    kernel_name=kernel_name,
                )

            # Assure that the path length unit is in the default length units
            # for the dataset by scaling the units of the smoothing length,
            # which in the above calculation is set to be code_length
            path_length_unit = Unit(
                "code_length", registry=data_source.ds.unit_registry
            )
            default_path_length_unit = data_source.ds.unit_system["length"]
            buf *= data_source.ds.quan(1, path_length_unit).in_units(
                default_path_length_unit
            )
            item_unit = data_source.ds._get_field_info(item).units
            item_unit = Unit(item_unit, registry=data_source.ds.unit_registry)
            funits = item_unit * default_path_length_unit

        else:
            # if there is a weight field, take two projections:
            # one of field*weight, the other of just weight, and divide them
            weight_buff = np.zeros((resolution[0], resolution[1]), dtype="float64")
            wounits = data_source.ds.field_info[weight].output_units

            for chunk in data_source.chunks([], "io"):
                off_axis_projection_SPH(
                    chunk[ptype, ppos[0]].to("code_length").d,
                    chunk[ptype, ppos[1]].to("code_length").d,
                    chunk[ptype, ppos[2]].to("code_length").d,
                    chunk[ptype, "mass"].to("code_mass").d,
                    chunk[ptype, "density"].to("code_density").d,
                    chunk[ptype, "smoothing_length"].to("code_length").d,
                    bounds,
                    center.to("code_length").d,
                    _width,
                    periodic,
                    chunk[item].in_units(ounits),
                    buf,
                    mask,
                    normal_vector,
                    north,
                    weight_field=chunk[weight].in_units(wounits),
                    depth=depth,
                    kernel_name=kernel_name,
                )

            for chunk in data_source.chunks([], "io"):
                off_axis_projection_SPH(
                    chunk[ptype, ppos[0]].to("code_length").d,
                    chunk[ptype, ppos[1]].to("code_length").d,
                    chunk[ptype, ppos[2]].to("code_length").d,
                    chunk[ptype, "mass"].to("code_mass").d,
                    chunk[ptype, "density"].to("code_density").d,
                    chunk[ptype, "smoothing_length"].to("code_length").d,
                    bounds,
                    center.to("code_length").d,
                    _width,
                    periodic,
                    chunk[weight].to(wounits),
                    weight_buff,
                    mask,
                    normal_vector,
                    north,
                    depth=depth,
                    kernel_name=kernel_name,
                )

            normalization_2d_utility(buf, weight_buff)
            item_unit = data_source.ds._get_field_info(item).units
            item_unit = Unit(item_unit, registry=data_source.ds.unit_registry)
            funits = item_unit

        myinfo = {
            "field": item,
            "east_vector": east_vector,
            "north_vector": north_vector,
            "normal_vector": normal_vector,
            "width": width,
            "depth": depth,
            "units": funits,
            "type": "SPH smoothed projection",
        }

        return ImageArray(
            buf, funits, registry=data_source.ds.unit_registry, info=myinfo
        )

    sc = Scene()
    data_source.ds.index
    if item is None:
        field = data_source.ds.field_list[0]
        mylog.info("Setting default field to %s", field.__repr__())

    funits = data_source.ds._get_field_info(item).units

    vol = KDTreeVolumeSource(data_source, item)
    vol.num_threads = num_threads
    if weight is None:
        vol.set_field(item)
    else:
        # This is a temporary field, which we will remove at the end.
        weightfield = ("index", "temp_weightfield")

        def _make_wf(f, w):
            def temp_weightfield(field, data):
                tr = data[f].astype("float64") * data[w]
                return tr.d

            return temp_weightfield

        data_source.ds.field_info.add_field(
            weightfield,
            sampling_type="cell",
            function=_make_wf(item, weight),
            units="",
        )
        # Now we have to tell the dataset to add it and to calculate
        # its dependencies..
        deps, _ = data_source.ds.field_info.check_derived_fields([weightfield])
        data_source.ds.field_dependencies.update(deps)
        vol.set_field(weightfield)
        vol.set_weight_field(weight)
    ptf = ProjectionTransferFunction()
    vol.set_transfer_function(ptf)
    camera = sc.add_camera(data_source)
    camera.set_width(width)
    if not is_sequence(resolution):
        resolution = [resolution] * 2
    camera.resolution = resolution
    if not is_sequence(width):
        width = data_source.ds.arr([width] * 3)
    normal = np.array(normal_vector)
    normal = normal / np.linalg.norm(normal)

    camera.position = center - width[2] * normal
    camera.focus = center

    # If north_vector is None, we set the default here.
    # This is chosen so that if normal_vector is one of the
    # cartesian coordinate axes, the projection will match
    # the corresponding on-axis projection.
    if north_vector is None:
        vecs = np.identity(3)
        t = np.cross(vecs, normal).sum(axis=1)
        ax = t.argmax()
        east_vector = np.cross(vecs[ax, :], normal).ravel()
        north = np.cross(normal, east_vector).ravel()
    else:
        north = np.array(north_vector)
        north = north / np.linalg.norm(north)
    camera.switch_orientation(normal, north)

    sc.add_source(vol)

    vol.set_sampler(camera, interpolated=False)
    assert vol.sampler is not None

    fields = [vol.field]
    if vol.weight_field is not None:
        fields.append(vol.weight_field)

    mylog.debug("Casting rays")
    index = data_source.ds.index
    lens = camera.lens

    # This implementation is optimized for octrees with plane-parallel lenses
    # and implicitely assumes that the cells are cubic.
    # NOTE: we should be able to relax the cubic assumption to a rectangular
    #       assumption (if all cells have the same aspect ratio) with some
    #       renormalization of the coordinates and the projection axes.
    #       This is NOT done in the following.
    dom_width = data_source.ds.domain_width
    cubic_domain = dom_width.max() == dom_width.min()
    if (
        isinstance(index, OctreeIndex)
        and isinstance(lens, PlaneParallelLens)
        and cubic_domain
    ):
        fields.extend(("index", k) for k in "xyz")
        fields.append(("index", "dx"))

        data_source.get_data(fields)
        # We need the width of the plot window in projected coordinates,
        # i.e. we ignore the z-component
        wmax = width[:2].max().to("code_length")
        xyz = data_source.ds.arr(
            np.zeros((len(data_source[vol.field]), 3)), "code_length"
        )

        for idim, periodic in enumerate(data_source.ds.periodicity):
            axis = data_source.ds.coordinates.axis_order[idim]
            # Recenter positions w.r.t. center of the plot window
            xyz[..., idim] = (data_source["index", axis] - center[idim]).to(
                "code_length"
            )
            if not periodic:
                continue
            # If we have periodic boundaries, we need to wrap the corresponding
            # coordinates into [-w/2, +w/2]
            w = data_source.ds.domain_width[idim].to("code_length")
            xyz[..., idim] = (xyz[..., idim] + w / 2) % w - w / 2

        # Rescale to [-0.5, +0.5]
        xyz = (xyz / wmax).to("1").d

        dx = (data_source["index", "dx"] / wmax).to("1").d

        if vol.weight_field is None:
            weight_field = np.ones_like(dx)
        else:
            weight_field = data_source[vol.weight_field]

        projected_weighted_qty = np.zeros(resolution)
        projected_weight = np.zeros(resolution)

        add_cells_to_image_offaxis(
            Xp=xyz,
            dXp=dx,
            qty=data_source[vol.field],
            weight=weight_field,
            rotation=camera.inv_mat.T,
            buffer=projected_weighted_qty,
            buffer_weight=projected_weight,
            Nx=resolution[0],
            Ny=resolution[1],
        )
        # Note: since dx was divided by wmax, we need to rescale by it
        projected_weighted_qty *= wmax.d / np.sqrt(3)
        projected_weight *= wmax.d / np.sqrt(3)

        image = ImageArray(
            data_source.ds.arr(
                np.stack([projected_weighted_qty, projected_weight], axis=-1),
                "dimensionless",
            ),
            funits,
            registry=data_source.ds.unit_registry,
            info={"imtype": "rendering"},
        )

    else:
        for grid, mask in data_source.blocks:
            data = []
            for f in fields:
                # strip units before multiplying by mask for speed
                grid_data = grid[f]
                units = grid_data.units
                data.append(
                    data_source.ds.arr(grid_data.d * mask, units, dtype="float64")
                )
            pg = PartitionedGrid(
                grid.id,
                data,
                mask.astype("uint8"),
                grid.LeftEdge,
                grid.RightEdge,
                grid.ActiveDimensions.astype("int64"),
            )
            grid.clear_data()
            vol.sampler(pg, num_threads=num_threads)

        image = vol.finalize_image(camera, vol.sampler.aimage)

        image = ImageArray(
            image, funits, registry=data_source.ds.unit_registry, info=image.info
        )

    # Remove the temporary weight field
    if weight is not None:
        data_source.ds.field_info.pop(("index", "temp_weightfield"))

    if method == "integrate":
        if weight is None:
            dl = width[2].in_units(data_source.ds.unit_system["length"])
            image *= dl
        else:
            mask = image[:, :, 1] == 0
            nmask = np.logical_not(mask)
            image[:, :, 0][nmask] /= image[:, :, 1][nmask]
            image[mask] = 0

    return image[:, :, 0]
