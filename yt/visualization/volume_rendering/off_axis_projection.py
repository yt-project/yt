"""
Volume rendering

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------


from .scene import Scene
from .render_source import VolumeSource
from .transfer_functions import ProjectionTransferFunction
from .utils import data_source_or_all
from yt.funcs import mylog, iterable
from yt.utilities.lib.partitioned_grid import \
    PartitionedGrid
from yt.data_objects.api import ImageArray
import numpy as np


def off_axis_projection(data_source, center, normal_vector,
                        width, resolution, item,
                        weight=None, volume=None,
                        no_ghost=False, interpolated=False,
                        north_vector=None, num_threads=1, method='integrate'):
    r"""Project through a dataset, off-axis, and return the image plane.

    This function will accept the necessary items to integrate through a volume
    at an arbitrary angle and return the integrated field of view to the user.
    Note that if a weight is supplied, it will multiply the pre-interpolated
    values together, then create cell-centered values, then interpolate within
    the cell to conduct the integration.

    Parameters
    ----------
    data_source : `~yt.data_objects.static_output.Dataset` or `~yt.data_objects.data_containers.YTSelectionDataContainer`
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

    >>> image = off_axis_projection(ds, [0.5, 0.5, 0.5], [0.2,0.3,0.4],
    ...                             0.2, N, "temperature", "density")
    >>> write_image(np.log10(image), "offaxis.png")

    """

    if method not in ['integrate','sum']:
        raise NotImplementedError("Only 'integrate' or 'sum' methods are valid for off-axis-projections")

    if interpolated is True:
        raise NotImplementedError("Only interpolated=False methods are currently implemented for off-axis-projections")


    data_source = data_source_or_all(data_source)
    # Sanitize units
    if not hasattr(center, "units"):
        center = data_source.ds.arr(center, 'code_length')
    if not hasattr(width, "units"):
        width = data_source.ds.arr(width, 'code_length')
    sc = Scene()
    data_source.ds.index
    if item is None:
        field = data_source.ds.field_list[0]
        mylog.info('Setting default field to %s' % field.__repr__())

    funits = data_source.ds._get_field_info(item).units

    vol = VolumeSource(data_source, item)
    if weight is None:
        vol.set_field(item)
    else:
        # This is a temporary field, which we will remove at the end.
        weightfield = ("index", "temp_weightfield")
        def _make_wf(f, w):
            def temp_weightfield(a, b):
                tr = b[f].astype("float64") * b[w]
                return b.apply_units(tr, a.units)
                return tr
            return temp_weightfield
        data_source.ds.field_info.add_field(weightfield,
            function=_make_wf(item, weight))
        # Now we have to tell the dataset to add it and to calculate
        # its dependencies..
        deps, _ = data_source.ds.field_info.check_derived_fields([weightfield])
        data_source.ds.field_dependencies.update(deps)
        fields = [weightfield, weight]
        vol.set_field(weightfield)
        vol.set_weight_field(weight)
    ptf = ProjectionTransferFunction()
    vol.set_transfer_function(ptf)
    camera = sc.add_camera(data_source)
    camera.set_width(width)
    if not iterable(resolution):
        resolution = [resolution]*2
    camera.resolution = resolution
    if not iterable(width):
        width = data_source.ds.arr([width]*3)
    normal = np.array(normal_vector)
    normal = normal / np.linalg.norm(normal)

    camera.position = center - width[2]*normal
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
    assert (vol.sampler is not None)

    fields = [vol.field]
    if vol.weight_field is not None:
        fields.append(vol.weight_field)

    mylog.debug("Casting rays")

    for i, (grid, mask) in enumerate(data_source.blocks):
        data = [(grid[f] * mask).astype("float64") for f in fields]
        pg = PartitionedGrid(
            grid.id, data,
            mask.astype('uint8'),
            grid.LeftEdge, grid.RightEdge, grid.ActiveDimensions.astype("int64"))
        grid.clear_data()
        vol.sampler(pg, num_threads = num_threads)

    image = vol.finalize_image(camera, vol.sampler.aimage)
    image = ImageArray(image, funits, registry=data_source.ds.unit_registry, info=image.info)

    if method == "integrate":
        if weight is None:
            dl = width[2].in_units(data_source.ds.unit_system["length"])
            image *= dl
        else:
            mask = image[:,:,1] == 0
            image[:,:,0] /= image[:,:,1]
            image[mask] = 0

    return image[:,:,0]
