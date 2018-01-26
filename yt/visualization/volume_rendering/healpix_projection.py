"""
HEALPix projection

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2018, yt Development Team.
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


def healpix_projection(data_source, center, radius, nside, item,
                       weight=None, volume=None,
                       no_ghost=False, interpolated=False,
                       normal_vector=None, north_vector=None,
                       num_threads=1, method='integrate'):
    if method not in ['integrate','sum']:
        raise NotImplementedError("Only 'integrate' or 'sum' methods are valid for off-axis-projections")

    if interpolated is True:
        raise NotImplementedError("Only interpolated=False methods are currently implemented for off-axis-projections")


    data_source = data_source_or_all(data_source)
    # Sanitize units
    if not hasattr(center, "units"):
        center = data_source.ds.arr(center, 'code_length')
    if not hasattr(radius, "units"):
        radius = data_source.ds.arr(radius, 'code_length')
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
            return temp_weightfield
        data_source.ds.field_info.add_field(weightfield, sampling_type="cell",
            function=_make_wf(item, weight))
        # Now we have to tell the dataset to add it and to calculate
        # its dependencies..
        deps, _ = data_source.ds.field_info.check_derived_fields([weightfield])
        data_source.ds.field_dependencies.update(deps)
        vol.set_field(weightfield)
        vol.set_weight_field(weight)
    ptf = ProjectionTransferFunction()
    vol.set_transfer_function(ptf)
    camera = sc.add_camera(data_source, lens_type='healpix')
    camera.set_width(radius)
    camera.lens.nside = nside

    if normal_vector is None:
        normal_vector = (1, 0, 0)
    normal = np.array(normal_vector)
    normal = normal / np.linalg.norm(normal)

    camera.focus = 0  # This is a hack. Need to be refined.
    camera.position = center
    camera.focus = center + radius * normal

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

    if weight is not None:
        data_source.ds.field_info.pop(("index", "temp_weightfield"))

    if method == "integrate":
        if weight is None:
            dl = radius.in_units(data_source.ds.unit_system["length"])
            image *= dl
        else:
            mask = image[:,0,1] == 0
            image[:,0,0] /= image[:,0,1]
            image[mask] = 0

    return image[:,0,0]
