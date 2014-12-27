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
from .camera import Camera
from .render_source import VolumeSource
from .transfer_functions import ProjectionTransferFunction
from .utils import data_source_or_all
from yt.funcs import mylog, iterable
from yt.utilities.lib.grid_traversal import \
        PartitionedGrid
import numpy as np


def off_axis_projection(data_source, center, normal_vector,
                        width, resolution, item,
                        weight=None, volume=None,
                        no_ghost=False, interpolated=False,
                        north_vector=None, num_threads=1):
    data_source = data_source_or_all(data_source)
    sc = Scene()
    data_source.ds.index
    if item is None:
        field = data_source.pf.field_list[0]
        mylog.info('Setting default field to %s' % field.__repr__())

    vol = VolumeSource(data_source, item)
    ptf = ProjectionTransferFunction()
    vol.set_transfer_function(ptf)
    if weight is None:
        vol.set_fields([item])
    else:
        vol.set_fields([item, weight])
    camera = Camera(data_source)
    sc.set_default_camera(camera)
    sc.add_source(vol)

    camera.lens.set_camera(camera)
    vol.set_sampler(camera)
    assert (vol.sampler is not None)

    mylog.debug("Casting rays")
    total_cells = 0
    double_check = False
    if double_check:
        for brick in vol.volume.bricks:
            for data in brick.my_data:
                if np.any(np.isnan(data)):
                    raise RuntimeError

    ds = data_source.ds
    north_vector = camera.unit_vectors[0]
    east_vector = camera.unit_vectors[1]
    normal_vector = camera.unit_vectors[2]
    fields = vol.field 
    if not iterable(width):
        width = data_source.ds.arr([width]*3) 

    mi = ds.domain_right_edge.copy()
    ma = ds.domain_left_edge.copy()
    for off1 in [-1, 1]:
        for off2 in [-1, 1]:
            for off3 in [-1, 1]:
                this_point = (center + width[0]/2. * off1 * north_vector
                                     + width[1]/2. * off2 * east_vector
                                     + width[2]/2. * off3 * normal_vector)
                np.minimum(mi, this_point, mi)
                np.maximum(ma, this_point, ma)
    # Now we have a bounding box.
    data_source = ds.region(center, mi, ma)

    for i, (grid, mask) in enumerate(data_source.blocks):
        data = [(grid[field] * mask).astype("float64") for field in fields]
        pg = PartitionedGrid(
            grid.id, data,
            mask.astype('uint8'),
            grid.LeftEdge, grid.RightEdge, grid.ActiveDimensions.astype("int64"))
        grid.clear_data()
        vol.sampler(pg, num_threads = num_threads)

    image = vol.finalize_image(camera, vol.sampler.aimage)

    if weight is None:
        dl = width[2]
        image *= dl
    else:
        image[:,:,0] /= image[:,:,1]

    return image[:,:,0], sc
