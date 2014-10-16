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
from yt.funcs import mylog


def off_axis_projection(data_source, center, normal_vector,
                        width, resolution, item,
                        weight=None, volume=None,
                        no_ghost=False, interpolated=False,
                        north_vector=None):
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
    cam = Camera(data_source)
    sc.set_default_camera(cam)
    sc.add_source(vol)
    im = sc.render()
    return im
