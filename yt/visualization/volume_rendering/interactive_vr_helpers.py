# encoding: utf-8
"""
Helper routines for Interactive Data Visualization
"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

# This part of the experimental Interactive Data Visualization

import numpy as np
from yt.funcs import mylog
from yt.data_objects.static_output import Dataset
from yt.utilities.exceptions import YTSceneFieldNotFound
from yt.utilities.on_demand_imports import NotAModule

def _render_opengl(data_source, field=None, window_size=None, cam_position=None,
                   cam_focus=None):
    '''High level wrapper for Interactive Data Visualization
    
    Parameters
    ----------
    data_source : :class:`yt.data_objects.data_containers.AMR3DData`,
                  :class:`yt.data_objects.static_output.Dataset`
        This is the source to be rendered, which can be any arbitrary yt
        3D object
    field : string, tuple, optional
        The field to be rendered. If unspecified, this will use the
        default_field for your dataset's frontend--usually ('gas', 'density').
    window_size : 2 element tuple of ints
        The width and the height of the Interactive Data Visualization window.
        For performance reasons it is recommended to use values that are natural
        powers of 2.
    cam_position : 3 element YTArray, optional
        The camera position in physical coordinates. If unspecified,
        data_source's domain right edge will be used.
    cam_focus: 3 element YTArray, optional
        The focus defines the point the camera is pointed at. If unspecified,
        data_source's domain center will be used.

    Examples
    --------

    >>> import yt
    >>> ds = yt.load("Enzo_64/DD0046/DD0046")
    >>> yt.interactive_render(ds)

    '''

    from .interactive_vr import SceneGraph, BlockCollection, TrackballCamera
    from .interactive_loop import RenderingContext

    if isinstance(data_source, Dataset):
        dobj = data_source.all_data()
    else:
        dobj = data_source
    if field is None:
        field = dobj.ds.default_field
        if field not in dobj.ds.derived_field_list:
            raise YTSceneFieldNotFound("""Could not find field '%s' in %s.
                  Please specify a field in create_scene()""" % (field, data_source.ds))
        mylog.info('Setting default field to %s' % field.__repr__())
    if window_size is None:
        window_size = (1024, 1024)
    if cam_position is None:
        cam_position = dobj.ds.domain_right_edge
    if cam_focus is None:
        cam_focus = dobj.ds.domain_center

    aspect_ratio = window_size[1] / window_size[0]
    far_plane = np.linalg.norm(cam_focus - cam_position) * 2.0
    near_plane = 0.01 * far_plane

    rc = RenderingContext(*window_size)
    scene = SceneGraph()
    collection = BlockCollection()
    collection.add_data(dobj, field)
    scene.add_collection(collection)

    c = TrackballCamera(position=cam_position, focus=cam_focus, near_plane=near_plane,
                        far_plane=far_plane, aspect_ratio=aspect_ratio)
    rc.start_loop(scene, c)


try:
    import cyglfw3 as glfw  # NOQA
    import OpenGL.GL as GL  # NOQA
    interactive_render = _render_opengl
except ImportError:
    interactive_render = NotAModule("opengl/cyglfw3")
