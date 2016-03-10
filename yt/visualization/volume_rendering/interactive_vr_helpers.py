"""
Helper routines for OpenGL VR
"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
from yt.utilities.on_demand_imports import NotAModule
from yt.data_objects.static_output import Dataset
import numpy as np

def _render_opengl(data_source, field=None, window_size=None, cam_position=None,
                   cam_focus=None):
    from .interactive_vr import SceneGraph, BlockCollection, TrackballCamera
    from .interactive_loop import RenderingContext

    if isinstance(data_source, Dataset):
        dobj = data_source.all_data()
    else:
        dobj = data_source
    if field is None:
        field = ("gas", "density")
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
