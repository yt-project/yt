"""
Test Surface Mesh Rendering

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import yt
import numpy as np
from yt.testing import requires_module
from yt.visualization.volume_rendering.render_source import MeshSource
from yt.visualization.volume_rendering.camera import Camera
from yt.frontends.stream.sample_data.unstructured_mesh import \
    _connectivity, \
    _coordinates

@requires_module("pyembree")
def test_surface_mesh_render():
    data ={}
    data[('connect1', 'test')] = np.ones_like(_connectivity)
    ds = yt.load_unstructured_mesh(_connectivity, _coordinates, node_data=data)
    ms = MeshSource(ds, ('connect1', 'test'))
    cam = Camera(ds)
    im = ms.render(cam)
    return im


if __name__ == "__main__":
    im = test_surface_mesh_render()
