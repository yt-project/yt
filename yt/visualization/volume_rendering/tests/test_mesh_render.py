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
from yt.testing import fake_tetrahedral_ds
from yt.testing import fake_hexahedral_ds
from yt.testing import requires_module
from yt.visualization.volume_rendering.render_source import MeshSource
from yt.visualization.volume_rendering.camera import Camera


@requires_module("pyembree")
def test_surface_mesh_render():
    ds = fake_tetrahedral_ds()
    ms = MeshSource(ds, ('connect1', 'test'))
    cam = Camera(ds)
    im1 = ms.render(cam)

    ds = fake_hexahedral_ds()
    ms = MeshSource(ds, ('connect1', 'test'))
    cam = Camera(ds)
    im2 = ms.render(cam)
    
    return im1, im2

if __name__ == "__main__":
    im1, im2 = test_surface_mesh_render()
