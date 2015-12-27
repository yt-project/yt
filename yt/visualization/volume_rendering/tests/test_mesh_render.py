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

from yt.testing import fake_tetrahedral_ds
from yt.testing import fake_hexahedral_ds
from yt.testing import requires_module
from yt.visualization.volume_rendering.render_source import MeshSource
from yt.visualization.volume_rendering.camera import Camera


@requires_module("pyembree")
def test_surface_mesh_render():

    images = []

    ds = fake_tetrahedral_ds()
    for field in ds.field_list:
        ms = MeshSource(ds, field)
        cam = Camera(ds)
        im = ms.render(cam)
        images.append(im)

    ds = fake_hexahedral_ds()
    for field in ds.field_list:
        ms = MeshSource(ds, field)
        cam = Camera(ds)
        im = ms.render(cam)
        images.append(im)

    return images
