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
from yt.visualization.volume_rendering.scene import Scene


@requires_module("pyembree")
def test_surface_mesh_render():

    images = []

    ds = fake_tetrahedral_ds()
    sc = Scene()
    for field in ds.field_list:
        sc.add_source(MeshSource(ds, field))
        sc.add_camera()
        im = sc.render()
        images.append(im)

    ds = fake_hexahedral_ds()
    for field in ds.field_list:
        sc.add_source(MeshSource(ds, field))
        sc.add_camera()
        im = sc.render()
        images.append(im)

    return images
