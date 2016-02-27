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

from yt.testing import \
    fake_tetrahedral_ds, \
    fake_hexahedral_ds, \
    requires_module
from yt.utilities.answer_testing.framework import \
    requires_ds, \
    data_dir_load, \
    GenericImageTest
from yt.visualization.volume_rendering.api import \
    MeshSource, \
    Camera, \
    create_scene


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


hex8 = "MOOSE_sample_data/out.e-s010"
hex8_fields = [('connect1', 'diffused'), ('connect2', 'convected')]

@requires_ds(hex8)
@requires_module("pyembree")
def test_hex8_render():
    for field in hex8_fields:
        ds = data_dir_load(hex8, kwargs={'step':-1})
        sc = create_scene(ds, field)
        im = sc.render()

        def mesh_render_image_func(filename_prefix):
            return im.write_image(filename_prefix)

        test = GenericImageTest(ds, mesh_render_image_func, 12)
        test.prefix = "render_answers_hex8_%s_%s" % field
        yield test


tet4 = "MOOSE_sample_data/high_order_elems_tet4_refine_out.e"
tet4_fields = [("connect1", "u")]

@requires_ds(tet4)
@requires_module("pyembree")
def test_tet4_render():
    for field in tet4_fields:
        ds = data_dir_load(tet4, kwargs={'step':-1})
        sc = create_scene(ds, field)
        im = sc.render()

        def mesh_render_image_func(filename_prefix):
            return im.write_image(filename_prefix)

        test = GenericImageTest(ds, mesh_render_image_func, 12)
        test.prefix = "render_answers_tet4_%s_%s" % field
        yield test


hex20 = "MOOSE_sample_data/mps_out.e"
hex20_fields = [('connect2', 'temp')]

@requires_ds(hex20)
@requires_module("pyembree")
def test_hex20_render():
    for field in hex20_fields:
        ds = data_dir_load(hex20, kwargs={'step':-1})
        sc = create_scene(ds, field)
        im = sc.render()

        def mesh_render_image_func(filename_prefix):
            return im.write_image(filename_prefix)

        test = GenericImageTest(ds, mesh_render_image_func, 12)
        test.prefix = "render_answers_hex20_%s_%s" % field
        yield test
