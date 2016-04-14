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
    Scene, \
    create_scene


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


def compare(ds, im, test_prefix, decimals=12):
    def mesh_render_image_func(filename_prefix):
        return im.write_image(filename_prefix)

    mesh_render_image_func.__name__ = "func_{}".format(test_prefix)
    test = GenericImageTest(ds, mesh_render_image_func, decimals)
    test.prefix = test_prefix
    return test

hex8 = "MOOSE_sample_data/out.e-s010"
hex8_fields = [('connect1', 'diffused'), ('connect2', 'convected')]

@requires_ds(hex8)
@requires_module("pyembree")
def test_hex8_render():
    for field in hex8_fields:
        ds = data_dir_load(hex8, kwargs={'step':-1})
        sc = create_scene(ds, field)
        im = sc.render()
        yield compare(ds, im, "render_answers_hex8_%s_%s" % field)


tet4 = "MOOSE_sample_data/high_order_elems_tet4_refine_out.e"
tet4_fields = [("connect1", "u")]

@requires_ds(tet4)
@requires_module("pyembree")
def test_tet4_render():
    for field in tet4_fields:
        ds = data_dir_load(tet4, kwargs={'step':-1})
        sc = create_scene(ds, field)
        im = sc.render()
        yield compare(ds, im, "render_answers_tet4_%s_%s" % field)


hex20 = "MOOSE_sample_data/mps_out.e"
hex20_fields = [('connect2', 'temp')]

@requires_ds(hex20)
@requires_module("pyembree")
def test_hex20_render():
    for field in hex20_fields:
        ds = data_dir_load(hex20, kwargs={'step':-1})
        sc = create_scene(ds, field)
        im = sc.render()
        yield compare(ds, im, "render_answers_hex20_%s_%s" % field)


wedge6 = "MOOSE_sample_data/wedge_out.e"
wedge6_fields = [('connect1', 'diffused')]

@requires_ds(wedge6)
@requires_module("pyembree")
def test_wedge6_render():
    for field in wedge6_fields:
        ds = data_dir_load(wedge6, kwargs={'step':-1})
        sc = create_scene(ds, field)
        im = sc.render()
        yield compare(ds, im, "render_answers_wedge6_%s_%s" % field)


@requires_ds(hex8)
@requires_module("pyembree")
def test_perspective_mesh_render():
    ds = data_dir_load(hex8)
    sc = create_scene(ds, ("connect2", "diffused"))

    cam = sc.add_camera(ds, lens_type='perspective')
    cam.focus = ds.arr([0.0, 0.0, 0.0], 'code_length')
    cam_pos = ds.arr([-4.5, 4.5, -4.5], 'code_length')
    north_vector = ds.arr([0.0, -1.0, -1.0], 'dimensionless')
    cam.set_position(cam_pos, north_vector)
    cam.resolution = (800, 800)
    im = sc.render()
    yield compare(ds, im, "perspective_mesh_render")


@requires_ds(hex8)
@requires_module("pyembree")
def test_composite_mesh_render():
    ds = data_dir_load(hex8)
    sc = Scene()
    cam = sc.add_camera(ds)
    cam.focus = ds.arr([0.0, 0.0, 0.0], 'code_length')
    cam.set_position(ds.arr([-3.0, 3.0, -3.0], 'code_length'),
                     ds.arr([0.0, -1.0, 0.0], 'dimensionless'))
    cam.set_width = ds.arr([8.0, 8.0, 8.0], 'code_length')
    cam.resolution = (800, 800)

    ms1 = MeshSource(ds, ('connect1', 'diffused'))
    ms2 = MeshSource(ds, ('connect2', 'diffused'))

    sc.add_source(ms1)
    sc.add_source(ms2)

    im = sc.render()
    yield compare(ds, im, "composite_mesh_render")
