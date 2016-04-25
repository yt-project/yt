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
from yt.config import \
    ytcfg

def compare(ds, im, test_prefix, decimals=12):
    def mesh_render_image_func(filename_prefix):
        return im.write_image(filename_prefix)

    mesh_render_image_func.__name__ = "func_{}".format(test_prefix)
    test = GenericImageTest(ds, mesh_render_image_func, decimals)
    test.prefix = test_prefix
    return test

def surface_mesh_render():
    images = []

    ds = fake_tetrahedral_ds()
    for field in ds.field_list:
        sc = Scene()
        sc.add_source(MeshSource(ds, field))
        sc.add_camera()
        im = sc.render()
        images.append(im)

    ds = fake_hexahedral_ds()
    for field in ds.field_list:
        sc = Scene()
        sc.add_source(MeshSource(ds, field))
        sc.add_camera()
        im = sc.render()
        images.append(im)

    return images
    
@requires_module("pyembree")
def test_surface_mesh_render_pyembree():
    ytcfg["yt", "ray_tracing_engine"] = "embree"
    surface_mesh_render()

def test_surface_mesh_render():
    ytcfg["yt", "ray_tracing_engine"] = "yt"
    surface_mesh_render()


hex8 = "MOOSE_sample_data/out.e-s010"
hex8_fields = [('connect1', 'diffused'), ('connect2', 'convected')]

def hex8_render(engine, field):
    ytcfg["yt", "ray_tracing_engine"] = engine
    ds = data_dir_load(hex8, kwargs={'step':-1})
    sc = create_scene(ds, field)
    im = sc.render()
    return compare(ds, im, "%s_render_answers_hex8_%s_%s" %
                   (engine, field[0], field[1]))

@requires_ds(hex8)
@requires_module("pyembree")
def test_hex8_render_pyembree():
    for field in hex8_fields:
        yield hex8_render("embree", field)

@requires_ds(hex8)
def test_hex8_render():
    for field in hex8_fields:
        yield hex8_render("yt", field)


tet4 = "MOOSE_sample_data/high_order_elems_tet4_refine_out.e"
tet4_fields = [("connect1", "u")]

def tet4_render(engine, field):
    ytcfg["yt", "ray_tracing_engine"] = engine
    ds = data_dir_load(tet4, kwargs={'step':-1})
    sc = create_scene(ds, field)
    im = sc.render()
    return compare(ds, im, "%s_render_answers_tet4_%s_%s" %
                   (engine, field[0], field[1]))

@requires_ds(tet4)
@requires_module("pyembree")
def test_tet4_render_pyembree():
    for field in tet4_fields:
        yield tet4_render("embree", field)

@requires_ds(tet4)
def test_tet4_render():
    for field in tet4_fields:
        yield tet4_render("yt", field)


hex20 = "MOOSE_sample_data/mps_out.e"
hex20_fields = [('connect2', 'temp')]

def hex20_render(engine, field):
    ytcfg["yt", "ray_tracing_engine"] = engine
    ds = data_dir_load(hex20, kwargs={'step':-1})
    sc = create_scene(ds, field)
    im = sc.render()
    return compare(ds, im, "%s_render_answers_hex20_%s_%s" %
                   (engine, field[0], field[1]))

@requires_ds(hex20)
@requires_module("pyembree")
def test_hex20_render_pyembree():
    for field in hex20_fields:
        yield hex20_render("embree", field)

@requires_ds(hex20)
def test_hex20_render():
    for field in hex20_fields:
        yield hex20_render("yt", field)


wedge6 = "MOOSE_sample_data/wedge_out.e"
wedge6_fields = [('connect1', 'diffused')]

def wedge6_render(engine, field):
    ytcfg["yt", "ray_tracing_engine"] = engine
    ds = data_dir_load(wedge6, kwargs={'step':-1})
    sc = create_scene(ds, field)
    im = sc.render()
    return compare(ds, im, "%s_render_answers_wedge6_%s_%s" %
                   (engine, field[0], field[1]))
    
@requires_ds(wedge6)
@requires_module("pyembree")
def test_wedge6_render_pyembree():
    for field in wedge6_fields:
        yield wedge6_render("embree", field)

@requires_ds(wedge6)
def test_wedge6_render():
    for field in wedge6_fields:
        yield wedge6_render("yt", field)

def perspective_mesh_render(engine):
    ytcfg["yt", "ray_tracing_engine"] = engine
    ds = data_dir_load(hex8)
    sc = create_scene(ds, ("connect2", "diffused"))
    cam = sc.add_camera(ds, lens_type='perspective')
    cam.focus = ds.arr([0.0, 0.0, 0.0], 'code_length')
    cam_pos = ds.arr([-4.5, 4.5, -4.5], 'code_length')
    north_vector = ds.arr([0.0, -1.0, -1.0], 'dimensionless')
    cam.set_position(cam_pos, north_vector)
    cam.resolution = (800, 800)
    im = sc.render()
    return compare(ds, im, "%s_perspective_mesh_render" % engine)
    
@requires_ds(hex8)
@requires_module("pyembree")
def test_perspective_mesh_render_pyembree():
    yield perspective_mesh_render("embree")

@requires_ds(hex8)
def test_perspective_mesh_render():
    yield perspective_mesh_render("yt")

def composite_mesh_render(engine):
    ytcfg["yt", "ray_tracing_engine"] = engine
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
    return compare(ds, im, "%s_composite_mesh_render" % engine)
    
@requires_ds(hex8)
@requires_module("pyembree")
def test_composite_mesh_render_pyembree():
    yield composite_mesh_render("embree")

@requires_ds(hex8)
def test_composite_mesh_render():
    yield composite_mesh_render("yt")
