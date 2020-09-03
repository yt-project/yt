from nose.plugins.attrib import attr

from yt.config import ytcfg
from yt.testing import (
    ANSWER_TEST_TAG,
    fake_hexahedral_ds,
    fake_tetrahedral_ds,
    requires_module,
)
from yt.utilities.answer_testing.framework import (
    GenericImageTest,
    data_dir_load,
    requires_ds,
)
from yt.visualization.volume_rendering.api import MeshSource, Scene, create_scene


def compare(ds, im, test_prefix, test_name=None, decimals=12):
    def mesh_render_image_func(filename_prefix):
        return im.write_image(filename_prefix)

    mesh_render_image_func.__name__ = f"func_{test_prefix}"
    test = GenericImageTest(ds, mesh_render_image_func, decimals)
    test.prefix = test_prefix
    test.answer_name = test_name
    return test


def surface_mesh_render():
    images = []

    ds = fake_tetrahedral_ds()
    for field in ds.field_list:
        if field[0] == "all":
            continue
        sc = Scene()
        sc.add_source(MeshSource(ds, field))
        sc.add_camera()
        im = sc.render()
        images.append(im)

    ds = fake_hexahedral_ds()
    for field in ds.field_list:
        if field[0] == "all":
            continue
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


@attr(ANSWER_TEST_TAG)
def test_fake_hexahedral_ds_render():
    ds = fake_hexahedral_ds()
    field_list = [("connect1", "elem"), ("connect1", "test")]
    for field in field_list:
        sc = create_scene(ds, field)
        im = sc.render()
        test_prefix = f"yt_render_fake_hexahedral_{field[0]}_{field[1]}"
        yield compare(
            ds, im, test_prefix=test_prefix, test_name="fake_hexahedral_ds_render"
        )


hex8 = "MOOSE_sample_data/out.e-s010"
hex8_fields = [("connect1", "diffused"), ("connect2", "convected")]


def hex8_render(engine, field):
    ytcfg["yt", "ray_tracing_engine"] = engine
    ds = data_dir_load(hex8, kwargs={"step": -1})
    sc = create_scene(ds, field)
    im = sc.render()
    return compare(ds, im, f"{engine}_render_answers_hex8_{field[0]}_{field[1]}")


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
    ds = data_dir_load(tet4, kwargs={"step": -1})
    sc = create_scene(ds, field)
    im = sc.render()
    return compare(ds, im, f"{engine}_render_answers_tet4_{field[0]}_{field[1]}")


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
hex20_fields = [("connect2", "temp")]


def hex20_render(engine, field):
    ytcfg["yt", "ray_tracing_engine"] = engine
    ds = data_dir_load(hex20, kwargs={"step": -1})
    sc = create_scene(ds, field)
    im = sc.render()
    return compare(ds, im, f"{engine}_render_answers_hex20_{field[0]}_{field[1]}")


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
wedge6_fields = [("connect1", "diffused")]


def wedge6_render(engine, field):
    ytcfg["yt", "ray_tracing_engine"] = engine
    ds = data_dir_load(wedge6, kwargs={"step": -1})
    sc = create_scene(ds, field)
    im = sc.render()
    return compare(ds, im, f"{engine}_render_answers_wedge6_{field[0]}_{field[1]}")


@requires_ds(wedge6)
@requires_module("pyembree")
def test_wedge6_render_pyembree():
    for field in wedge6_fields:
        yield wedge6_render("embree", field)


@requires_ds(wedge6)
def test_wedge6_render():
    for field in wedge6_fields:
        yield wedge6_render("yt", field)


tet10 = "SecondOrderTets/tet10_unstructured_out.e"
tet10_fields = [("connect1", "uz")]


def tet10_render(engine, field):
    ytcfg["yt", "ray_tracing_engine"] = engine
    ds = data_dir_load(tet10, kwargs={"step": -1})
    sc = create_scene(ds, field)
    ms = sc.get_source(0)
    ms.color_bounds = (-0.01, 0.2)
    im = sc.render()
    return compare(ds, im, f"{engine}_render_answers_tet10_{field[0]}_{field[1]}")


@requires_ds(tet10)
@requires_module("pyembree")
def test_tet10_render_pyembree():
    for field in tet10_fields:
        yield tet10_render("embree", field)


@requires_ds(tet10)
def test_tet10_render():
    for field in tet10_fields:
        yield tet10_render("yt", field)


def perspective_mesh_render(engine):
    ytcfg["yt", "ray_tracing_engine"] = engine
    ds = data_dir_load(hex8)
    sc = create_scene(ds, ("connect2", "diffused"))
    cam = sc.add_camera(ds, lens_type="perspective")
    cam.focus = ds.arr([0.0, 0.0, 0.0], "code_length")
    cam_pos = ds.arr([-4.5, 4.5, -4.5], "code_length")
    north_vector = ds.arr([0.0, -1.0, -1.0], "dimensionless")
    cam.set_position(cam_pos, north_vector)
    cam.resolution = (800, 800)
    im = sc.render()
    return compare(ds, im, f"{engine}_perspective_mesh_render")


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
    cam.focus = ds.arr([0.0, 0.0, 0.0], "code_length")
    cam.set_position(
        ds.arr([-3.0, 3.0, -3.0], "code_length"),
        ds.arr([0.0, -1.0, 0.0], "dimensionless"),
    )
    cam.set_width = ds.arr([8.0, 8.0, 8.0], "code_length")
    cam.resolution = (800, 800)
    ms1 = MeshSource(ds, ("connect1", "diffused"))
    ms2 = MeshSource(ds, ("connect2", "diffused"))
    sc.add_source(ms1)
    sc.add_source(ms2)
    im = sc.render()
    return compare(ds, im, f"{engine}_composite_mesh_render")


@requires_ds(hex8)
@requires_module("pyembree")
def test_composite_mesh_render_pyembree():
    yield composite_mesh_render("embree")


@requires_ds(hex8)
def test_composite_mesh_render():
    yield composite_mesh_render("yt")
