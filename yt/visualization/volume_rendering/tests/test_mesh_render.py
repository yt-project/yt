import matplotlib.pyplot as plt
import pytest

from yt.config import ytcfg
from yt.testing import (
    fake_hexahedral_ds,
    fake_tetrahedral_ds,
    requires_file,
    requires_module,
)
from yt.utilities.answer_testing.framework import data_dir_load, data_dir_load_v2
from yt.visualization.volume_rendering.api import MeshSource, Scene, create_scene
from yt.visualization.volume_rendering.render_source import set_raytracing_engine


@pytest.fixture
@requires_module("pyembree")
def with_pyembree_ray_tracing_engine():
    old = ytcfg["yt", "ray_tracing_engine"]

    # the @requires_module decorator only guards against pyembree not being
    # available for import, but it might not be installed properly regardless
    # so we need to be extra careful not to run tests with the default engine
    try:
        set_raytracing_engine("embree")
    except UserWarning as exc:
        pytest.skip(str(exc))
    else:
        if ytcfg["yt", "ray_tracing_engine"] != "embree":
            pytest.skip("Error while setting embree raytracing engine")

    yield

    set_raytracing_engine(old)


@pytest.fixture
def with_default_ray_tracing_engine():
    old = ytcfg["yt", "ray_tracing_engine"]
    set_raytracing_engine("yt")
    yield
    set_raytracing_engine(old)


@pytest.mark.usefixtures("with_default_ray_tracing_engine")
class TestMeshRenderDefaultEngine:
    @classmethod
    def setup_class(cls):
        cls.images = {}

        cls.ds_t = fake_tetrahedral_ds()
        for field in cls.ds_t.field_list:
            if field[0] == "all":
                continue
            sc = Scene()
            sc.add_source(MeshSource(cls.ds_t, field))
            sc.add_camera()
            im = sc.render()
            cls.images["tetrahedral_" + "_".join(field)] = im

        cls.ds_h = fake_hexahedral_ds()
        for field in cls.ds_t.field_list:
            if field[0] == "all":
                continue
            sc = Scene()
            sc.add_source(MeshSource(cls.ds_h, field))
            sc.add_camera()
            im = sc.render()
            cls.images["hexahedral_" + "_".join(field)] = im

    @pytest.mark.parametrize("kind", ["tetrahedral", "hexahedral"])
    @pytest.mark.parametrize("fname", ["elem", "test"])
    @pytest.mark.mpl_image_compare(remove_text=True)
    def test_mesh_render_default_engine(self, kind, fname):
        fig, ax = plt.subplots()
        ax.imshow(self.images[f"{kind}_connect1_{fname}"])
        return fig


@pytest.mark.usefixtures("with_pyembree_ray_tracing_engine")
class TestMeshRenderPyembreeEngine:
    @classmethod
    def setup_class(cls):
        cls.images = {}

        cls.ds_t = fake_tetrahedral_ds()
        for field in cls.ds_t.field_list:
            if field[0] == "all":
                continue
            sc = Scene()
            sc.add_source(MeshSource(cls.ds_t, field))
            sc.add_camera()
            im = sc.render()
            cls.images["tetrahedral_" + "_".join(field)] = im

        cls.ds_h = fake_hexahedral_ds()
        for field in cls.ds_t.field_list:
            if field[0] == "all":
                continue
            sc = Scene()
            sc.add_source(MeshSource(cls.ds_h, field))
            sc.add_camera()
            im = sc.render()
            cls.images["hexahedral_" + "_".join(field)] = im

    @pytest.mark.parametrize("kind", ["tetrahedral", "hexahedral"])
    @pytest.mark.parametrize("fname", ["elem", "test"])
    @pytest.mark.mpl_image_compare(remove_text=True)
    def test_mesh_render_pyembree_engine(self, kind, fname):
        fig, ax = plt.subplots()
        ax.imshow(self.images[f"{kind}_connect1_{fname}"])
        return fig


hex8 = "MOOSE_sample_data/out.e-s010"


@pytest.mark.usefixtures("with_default_ray_tracing_engine")
class TestHex8DefaultEngine:
    @requires_file(hex8)
    @classmethod
    def setup_class(cls):
        cls.ds = data_dir_load_v2(hex8, step=-1)
        cls.images = {}
        for field in [("connect1", "diffused"), ("connect2", "convected")]:
            sc = create_scene(cls.ds, field)
            im = sc.render()
            cls.images["_".join(field)] = im

    @pytest.mark.mpl_image_compare(remove_text=True)
    def test_mesh_render_default_engine_hex8_connect1_diffused(self):
        fig, ax = plt.subplots()
        ax.imshow(self.images["connect1_diffused"])
        return fig

    @pytest.mark.mpl_image_compare(remove_text=True)
    def test_mesh_render_default_engine_hex8_connect2_convected(self):
        fig, ax = plt.subplots()
        ax.imshow(self.images["connect2_convected"])
        return fig


@pytest.mark.usefixtures("with_pyembree_ray_tracing_engine")
class TestHex8PyembreeEngine:
    @requires_file(hex8)
    @classmethod
    def setup_class(cls):
        cls.ds = data_dir_load_v2(hex8, step=-1)
        cls.images = {}
        for field in [("connect1", "diffused"), ("connect2", "convected")]:
            sc = create_scene(cls.ds, field)
            im = sc.render()
            cls.images["_".join(field)] = im

    @pytest.mark.mpl_image_compare(remove_text=True)
    def test_mesh_render_pyembree_engine_hex8_connect1_diffused(self):
        fig, ax = plt.subplots()
        ax.imshow(self.images["connect1_diffused"])
        return fig

    @pytest.mark.mpl_image_compare(remove_text=True)
    def test_mesh_render_pyembree_engine_hex8_connect2_convected(self):
        fig, ax = plt.subplots()
        ax.imshow(self.images["connect2_convected"])
        return fig


tet4 = "MOOSE_sample_data/high_order_elems_tet4_refine_out.e"


@pytest.mark.usefixtures("with_default_ray_tracing_engine")
class TestTet4DefaultEngine:
    @requires_file(tet4)
    @classmethod
    def setup_class(cls):
        cls.ds = data_dir_load_v2(tet4, step=-1)
        cls.images = {}
        for field in [("connect1", "u")]:
            sc = create_scene(cls.ds, field)
            im = sc.render()
            cls.images["_".join(field)] = im

    @pytest.mark.mpl_image_compare(remove_text=True)
    def test_mesh_render_default_engine_tet4(self):
        fig, ax = plt.subplots()
        ax.imshow(self.images["connect1_u"])
        return fig


@pytest.mark.usefixtures("with_pyembree_ray_tracing_engine")
class TestTet4PyembreeEngine:
    @requires_file(tet4)
    @classmethod
    def setup_class(cls):
        cls.ds = data_dir_load_v2(tet4, step=-1)
        cls.images = {}
        for field in [("connect1", "u")]:
            sc = create_scene(cls.ds, field)
            im = sc.render()
            cls.images["_".join(field)] = im

    @pytest.mark.mpl_image_compare(remove_text=True)
    def test_mesh_render_pyembree_engine_tet4(self):
        fig, ax = plt.subplots()
        ax.imshow(self.images["connect1_u"])
        return fig


hex20 = "MOOSE_sample_data/mps_out.e"


@pytest.mark.usefixtures("with_default_ray_tracing_engine")
class TestHex20DefaultEngine:
    @requires_file(hex20)
    @classmethod
    def setup_class(cls):
        cls.ds = data_dir_load_v2(hex20, step=-1)
        cls.images = {}
        for field in [("connect2", "temp")]:
            sc = create_scene(cls.ds, field)
            im = sc.render()
            cls.images["_".join(field)] = im

    @pytest.mark.mpl_image_compare(remove_text=True)
    def test_mesh_render_default_engine_hex20(self):
        fig, ax = plt.subplots()
        ax.imshow(self.images["connect2_temp"])
        return fig


@pytest.mark.usefixtures("with_pyembree_ray_tracing_engine")
class TestHex20PyembreeEngine:
    @requires_file(hex20)
    @classmethod
    def setup_class(cls):
        cls.ds = data_dir_load_v2(hex20, step=-1)
        cls.images = {}
        for field in [("connect2", "temp")]:
            sc = create_scene(cls.ds, field)
            im = sc.render()
            cls.images["_".join(field)] = im

    @pytest.mark.mpl_image_compare(remove_text=True)
    def test_mesh_render_pyembree_engine_hex20(self):
        fig, ax = plt.subplots()
        ax.imshow(self.images["connect2_temp"])
        return fig


wedge6 = "MOOSE_sample_data/wedge_out.e"


@pytest.mark.usefixtures("with_default_ray_tracing_engine")
class TestWedge6DefaultEngine:
    @requires_file(wedge6)
    @classmethod
    def setup_class(cls):
        cls.ds = data_dir_load_v2(wedge6, step=-1)
        cls.images = {}
        for field in [("connect1", "diffused")]:
            sc = create_scene(cls.ds, field)
            im = sc.render()
            cls.images["_".join(field)] = im

    @pytest.mark.mpl_image_compare(remove_text=True)
    def test_mesh_render_default_engine_wedge6(self):
        fig, ax = plt.subplots()
        ax.imshow(self.images["connect1_diffused"])
        return fig


@pytest.mark.usefixtures("with_pyembree_ray_tracing_engine")
class TestWedge6PyembreeEngine:
    @requires_file(wedge6)
    @classmethod
    def setup_class(cls):
        cls.ds = data_dir_load_v2(wedge6, step=-1)
        cls.images = {}
        for field in [("connect1", "diffused")]:
            sc = create_scene(cls.ds, field)
            im = sc.render()
            cls.images["_".join(field)] = im

    @pytest.mark.mpl_image_compare(remove_text=True)
    def test_mesh_render_pyembree_engine_wedge6(self):
        fig, ax = plt.subplots()
        ax.imshow(self.images["connect1_diffused"])
        return fig


tet10 = "SecondOrderTets/tet10_unstructured_out.e"


@pytest.mark.usefixtures("with_default_ray_tracing_engine")
class TestTet10DefaultEngine:
    @requires_file(tet10)
    @classmethod
    def setup_class(cls):
        cls.ds = data_dir_load_v2(tet10, step=-1)
        cls.images = {}
        for field in [("connect1", "uz")]:
            sc = create_scene(cls.ds, field)
            im = sc.render()
            cls.images["_".join(field)] = im

    @pytest.mark.mpl_image_compare(remove_text=True)
    def test_mesh_render_default_engine_tet10(self):
        fig, ax = plt.subplots()
        ax.imshow(self.images["connect1_uz"])
        return fig


@pytest.mark.usefixtures("with_pyembree_ray_tracing_engine")
class TestTet10PyembreeEngine:
    @requires_file(tet10)
    @classmethod
    def setup_class(cls):
        cls.ds = data_dir_load_v2(tet10, step=-1)
        cls.images = {}
        for field in [("connect1", "uz")]:
            sc = create_scene(cls.ds, field)
            im = sc.render()
            cls.images["_".join(field)] = im

    @pytest.mark.mpl_image_compare(remove_text=True)
    def test_mesh_render_pyembree_engine_tet10(self):
        fig, ax = plt.subplots()
        ax.imshow(self.images["connect1_uz"])
        return fig


@requires_file(hex8)
@pytest.mark.usefixtures("with_default_ray_tracing_engine")
@pytest.mark.mpl_image_compare
def test_perspective_mesh_render_default():
    ds = data_dir_load(hex8)
    sc = create_scene(ds, ("connect2", "diffused"))
    cam = sc.add_camera(ds, lens_type="perspective")
    cam.focus = ds.arr([0.0, 0.0, 0.0], "code_length")
    cam_pos = ds.arr([-4.5, 4.5, -4.5], "code_length")
    north_vector = ds.arr([0.0, -1.0, -1.0], "dimensionless")
    cam.set_position(cam_pos, north_vector)
    cam.resolution = (800, 800)
    im = sc.render()

    fig, ax = plt.subplots()
    ax.imshow(im)
    return fig


@requires_file(hex8)
@pytest.mark.usefixtures("with_pyembree_ray_tracing_engine")
@pytest.mark.mpl_image_compare
def test_perspective_mesh_render_pyembree():
    ds = data_dir_load(hex8)
    sc = create_scene(ds, ("connect2", "diffused"))
    cam = sc.add_camera(ds, lens_type="perspective")
    cam.focus = ds.arr([0.0, 0.0, 0.0], "code_length")
    cam_pos = ds.arr([-4.5, 4.5, -4.5], "code_length")
    north_vector = ds.arr([0.0, -1.0, -1.0], "dimensionless")
    cam.set_position(cam_pos, north_vector)
    cam.resolution = (800, 800)
    im = sc.render()

    fig, ax = plt.subplots()
    ax.imshow(im)
    return fig


@requires_file(hex8)
@pytest.mark.usefixtures("with_default_ray_tracing_engine")
@pytest.mark.mpl_image_compare
def test_composite_mesh_render_default():
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

    fig, ax = plt.subplots()
    ax.imshow(im)
    return fig


@requires_file(hex8)
@pytest.mark.usefixtures("with_pyembree_ray_tracing_engine")
@pytest.mark.mpl_image_compare
def test_composite_mesh_render_pyembree():
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

    fig, ax = plt.subplots()
    ax.imshow(im)
    return fig
