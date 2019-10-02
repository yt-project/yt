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
from collections import OrderedDict
import os
import tempfile

import pytest

from yt.config import \
    ytcfg
from yt.testing import \
    fake_tetrahedral_ds, \
    fake_hexahedral_ds, \
    requires_module
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils
from yt.visualization.volume_rendering.api import \
    MeshSource, \
    Scene, \
    create_scene


hex8 = "MOOSE_sample_data/out.e-s010"
hex8_fields = [('connect1', 'diffused'), ('connect2', 'convected')]
tet4 = "MOOSE_sample_data/high_order_elems_tet4_refine_out.e"
tet4_fields = [("connect1", "u")]
hex20 = "MOOSE_sample_data/mps_out.e"
hex20_fields = [('connect2', 'temp')]
wedge6 = "MOOSE_sample_data/wedge_out.e"
wedge6_fields = [('connect1', 'diffused')]
tet10 = "SecondOrderTets/tet10_unstructured_out.e"
tet10_fields = [('connect1', 'uz')]


def surface_mesh_render():
    images = []
    ds = fake_tetrahedral_ds()
    for field in ds.field_list:
        if field[0] == 'all':
            continue
        sc = Scene()
        sc.add_source(MeshSource(ds, field))
        sc.add_camera()
        im = sc.render()
        images.append(im)
    ds = fake_hexahedral_ds()
    for field in ds.field_list:
        if field[0] == 'all':
            continue
        sc = Scene()
        sc.add_source(MeshSource(ds, field))
        sc.add_camera()
        im = sc.render()
        images.append(im)
    return images

def hex8_render(engine, field):
    ytcfg["yt", "ray_tracing_engine"] = engine
    ds = utils.data_dir_load(hex8, kwargs={'step':-1})
    sc = create_scene(ds, field)
    im = sc.render()
    fd, im_name = tempfile.mkstemp(suffix='.png', prefix='tmp', dir=os.getcwd())
    im.save(im_name)
    return im_name


def tet4_render(engine, field):
    ytcfg["yt", "ray_tracing_engine"] = engine
    ds = utils.data_dir_load(tet4, kwargs={'step':-1})
    sc = create_scene(ds, field)
    im = sc.render()
    fd, im_name = tempfile.mkstemp(suffix='.png', prefix='tmp', dir=os.getcwd())
    im.save(im_name)
    return im_name


def hex20_render(engine, field):
    ytcfg["yt", "ray_tracing_engine"] = engine
    ds = utils.data_dir_load(hex20, kwargs={'step':-1})
    sc = create_scene(ds, field)
    im = sc.render()
    fd, im_name = tempfile.mkstemp(suffix='.png', prefix='tmp', dir=os.getcwd())
    im.save(im_name)
    return im_name


def wedge6_render(engine, field):
    ytcfg["yt", "ray_tracing_engine"] = engine
    ds = utils.data_dir_load(wedge6, kwargs={'step':-1})
    sc = create_scene(ds, field)
    im = sc.render()
    fd, im_name = tempfile.mkstemp(suffix='.png', prefix='tmp', dir=os.getcwd())
    im.save(im_name)
    return im_name


def tet10_render(engine, field):
    ytcfg["yt", "ray_tracing_engine"] = engine
    ds = utils.data_dir_load(tet10, kwargs={'step':-1})
    sc = create_scene(ds, field)
    ms = sc.get_source(0)
    ms.color_bounds = (-.01, .2)
    im = sc.render()
    fd, im_name = tempfile.mkstemp(suffix='.png', prefix='tmp', dir=os.getcwd())
    im.save(im_name)
    return im_name


def perspective_mesh_render(engine):
    ytcfg["yt", "ray_tracing_engine"] = engine
    ds = utils.data_dir_load(hex8)
    sc = create_scene(ds, ("connect2", "diffused"))
    cam = sc.add_camera(ds, lens_type='perspective')
    cam.focus = ds.arr([0.0, 0.0, 0.0], 'code_length')
    cam_pos = ds.arr([-4.5, 4.5, -4.5], 'code_length')
    north_vector = ds.arr([0.0, -1.0, -1.0], 'dimensionless')
    cam.set_position(cam_pos, north_vector)
    cam.resolution = (800, 800)
    im = sc.render()
    fd, im_name = tempfile.mkstemp(suffix='.png', prefix='tmp', dir=os.getcwd())
    im.save(im_name)
    return im_name


def composite_mesh_render(engine):
    ytcfg["yt", "ray_tracing_engine"] = engine
    ds = utils.data_dir_load(hex8)
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
    fd, im_name = tempfile.mkstemp(suffix='.png', prefix='tmp', dir=os.getcwd())
    im.save(im_name)
    return im_name

@pytest.mark.skipif(not pytest.config.getvalue('--with-answer-testing'),
    reason="--with-answer-testing not set.")
@pytest.mark.usefixtures('temp_dir', 'answer_file')
class TestVolumeRenderMesh(fw.AnswerTest):
    @pytest.mark.usefixtures('hashing')
    def test_fake_hexahedral_ds_render(self):
        ds = fake_hexahedral_ds()
        field_list = [('connect1', 'elem'), ('connect1', 'test')]
        self.hashes['generic_image'] = OrderedDict()
        for field in field_list:
            fd, im_name = tempfile.mkstemp(suffix='.png', prefix='tmp', dir=os.getcwd())
            sc = create_scene(ds, field)
            im = sc.render()
            im.save(im_name)
            gi_hd = self.generic_image_test(im_name)
            self.hashes['generic_image'][field] = gi_hd

    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(hex8)
    @requires_module("pyembree")
    def test_composite_mesh_render_pyembree(self):
        im_name = composite_mesh_render("embree")
        gi_hd = self.generic_image_test(im_name)
        self.hashes['generic_image'] = gi_hd

    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(hex8)
    def test_composite_mesh_render(self):
        im_name = composite_mesh_render("yt")
        gi_hd = self.generic_image_test(im_name)
        self.hashes['generic_image'] = gi_hd

    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(hex20)
    @requires_module("pyembree")
    def test_hex20_render_pyembree(self):
        self.hashes['generic_image'] = OrderedDict()
        for field in hex20_fields:
            im_name = hex20_render("embree", field)
            gi_hd = self.generic_image_test(im_name)
            self.hashes['generic_image'][field] = gi_hd

    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(hex20)
    def test_hex20_render(self):
        self.hashes['generic_image'] = OrderedDict()
        for field in hex20_fields:
            im_name = hex20_render("yt", field)
            gi_hd = self.generic_image_test(im_name)
            self.hashes['generic_image'][field] = gi_hd

    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(hex8)
    @requires_module("pyembree")
    def test_hex8_render_pyembree(self):
        self.hashes['generic_image'] = OrderedDict()
        for field in hex8_fields:
            im_name = hex8_render("embree", field)
            gi_hd = self.generic_image_test(im_name)
            self.hashes['generic_image'][field] = gi_hd

    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(hex8)
    @utils.requires_ds(hex8)
    def test_hex8_render(self):
        self.hashes['generic_image'] = OrderedDict()
        for field in hex8_fields:
            im_name = hex8_render("yt", field)
            gi_hd = self.generic_image_test(im_name)
            self.hashes['generic_image'][field] = gi_hd

    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(hex8)
    @requires_module("pyembree")
    def test_perspective_mesh_render_pyembree(self):
        im_name = perspective_mesh_render("embree")
        gi_hd = self.generic_image_test(im_name)
        self.hashes['generic_image'] = gi_hd

    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(hex8)
    def test_perspective_mesh_render(self):
        im_name = perspective_mesh_render("yt")
        gi_hd = self.generic_image_test(im_name)
        self.hashes['generic_image'] = gi_hd

    @requires_module("pyembree")
    def test_surface_mesh_render_pyembree(self):
        ytcfg["yt", "ray_tracing_engine"] = "embree"
        surface_mesh_render()

    def test_surface_mesh_render(self):
        ytcfg["yt", "ray_tracing_engine"] = "yt"
        surface_mesh_render()

    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(tet10)
    @requires_module("pyembree")
    def test_tet10_render_pyembree(self):
        self.hashes['generic_image'] = OrderedDict()
        for field in tet10_fields:
            im_name = tet10_render("embree", field)
            gi_hd = self.generic_image_test(im_name)
            self.hashes['generic_image'][field] = gi_hd

    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(tet10)
    def test_tet10_render(self):
        self.hashes['generic_image'] = OrderedDict()
        for field in tet10_fields:
            im_name = tet10_render("yt", field)
            gi_hd = self.generic_image_test(im_name)
            self.hashes['generic_image'][field] = gi_hd

    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(tet4)
    @requires_module("pyembree")
    def test_tet4_render_pyembree(self):
        self.hashes['generic_image'] = OrderedDict()
        for field in tet4_fields:
            im_name = tet4_render("embree", field)
            gi_hd = self.generic_image_test(im_name)
            self.hashes['generic_image'][field] = gi_hd

    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(tet4)
    def test_tet4_render(self):
        self.hashes['generic_image'] = OrderedDict()
        for field in tet4_fields:
            im_name = tet4_render("yt", field)
            gi_hd = self.generic_image_test(im_name)
            self.hashes['generic_image'][field] = gi_hd

    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(wedge6)
    @requires_module("pyembree")
    def test_wedge6_render_pyembree(self):
        self.hashes['generic_image'] = OrderedDict()
        for field in wedge6_fields:
            im_name = wedge6_render("embree", field)
            gi_hd = self.generic_image_test(im_name)
            self.hashes['generic_image'][field] = gi_hd

    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(wedge6)
    def test_wedge6_render(self):
        self.hashes['generic_image'] = OrderedDict()
        for field in wedge6_fields:
            im_name = wedge6_render("yt", field)
            gi_hd = self.generic_image_test(im_name)
            self.hashes['generic_image'][field] = gi_hd
