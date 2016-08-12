"""
Testsuite for PlotWindow class



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
import itertools
import numpy as np
import os
import tempfile
import shutil
import unittest
from nose.tools import assert_true
from yt.extern.parameterized import parameterized, param
from yt.testing import \
    fake_random_ds, assert_equal, assert_rel_equal, assert_array_equal, \
    assert_array_almost_equal, assert_raises
from yt.utilities.answer_testing.framework import \
    requires_ds, data_dir_load, PlotWindowAttributeTest
from yt.utilities.exceptions import \
    YTInvalidFieldType
from yt.visualization.api import \
    SlicePlot, ProjectionPlot, OffAxisSlicePlot, OffAxisProjectionPlot
from yt.units.yt_array import YTArray, YTQuantity
from yt.frontends.stream.api import load_uniform_grid
from collections import OrderedDict

def setup():
    """Test specific setup."""
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"


def assert_fname(fname):
    """Function that checks file type using libmagic"""
    if fname is None:
        return

    with open(fname, 'rb') as fimg:
        data = fimg.read()
    image_type = ''

    # see http://www.w3.org/TR/PNG/#5PNG-file-signature
    if data.startswith(b'\211PNG\r\n\032\n'):
        image_type = '.png'
    # see http://www.mathguide.de/info/tools/media-types/image/jpeg
    elif data.startswith(b'\377\330'):
        image_type = '.jpeg'
    elif data.startswith(b'%!PS-Adobe'):
        data_str = data.decode("utf-8", "ignore")
        if 'EPSF' in data_str[:data_str.index('\n')]:
            image_type = '.eps'
        else:
            image_type = '.ps'
    elif data.startswith(b'%PDF'):
        image_type = '.pdf'

    return image_type == os.path.splitext(fname)[1]


TEST_FLNMS = [None, 'test', 'test.png', 'test.eps',
              'test.ps', 'test.pdf']
M7 = "DD0010/moving7_0010"
WT = "WindTunnel/windtunnel_4lev_hdf5_plt_cnt_0030"

FPROPS = {'family': 'sans-serif', 'style': 'italic',
          'weight': 'bold', 'size': 24}

ATTR_ARGS = {"pan": [(((0.1, 0.1), ), {})],
             "pan_rel": [(((0.1, 0.1), ), {})],
             "set_axes_unit": [(("kpc", ), {}),
                               (("Mpc", ), {}),
                               ((("kpc", "kpc"),), {}),
                               ((("kpc", "Mpc"),), {})],
             "set_buff_size": [((1600, ), {}),
                               (((600, 800), ), {})],
             "set_center": [(((0.4, 0.3), ), {})],
             "set_cmap": [(('density', 'RdBu'), {}),
                          (('density', 'kamae'), {})],
             "set_font": [((OrderedDict(sorted(FPROPS.items(), key=lambda t: t[0])), ),
                           {})],
             "set_log": [(('density', False), {})],
             "set_window_size": [((7.0, ), {})],
             "set_zlim": [(('density', 1e-25, 1e-23), {}),
                          (('density', 1e-25, None), {'dynamic_range': 4})],
             "zoom": [((10, ), {})]}


CENTER_SPECS = (
    "m",
    "M",
    "max",
    "Max",
    "c",
    "C",
    "center",
    "Center",
    [0.5, 0.5, 0.5],
    [[0.2, 0.3, 0.4], "cm"],
    YTArray([0.3, 0.4, 0.7], "cm")
)

WIDTH_SPECS = {
    # Width choices map to xlim, ylim, width, axes_unit_name 4-tuples
    None   : (
        ((0, 'code_length'), (1, 'code_length')),
        ((0, 'code_length'), (1, 'code_length')),
        ((1, 'code_length'), (1, 'code_length')),
        None,
    ),
    0.2 : (
        ((0.4, 'code_length'), (0.6, 'code_length')),
        ((0.4, 'code_length'), (0.6, 'code_length')),
        ((0.2, 'code_length'), (0.2, 'code_length')),
        None,
    ),
    (0.4, 0.3) : (
        ((0.3, 'code_length'), (0.7, 'code_length')),
        ((0.35, 'code_length'), (0.65, 'code_length')),
        ((0.4, 'code_length'), (0.3, 'code_length')),
        None,
    ),
    (1.2, 'cm') : (
        ((-0.1, 'code_length'), (1.1, 'code_length')),
        ((-0.1, 'code_length'), (1.1, 'code_length')),
        ((1.2,  'code_length'), (1.2, 'code_length')),
        ('cm', 'cm'),
    ),
    ((1.2, 'cm'), (2.0, 'cm')) : (
        ((-0.1, 'code_length'), (1.1, 'code_length')),
        ((-0.5, 'code_length'), (1.5, 'code_length')),
        ((1.2,  'code_length'), (2.0, 'code_length')),
        ('cm', 'cm'),
    ),
    ((1.2, 'cm'), (0.02, 'm')) : (
        ((-0.1, 'code_length'), (1.1, 'code_length')),
        ((-0.5, 'code_length'), (1.5, 'code_length')),
        ((1.2,  'code_length'), (2.0, 'code_length')),
        ('cm', 'm'),
    ),
}

WEIGHT_FIELDS = (
    None,
    'density',
    ('gas', 'density'),
)

PROJECTION_METHODS = (
    'integrate',
    'sum',
    'mip'
)

def simple_contour(test_obj, plot):
    plot.annotate_contour(test_obj.plot_field)

def simple_velocity(test_obj, plot):
    plot.annotate_velocity()

def simple_streamlines(test_obj, plot):
    ax = test_obj.plot_axis
    xax = test_obj.ds.coordinates.x_axis[ax]
    yax = test_obj.ds.coordinates.y_axis[ax]
    xn = test_obj.ds.coordinates.axis_name[xax]
    yn = test_obj.ds.coordinates.axis_name[yax]
    plot.annotate_streamlines("velocity_%s" % xn,
                              "velocity_%s" % yn)

CALLBACK_TESTS = (
    ("simple_contour", (simple_contour,)),
    ("simple_velocity", (simple_velocity,)),
    #("simple_streamlines", (simple_streamlines,)),
    #("simple_all", (simple_contour, simple_velocity, simple_streamlines)),
)

@requires_ds(M7)
def test_attributes():
    """Test plot member functions that aren't callbacks"""
    plot_field = 'density'
    decimals = 12

    ds = data_dir_load(M7)
    for ax in 'xyz':
        for attr_name in ATTR_ARGS.keys():
            for args in ATTR_ARGS[attr_name]:
                test = PlotWindowAttributeTest(ds, plot_field, ax, attr_name,
                                               args, decimals)
                test_attributes.__name__ = test.description
                yield test
                for n, r in CALLBACK_TESTS:
                    yield PlotWindowAttributeTest(ds, plot_field, ax, attr_name,
                        args, decimals, callback_id=n, callback_runners=r)


@requires_ds(WT)
def test_attributes_wt():
    plot_field = 'density'
    decimals = 12

    ds = data_dir_load(WT)
    ax = 'z'
    for attr_name in ATTR_ARGS.keys():
        for args in ATTR_ARGS[attr_name]:
            yield PlotWindowAttributeTest(ds, plot_field, ax, attr_name,
                                          args, decimals)
            for n, r in CALLBACK_TESTS:
                yield PlotWindowAttributeTest(ds, plot_field, ax, attr_name,
                    args, decimals, callback_id=n, callback_runners=r)

class TestHideAxesColorbar(unittest.TestCase):

    ds = None

    def setUp(self):
        if self.ds is None:
            self.ds = fake_random_ds(64)
            self.slc = SlicePlot(self.ds, 0, "density")
        self.tmpdir = tempfile.mkdtemp()
        self.curdir = os.getcwd()
        os.chdir(self.tmpdir)

    def tearDown(self):
        os.chdir(self.curdir)
        shutil.rmtree(self.tmpdir)

    def test_hide_show_axes(self):
        self.slc.hide_axes()
        self.slc.save()
        self.slc.show_axes()
        self.slc.save()

    def test_hide_show_colorbar(self):
        self.slc.hide_colorbar()
        self.slc.save()
        self.slc.show_colorbar()
        self.slc.save()

    def test_hide_axes_colorbar(self):
        self.slc.hide_colorbar()
        self.slc.hide_axes()
        self.slc.save()

class TestSetWidth(unittest.TestCase):

    ds = None

    def setUp(self):
        if self.ds is None:
            self.ds = fake_random_ds(64)
            self.slc = SlicePlot(self.ds, 0, "density")

    def _assert_05cm(self):
        assert_array_equal([self.slc.xlim, self.slc.ylim, self.slc.width],
                         [(YTQuantity(0.25, 'cm'), YTQuantity(0.75, 'cm')),
                          (YTQuantity(0.25, 'cm'), YTQuantity(0.75, 'cm')),
                          (YTQuantity(0.5,  'cm'), YTQuantity(0.5,  'cm'))])

    def _assert_05_075cm(self):
        assert_array_equal([self.slc.xlim, self.slc.ylim, self.slc.width],
                         [(YTQuantity(0.25,  'cm'), YTQuantity(0.75,  'cm')),
                          (YTQuantity(0.125, 'cm'), YTQuantity(0.875, 'cm')),
                          (YTQuantity(0.5,   'cm'), YTQuantity(0.75,  'cm'))])

    def test_set_width_one(self):
        assert_equal([self.slc.xlim, self.slc.ylim, self.slc.width],
                     [(0.0, 1.0), (0.0, 1.0), (1.0, 1.0)])
        assert_true(self.slc._axes_unit_names is None)

    def test_set_width_nonequal(self):
        self.slc.set_width((0.5, 0.8))
        assert_rel_equal([self.slc.xlim, self.slc.ylim, self.slc.width],
                         [(0.25, 0.75), (0.1, 0.9), (0.5, 0.8)], 15)
        assert_true(self.slc._axes_unit_names is None)

    def test_twoargs_eq(self):
        self.slc.set_width(0.5, 'cm')
        self._assert_05cm()
        assert_true(self.slc._axes_unit_names == ('cm', 'cm'))

    def test_tuple_eq(self):
        self.slc.set_width((0.5, 'cm'))
        self._assert_05cm()
        assert_true(self.slc._axes_unit_names == ('cm', 'cm'))

    def test_tuple_of_tuples_neq(self):
        self.slc.set_width(((0.5, 'cm'), (0.75, 'cm')))
        self._assert_05_075cm()
        assert_true(self.slc._axes_unit_names == ('cm', 'cm'))


class TestPlotWindowSave(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        test_ds = fake_random_ds(64)
        normal = [1, 1, 1]
        ds_region = test_ds.region([0.5] * 3, [0.4] * 3, [0.6] * 3)
        projections = []
        projections_ds = []
        projections_c = []
        projections_wf = []
        projections_w = {}
        projections_m = []
        for dim in range(3):
            projections.append(ProjectionPlot(test_ds, dim, "density"))
            projections_ds.append(ProjectionPlot(test_ds, dim, "density",
                                                 data_source=ds_region))
        for center in CENTER_SPECS:
            projections_c.append(ProjectionPlot(test_ds, dim, "density",
                                                center=center))
        for width in WIDTH_SPECS:
            projections_w[width] = ProjectionPlot(test_ds, dim, 'density',
                                                  width=width)
        for wf in WEIGHT_FIELDS:
            projections_wf.append(ProjectionPlot(test_ds, dim, "density",
                                                 weight_field=wf))
        for m in PROJECTION_METHODS:
            projections_m.append(ProjectionPlot(test_ds, dim, "density",
                                                 method=m))

        cls.slices = [SlicePlot(test_ds, dim, "density") for dim in range(3)]
        cls.projections = projections
        cls.projections_ds = projections_ds
        cls.projections_c = projections_c
        cls.projections_wf = projections_wf
        cls.projections_w = projections_w
        cls.projections_m = projections_m
        cls.offaxis_slice = OffAxisSlicePlot(test_ds, normal, "density")
        cls.offaxis_proj = OffAxisProjectionPlot(test_ds, normal, "density")

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.curdir = os.getcwd()
        os.chdir(self.tmpdir)

    def tearDown(self):
        os.chdir(self.curdir)
        shutil.rmtree(self.tmpdir)

    @parameterized.expand(
        param.explicit(item)
        for item in itertools.product(range(3), TEST_FLNMS))
    def test_slice_plot(self, dim, fname):
        assert assert_fname(self.slices[dim].save(fname)[0])

    @parameterized.expand(
        param.explicit(item)
        for item in itertools.product(range(3), TEST_FLNMS))
    def test_projection_plot(self, dim, fname):
        assert assert_fname(self.projections[dim].save(fname)[0])

    @parameterized.expand([(0, ), (1, ), (2, )])
    def test_projection_plot_ds(self, dim):
        self.projections_ds[dim].save()

    @parameterized.expand([(i, ) for i in range(len(CENTER_SPECS))])
    def test_projection_plot_c(self, dim):
        self.projections_c[dim].save()

    @parameterized.expand([(i, ) for i in range(len(WEIGHT_FIELDS))])
    def test_projection_plot_wf(self, dim):
        self.projections_wf[dim].save()

    @parameterized.expand([(i, ) for i in range(len(PROJECTION_METHODS))])
    def test_projection_plot_m(self, dim):
        self.projections_m[dim].save()

    @parameterized.expand(
        param.explicit((fname, ))
        for fname in TEST_FLNMS)
    def test_offaxis_slice_plot(self, fname):
        assert assert_fname(self.offaxis_slice.save(fname)[0])

    @parameterized.expand(
        param.explicit((fname, ))
        for fname in TEST_FLNMS)
    def test_offaxis_projection_plot(self, fname):
        assert assert_fname(self.offaxis_proj.save(fname)[0])

    def test_ipython_repr(self):
        self.slices[0]._repr_html_()

    @parameterized.expand(
        param.explicit((width, ))
        for width in WIDTH_SPECS)
    def test_creation_with_width(self, width):
        xlim, ylim, pwidth, aun = WIDTH_SPECS[width]
        plot = self.projections_w[width]

        xlim = [plot.ds.quan(el[0], el[1]) for el in xlim]
        ylim = [plot.ds.quan(el[0], el[1]) for el in ylim]
        pwidth = [plot.ds.quan(el[0], el[1]) for el in pwidth]

        [assert_array_almost_equal(px, x, 14) for px, x in zip(plot.xlim, xlim)]
        [assert_array_almost_equal(py, y, 14) for py, y in zip(plot.ylim, ylim)]
        [assert_array_almost_equal(pw, w, 14) for pw, w in zip(plot.width, pwidth)]
        assert_true(aun == plot._axes_unit_names)

def test_on_off_compare():
    # fake density field that varies in the x-direction only
    den = np.arange(32**3) / 32**2 + 1
    den = den.reshape(32, 32, 32)
    den = np.array(den, dtype=np.float64)
    data = dict(density = (den, "g/cm**3"))
    bbox = np.array([[-1.5, 1.5], [-1.5, 1.5], [-1.5, 1.5]])
    ds = load_uniform_grid(data, den.shape, length_unit="Mpc", bbox=bbox, nprocs=64)

    sl_on = SlicePlot(ds, "z", ["density"])

    L = [0, 0, 1]
    north_vector = [0, 1, 0]
    sl_off = OffAxisSlicePlot(ds, L, 'density', center=[0,0,0], north_vector=north_vector)

    assert_array_almost_equal(sl_on.frb['density'], sl_off.frb['density'])

def test_plot_particle_field_error():
    ds = fake_random_ds(32, particles=100)

    field_names = [
        'particle_mass',
        ['particle_mass', 'density'],
        ['density', 'particle_mass'],
    ]

    objects_normals = [
        (SlicePlot, 2),
        (SlicePlot, [1, 1, 1]),
        (ProjectionPlot, 2),
        (OffAxisProjectionPlot, [1, 1, 1]),
    ]

    for object, normal in objects_normals:
        for field_name_list in field_names:
            assert_raises(
                YTInvalidFieldType, object, ds, normal, field_name_list)


def test_frb_regen():
    ds = fake_random_ds(32)
    slc = SlicePlot(ds, 2, 'density')
    slc.set_buff_size(1200)
    assert_equal(slc.frb['density'].shape, (1200, 1200))
