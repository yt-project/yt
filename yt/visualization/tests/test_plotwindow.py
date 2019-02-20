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
import matplotlib
import numpy as np
import os
import shutil
import tempfile
import unittest

from distutils.version import LooseVersion
from nose.tools import assert_true

from yt.testing import \
    fake_random_ds, assert_equal, assert_rel_equal, assert_array_equal, \
    assert_array_almost_equal, assert_raises, assert_fname
from yt.utilities.answer_testing.framework import \
    requires_ds, data_dir_load, PlotWindowAttributeTest
from yt.utilities.exceptions import \
    YTInvalidFieldType
from yt.visualization.api import \
    SlicePlot, ProjectionPlot, OffAxisSlicePlot, OffAxisProjectionPlot, \
    plot_2d
from yt.units.yt_array import YTArray, YTQuantity
from yt.units import kboltz
from yt.frontends.stream.api import load_uniform_grid
from collections import OrderedDict

def setup():
    """Test specific setup."""
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"


TEST_FLNMS = ['test.png']
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
             "zoom": [((10, ), {})],
             "toggle_right_handed": [((),{})]}


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
        del self.ds
        del self.slc

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

    def tearDown(self):
        del self.ds
        del self.slc

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
        
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.curdir = os.getcwd()
        os.chdir(self.tmpdir)

    def tearDown(self):
        os.chdir(self.curdir)
        shutil.rmtree(self.tmpdir)

    def test_slice_plot(self):
        test_ds = fake_random_ds(16)
        for dim in range(3):
            slc = SlicePlot(test_ds, dim, 'density')
            for fname in TEST_FLNMS:
                assert_fname(slc.save(fname)[0])

    def test_repr_html(self):
        test_ds = fake_random_ds(16)
        slc = SlicePlot(test_ds, 0, 'density')
        slc._repr_html_()

    def test_projection_plot(self):
        test_ds = fake_random_ds(16)
        for dim in range(3):
            proj = ProjectionPlot(test_ds, dim, 'density')
            for fname in TEST_FLNMS:
                assert_fname(proj.save(fname)[0])

    def test_projection_plot_ds(self):
        test_ds = fake_random_ds(16)
        reg = test_ds.region([0.5] * 3, [0.4] * 3, [0.6] * 3)
        for dim in range(3):
            proj = ProjectionPlot(test_ds, dim, 'density', data_source=reg)
            proj.save()

    def test_projection_plot_c(self):
        test_ds = fake_random_ds(16)
        for center in CENTER_SPECS:
            proj = ProjectionPlot(test_ds, 0, 'density', center=center)
            proj.save()

    def test_projection_plot_wf(self):
        test_ds = fake_random_ds(16)
        for wf in WEIGHT_FIELDS:
            proj = ProjectionPlot(test_ds, 0, 'density', weight_field=wf)
            proj.save()

    def test_projection_plot_m(self):
        test_ds = fake_random_ds(16)
        for method in PROJECTION_METHODS:
            proj = ProjectionPlot(test_ds, 0, 'density', method=method)
            proj.save()

    def test_offaxis_slice_plot(self):
        test_ds = fake_random_ds(16)
        slc = OffAxisSlicePlot(test_ds, [1, 1, 1], "density")
        for fname in TEST_FLNMS:
            assert_fname(slc.save(fname)[0])

    def test_offaxis_projection_plot(self):
        test_ds = fake_random_ds(16)
        prj = OffAxisProjectionPlot(test_ds, [1, 1, 1], "density")
        for fname in TEST_FLNMS:
            assert_fname(prj.save(fname)[0])

    def test_creation_with_width(self):
        test_ds = fake_random_ds(16)
        for width in WIDTH_SPECS:
            xlim, ylim, pwidth, aun = WIDTH_SPECS[width]
            plot = ProjectionPlot(test_ds, 0, 'density', width=width)

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

    sl_on.set_buff_size((800, 400))
    sl_on._recreate_frb()
    sl_off.set_buff_size((800, 400))
    sl_off._recreate_frb()

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

def test_setup_origin():
    origin_inputs = ('domain',
                     'left-window',
                     'center-domain',
                     'lower-right-window',
                     ('window',),
                     ('right', 'domain'),
                     ('lower', 'window'),
                     ('lower', 'right', 'window'),
                     (0.5, 0.5, 'domain'),
                     ((50, 'cm'), (50, 'cm'), 'domain'))
    w=(10, 'cm')

    ds = fake_random_ds(32, length_unit=100.0)
    generated_limits = []
    #lower limit -> llim
    #upper limit -> ulim
    #                 xllim xulim yllim yulim
    correct_limits = [45.0, 55.0, 45.0, 55.0,
                      0.0, 10.0, 0.0, 10.0,
                      -5.0, 5.0, -5.0, 5.0,
                      -10.0, 0, 0, 10.0,
                      0.0, 10.0, 0.0, 10.0,
                      -55.0, -45.0, -55.0, -45.0,
                      -5.0, 5.0, 0.0, 10.0,
                      -10.0, 0, 0, 10.0,
                      -5.0, 5.0, -5.0, 5.0,
                      -5.0, 5.0, -5.0, 5.0
                      ]
    for o in origin_inputs:
        slc = SlicePlot(ds, 2, 'density', width=w, origin=o)
        ax = slc.plots['density'].axes
        xlims = ax.get_xlim()
        ylims = ax.get_ylim()
        lims = [xlims[0], xlims[1], ylims[0], ylims[1]]
        for l in lims:
            generated_limits.append(l)
    assert_array_almost_equal(correct_limits, generated_limits)

def test_frb_regen():
    ds = fake_random_ds(32)
    slc = SlicePlot(ds, 2, 'density')
    slc.set_buff_size(1200)
    assert_equal(slc.frb['density'].shape, (1200, 1200))
    slc.set_buff_size((400., 200.7))
    assert_equal(slc.frb['density'].shape, (200, 400))


def test_set_background_color():
    ds = fake_random_ds(32)
    plot = SlicePlot(ds, 2, 'density')
    for field in ['density', ('gas', 'density')]:
        plot.set_background_color(field, 'red')
        plot._setup_plots()
        ax = plot.plots[field].axes
        if LooseVersion(matplotlib.__version__) < LooseVersion('2.0.0'):
            assert_equal(ax.get_axis_bgcolor(), 'red')
        else:
            assert_equal(ax.get_facecolor(), (1.0, 0.0, 0.0, 1.0))

def test_set_unit():
    ds = fake_random_ds(32, fields=('temperature',), units=('K',))
    slc = SlicePlot(ds, 2, 'temperature')

    orig_array = slc.frb['gas', 'temperature'].copy()

    slc.set_unit('temperature', 'degF')

    assert str(slc.frb['gas', 'temperature'].units) == 'degF'
    assert_array_almost_equal(np.array(slc.frb['gas', 'temperature']),
                              np.array(orig_array)*1.8 - 459.67)

    # test that a plot modifying function that destroys the frb preserves the
    # new unit
    slc.set_buff_size(1000)

    assert str(slc.frb['gas', 'temperature'].units) == 'degF'

    slc.set_buff_size(800)

    slc.set_unit('temperature', 'K')
    assert str(slc.frb['gas', 'temperature'].units) == 'K'
    assert_array_almost_equal(slc.frb['gas', 'temperature'], orig_array)

    slc.set_unit('temperature', 'keV', equivalency='thermal')
    assert str(slc.frb['gas', 'temperature'].units) == 'keV'
    assert_array_almost_equal(slc.frb['gas', 'temperature'],
                              (orig_array*kboltz).to('keV'))

    # test that a plot modifying function that destroys the frb preserves the
    # new unit with an equivalency
    slc.set_buff_size(1000)

    assert str(slc.frb['gas', 'temperature'].units) == 'keV'

    # test that destroying the FRB then changing the unit using an equivalency
    # doesn't error out, see issue #1316
    slc = SlicePlot(ds, 2, 'temperature')
    slc.set_buff_size(1000)
    slc.set_unit('temperature', 'keV', equivalency='thermal')
    assert str(slc.frb['gas', 'temperature'].units) == 'keV'

WD = "WDMerger_hdf5_chk_1000/WDMerger_hdf5_chk_1000.hdf5"

@requires_ds(WD)
def test_plot_2d():
    # Cartesian
    ds = fake_random_ds((32,32,1), fields=('temperature',), units=('K',))
    slc = SlicePlot(ds, "z", ["temperature"], width=(0.2,"unitary"),
                    center=[0.4, 0.3, 0.5])
    slc2 = plot_2d(ds, "temperature", width=(0.2,"unitary"),
                   center=[0.4, 0.3])
    slc3 = plot_2d(ds, "temperature", width=(0.2,"unitary"),
                   center=ds.arr([0.4, 0.3], "cm"))
    assert_array_equal(slc.frb['temperature'], slc2.frb['temperature'])
    assert_array_equal(slc.frb['temperature'], slc3.frb['temperature'])
    # Cylindrical
    ds = data_dir_load(WD)
    slc = SlicePlot(ds, "theta", ["density"], width=(30000.0, "km"))
    slc2 = plot_2d(ds, "density", width=(30000.0, "km"))
    assert_array_equal(slc.frb['density'], slc2.frb['density'])

def test_symlog_colorbar():
    ds = fake_random_ds(16)

    def _thresh_density(field, data):
        wh = data['density'] < 0.5
        ret = data['density']
        ret[wh] = 0
        return ret

    def _neg_density(field, data):
        return -data['threshold_density']

    ds.add_field('threshold_density', function=_thresh_density,
                 units='g/cm**3', sampling_type='cell')
    ds.add_field('negative_density', function=_neg_density, units='g/cm**3',
                 sampling_type='cell')

    for field in ['density', 'threshold_density', 'negative_density']:
        plot = SlicePlot(ds, 2, field)
        plot.set_log(field, True, linthresh=0.1)
        with tempfile.NamedTemporaryFile(suffix='png') as f:
            plot.save(f.name)

def test_nan_data():
    data = np.random.random((16, 16, 16)) - 0.5
    data[:9, :9, :9] = np.nan

    data = {'density': data}

    ds = load_uniform_grid(data, [16, 16, 16])

    plot = SlicePlot(ds, 'z', 'density')

    with tempfile.NamedTemporaryFile(suffix='png') as f:
        plot.save(f.name)
