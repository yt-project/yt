import os
import shutil
import tempfile
import unittest
from collections import OrderedDict

import numpy as np
from matplotlib.colors import LogNorm, Normalize, SymLogNorm
from nose.tools import assert_true
from unyt import unyt_array

from yt.loaders import load_uniform_grid
from yt.testing import (
    assert_allclose_units,
    assert_array_almost_equal,
    assert_array_equal,
    assert_equal,
    assert_fname,
    assert_raises,
    assert_rel_equal,
    fake_amr_ds,
    fake_random_ds,
    requires_file,
    requires_module,
)
from yt.units import kboltz
from yt.units.yt_array import YTArray, YTQuantity
from yt.utilities.answer_testing.framework import (
    PlotWindowAttributeTest,
    data_dir_load,
    requires_ds,
)
from yt.utilities.exceptions import YTInvalidFieldType
from yt.visualization.plot_window import (
    AxisAlignedProjectionPlot,
    AxisAlignedSlicePlot,
    NormalPlot,
    OffAxisProjectionPlot,
    OffAxisSlicePlot,
    ProjectionPlot,
    SlicePlot,
    plot_2d,
)


def setup():
    """Test specific setup."""
    from yt.config import ytcfg

    ytcfg["yt", "internals", "within_testing"] = True


TEST_FLNMS = ["test.png"]
M7 = "DD0010/moving7_0010"

FPROPS = {"family": "sans-serif", "style": "italic", "weight": "bold", "size": 24}

ATTR_ARGS = {
    "pan": [(((0.1, 0.1),), {})],
    "pan_rel": [(((0.1, 0.1),), {})],
    "set_axes_unit": [
        (("kpc",), {}),
        (("Mpc",), {}),
        ((("kpc", "kpc"),), {}),
        ((("kpc", "Mpc"),), {}),
    ],
    "set_buff_size": [((1600,), {}), (((600, 800),), {})],
    "set_center": [(((0.4, 0.3),), {})],
    "set_cmap": [(("density", "RdBu"), {}), (("density", "cmyt.pastel"), {})],
    "set_font": [((OrderedDict(sorted(FPROPS.items(), key=lambda t: t[0])),), {})],
    "set_log": [(("density", False), {})],
    "set_figure_size": [((7.0,), {})],
    "set_zlim": [
        (("density", 1e-25, 1e-23), {}),
        (("density",), {"zmin": 1e-25, "dynamic_range": 4}),
    ],
    "zoom": [((10,), {})],
    "toggle_right_handed": [((), {})],
}


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
    YTArray([0.3, 0.4, 0.7], "cm"),
)

WIDTH_SPECS = {
    # Width choices map to xlim, ylim, width, axes_unit_name 4-tuples
    None: (
        ((0, "code_length"), (1, "code_length")),
        ((0, "code_length"), (1, "code_length")),
        ((1, "code_length"), (1, "code_length")),
        None,
    ),
    0.2: (
        ((0.4, "code_length"), (0.6, "code_length")),
        ((0.4, "code_length"), (0.6, "code_length")),
        ((0.2, "code_length"), (0.2, "code_length")),
        None,
    ),
    (0.4, 0.3): (
        ((0.3, "code_length"), (0.7, "code_length")),
        ((0.35, "code_length"), (0.65, "code_length")),
        ((0.4, "code_length"), (0.3, "code_length")),
        None,
    ),
    (1.2, "cm"): (
        ((-0.1, "code_length"), (1.1, "code_length")),
        ((-0.1, "code_length"), (1.1, "code_length")),
        ((1.2, "code_length"), (1.2, "code_length")),
        ("cm", "cm"),
    ),
    ((1.2, "cm"), (2.0, "cm")): (
        ((-0.1, "code_length"), (1.1, "code_length")),
        ((-0.5, "code_length"), (1.5, "code_length")),
        ((1.2, "code_length"), (2.0, "code_length")),
        ("cm", "cm"),
    ),
    ((1.2, "cm"), (0.02, "m")): (
        ((-0.1, "code_length"), (1.1, "code_length")),
        ((-0.5, "code_length"), (1.5, "code_length")),
        ((1.2, "code_length"), (2.0, "code_length")),
        ("cm", "m"),
    ),
}

WEIGHT_FIELDS = (
    None,
    ("gas", "density"),
)

PROJECTION_METHODS = ("integrate", "sum", "min", "max")

BUFF_SIZES = [(800, 800), (1600, 1600), (1254, 1254), (800, 600)]


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
    plot.annotate_streamlines(("gas", f"velocity_{xn}"), ("gas", f"velocity_{yn}"))


CALLBACK_TESTS = (
    ("simple_contour", (simple_contour,)),
    ("simple_velocity", (simple_velocity,)),
    # ("simple_streamlines", (simple_streamlines,)),
    # ("simple_all", (simple_contour, simple_velocity, simple_streamlines)),
)


@requires_ds(M7)
def test_attributes():
    """Test plot member functions that aren't callbacks"""
    plot_field = ("gas", "density")
    decimals = 12

    ds = data_dir_load(M7)
    for ax in "xyz":
        for attr_name in ATTR_ARGS.keys():
            for args in ATTR_ARGS[attr_name]:
                test = PlotWindowAttributeTest(
                    ds, plot_field, ax, attr_name, args, decimals
                )
                test_attributes.__name__ = test.description
                yield test
                for n, r in CALLBACK_TESTS:
                    yield PlotWindowAttributeTest(
                        ds,
                        plot_field,
                        ax,
                        attr_name,
                        args,
                        decimals,
                        callback_id=n,
                        callback_runners=r,
                    )


class TestHideAxesColorbar(unittest.TestCase):

    ds = None

    def setUp(self):
        if self.ds is None:
            self.ds = fake_random_ds(64)
            self.slc = SlicePlot(self.ds, 0, ("gas", "density"))
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
            self.slc = SlicePlot(self.ds, 0, ("gas", "density"))

    def tearDown(self):
        del self.ds
        del self.slc

    def _assert_05cm(self):
        assert_array_equal(
            [self.slc.xlim, self.slc.ylim, self.slc.width],
            [
                (YTQuantity(0.25, "cm"), YTQuantity(0.75, "cm")),
                (YTQuantity(0.25, "cm"), YTQuantity(0.75, "cm")),
                (YTQuantity(0.5, "cm"), YTQuantity(0.5, "cm")),
            ],
        )

    def _assert_05_075cm(self):
        assert_array_equal(
            [self.slc.xlim, self.slc.ylim, self.slc.width],
            [
                (YTQuantity(0.25, "cm"), YTQuantity(0.75, "cm")),
                (YTQuantity(0.125, "cm"), YTQuantity(0.875, "cm")),
                (YTQuantity(0.5, "cm"), YTQuantity(0.75, "cm")),
            ],
        )

    def test_set_width_one(self):
        assert_equal(
            [self.slc.xlim, self.slc.ylim, self.slc.width],
            [(0.0, 1.0), (0.0, 1.0), (1.0, 1.0)],
        )
        assert_true(self.slc._axes_unit_names is None)

    def test_set_width_nonequal(self):
        self.slc.set_width((0.5, 0.8))
        assert_rel_equal(
            [self.slc.xlim, self.slc.ylim, self.slc.width],
            [(0.25, 0.75), (0.1, 0.9), (0.5, 0.8)],
            15,
        )
        assert_true(self.slc._axes_unit_names is None)

    def test_twoargs_eq(self):
        self.slc.set_width(0.5, "cm")
        self._assert_05cm()
        assert_true(self.slc._axes_unit_names == ("cm", "cm"))

    def test_tuple_eq(self):
        self.slc.set_width((0.5, "cm"))
        self._assert_05cm()
        assert_true(self.slc._axes_unit_names == ("cm", "cm"))

    def test_tuple_of_tuples_neq(self):
        self.slc.set_width(((0.5, "cm"), (0.75, "cm")))
        self._assert_05_075cm()
        assert_true(self.slc._axes_unit_names == ("cm", "cm"))


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
            slc = SlicePlot(test_ds, dim, ("gas", "density"))
            for fname in TEST_FLNMS:
                assert_fname(slc.save(fname)[0])

    def test_repr_html(self):
        test_ds = fake_random_ds(16)
        slc = SlicePlot(test_ds, 0, ("gas", "density"))
        slc._repr_html_()

    def test_projection_plot(self):
        test_ds = fake_random_ds(16)
        for dim in range(3):
            proj = ProjectionPlot(test_ds, dim, ("gas", "density"))
            for fname in TEST_FLNMS:
                assert_fname(proj.save(fname)[0])

    def test_projection_plot_ds(self):
        test_ds = fake_random_ds(16)
        reg = test_ds.region([0.5] * 3, [0.4] * 3, [0.6] * 3)
        for dim in range(3):
            proj = ProjectionPlot(test_ds, dim, ("gas", "density"), data_source=reg)
            proj.save()

    def test_projection_plot_c(self):
        test_ds = fake_random_ds(16)
        for center in CENTER_SPECS:
            proj = ProjectionPlot(test_ds, 0, ("gas", "density"), center=center)
            proj.save()

    def test_projection_plot_wf(self):
        test_ds = fake_random_ds(16)
        for wf in WEIGHT_FIELDS:
            proj = ProjectionPlot(test_ds, 0, ("gas", "density"), weight_field=wf)
            proj.save()

    def test_projection_plot_m(self):
        test_ds = fake_random_ds(16)
        for method in PROJECTION_METHODS:
            proj = ProjectionPlot(test_ds, 0, ("gas", "density"), method=method)
            proj.save()

    def test_projection_plot_bs(self):
        test_ds = fake_random_ds(16)
        for bf in BUFF_SIZES:
            proj = ProjectionPlot(test_ds, 0, ("gas", "density"), buff_size=bf)
            image = proj.frb["gas", "density"]

            # note that image.shape is inverted relative to the passed in buff_size
            assert_equal(image.shape[::-1], bf)

    def test_offaxis_slice_plot(self):
        test_ds = fake_random_ds(16)
        slc = OffAxisSlicePlot(test_ds, [1, 1, 1], ("gas", "density"))
        for fname in TEST_FLNMS:
            assert_fname(slc.save(fname)[0])

    def test_offaxis_projection_plot(self):
        test_ds = fake_random_ds(16)
        prj = OffAxisProjectionPlot(test_ds, [1, 1, 1], ("gas", "density"))
        for fname in TEST_FLNMS:
            assert_fname(prj.save(fname)[0])

    def test_creation_with_width(self):
        test_ds = fake_random_ds(16)
        for width in WIDTH_SPECS:
            xlim, ylim, pwidth, aun = WIDTH_SPECS[width]
            plot = ProjectionPlot(test_ds, 0, ("gas", "density"), width=width)

            xlim = [plot.ds.quan(el[0], el[1]) for el in xlim]
            ylim = [plot.ds.quan(el[0], el[1]) for el in ylim]
            pwidth = [plot.ds.quan(el[0], el[1]) for el in pwidth]

            [assert_array_almost_equal(px, x, 14) for px, x in zip(plot.xlim, xlim)]
            [assert_array_almost_equal(py, y, 14) for py, y in zip(plot.ylim, ylim)]
            [assert_array_almost_equal(pw, w, 14) for pw, w in zip(plot.width, pwidth)]
            assert_true(aun == plot._axes_unit_names)


class TestPerFieldConfig(unittest.TestCase):

    ds = None

    def setUp(self):
        from yt.config import ytcfg

        newConfig = {
            ("yt", "default_colormap"): "viridis",
            ("plot", "gas", "log"): False,
            ("plot", "gas", "density", "units"): "lb/yard**3",
            ("plot", "gas", "density", "path_length_units"): "mile",
            ("plot", "gas", "density", "cmap"): "plasma",
            ("plot", "gas", "temperature", "log"): True,
            ("plot", "gas", "temperature", "linthresh"): 100,
            ("plot", "gas", "temperature", "cmap"): "hot",
            ("plot", "gas", "pressure", "log"): True,
            ("plot", "index", "radius", "linthresh"): 1e3,
        }
        # Backup the old config
        oldConfig = {}
        for key in newConfig.keys():
            try:
                val = ytcfg[key]
                oldConfig[key] = val
            except KeyError:
                pass
        for key, val in newConfig.items():
            ytcfg[key] = val

        self.oldConfig = oldConfig
        self.newConfig = newConfig

        fields = [("gas", "density"), ("gas", "temperature"), ("gas", "pressure")]
        units = ["g/cm**3", "K", "dyn/cm**2"]
        fields_to_plot = fields + [("index", "radius")]
        if self.ds is None:
            self.ds = fake_random_ds(16, fields=fields, units=units)
            self.proj = ProjectionPlot(self.ds, 0, fields_to_plot)

    def tearDown(self):
        from yt.config import ytcfg

        del self.ds
        del self.proj
        for key in self.newConfig.keys():
            ytcfg.remove(*key)
        for key, val in self.oldConfig.items():
            ytcfg[key] = val

    def test_units(self):
        from unyt import Unit

        assert_equal(self.proj.frb["gas", "density"].units, Unit("mile*lb/yd**3"))
        assert_equal(self.proj.frb["gas", "temperature"].units, Unit("cm*K"))
        assert_equal(self.proj.frb["gas", "pressure"].units, Unit("dyn/cm"))

    def test_scale(self):

        assert_equal(
            self.proj.plots["gas", "density"].norm_handler.norm_type, Normalize
        )
        assert_equal(
            self.proj.plots["gas", "temperature"].norm_handler.norm_type, SymLogNorm
        )
        assert_allclose_units(
            self.proj.plots["gas", "temperature"].norm_handler.linthresh,
            unyt_array(100, "K*cm"),
        )
        assert_equal(self.proj.plots["gas", "pressure"].norm_handler.norm_type, LogNorm)
        assert_equal(
            self.proj.plots["index", "radius"].norm_handler.norm_type, SymLogNorm
        )

    def test_cmap(self):
        assert_equal(
            self.proj.plots["gas", "density"].colorbar_handler.cmap.name, "plasma"
        )
        assert_equal(
            self.proj.plots["gas", "temperature"].colorbar_handler.cmap.name, "hot"
        )
        assert_equal(
            self.proj.plots["gas", "pressure"].colorbar_handler.cmap.name, "viridis"
        )


def test_on_off_compare():
    # fake density field that varies in the x-direction only
    den = np.arange(32**3) / 32**2 + 1
    den = den.reshape(32, 32, 32)
    den = np.array(den, dtype=np.float64)
    data = dict(density=(den, "g/cm**3"))
    bbox = np.array([[-1.5, 1.5], [-1.5, 1.5], [-1.5, 1.5]])
    ds = load_uniform_grid(data, den.shape, length_unit="Mpc", bbox=bbox, nprocs=64)

    sl_on = SlicePlot(ds, "z", [("gas", "density")])

    L = [0, 0, 1]
    north_vector = [0, 1, 0]
    sl_off = OffAxisSlicePlot(
        ds, L, ("gas", "density"), center=[0, 0, 0], north_vector=north_vector
    )

    assert_array_almost_equal(
        sl_on.frb[("gas", "density")], sl_off.frb[("gas", "density")]
    )

    sl_on.set_buff_size((800, 400))
    sl_on._recreate_frb()
    sl_off.set_buff_size((800, 400))
    sl_off._recreate_frb()

    assert_array_almost_equal(
        sl_on.frb[("gas", "density")], sl_off.frb[("gas", "density")]
    )


def test_plot_particle_field_error():
    ds = fake_random_ds(32, particles=100)

    field_names = [
        ("all", "particle_mass"),
        [("all", "particle_mass"), ("gas", "density")],
        [("gas", "density"), ("all", "particle_mass")],
    ]

    objects_normals = [
        (SlicePlot, 2),
        (SlicePlot, [1, 1, 1]),
        (ProjectionPlot, 2),
        (OffAxisProjectionPlot, [1, 1, 1]),
    ]

    for object, normal in objects_normals:
        for field_name_list in field_names:
            assert_raises(YTInvalidFieldType, object, ds, normal, field_name_list)


def test_setup_origin():
    origin_inputs = (
        "domain",
        "left-window",
        "center-domain",
        "lower-right-window",
        ("window",),
        ("right", "domain"),
        ("lower", "window"),
        ("lower", "right", "window"),
        (0.5, 0.5, "domain"),
        ((50, "cm"), (50, "cm"), "domain"),
    )
    w = (10, "cm")

    ds = fake_random_ds(32, length_unit=100.0)
    generated_limits = []
    # lower limit -> llim
    # upper limit -> ulim
    #                 xllim xulim yllim yulim
    correct_limits = [
        45.0,
        55.0,
        45.0,
        55.0,
        0.0,
        10.0,
        0.0,
        10.0,
        -5.0,
        5.0,
        -5.0,
        5.0,
        -10.0,
        0,
        0,
        10.0,
        0.0,
        10.0,
        0.0,
        10.0,
        -55.0,
        -45.0,
        -55.0,
        -45.0,
        -5.0,
        5.0,
        0.0,
        10.0,
        -10.0,
        0,
        0,
        10.0,
        -5.0,
        5.0,
        -5.0,
        5.0,
        -5.0,
        5.0,
        -5.0,
        5.0,
    ]
    for o in origin_inputs:
        slc = SlicePlot(ds, 2, ("gas", "density"), width=w, origin=o)
        ax = slc.plots[("gas", "density")].axes
        xlims = ax.get_xlim()
        ylims = ax.get_ylim()
        lims = [xlims[0], xlims[1], ylims[0], ylims[1]]
        for l in lims:
            generated_limits.append(l)
    assert_array_almost_equal(correct_limits, generated_limits)


def test_frb_regen():
    ds = fake_random_ds(32)
    slc = SlicePlot(ds, 2, ("gas", "density"))
    slc.set_buff_size(1200)
    assert_equal(slc.frb[("gas", "density")].shape, (1200, 1200))
    slc.set_buff_size((400.0, 200.7))
    assert_equal(slc.frb[("gas", "density")].shape, (200, 400))


def test_set_background_color():
    ds = fake_random_ds(32)
    plot = SlicePlot(ds, 2, ("gas", "density"))
    plot.set_background_color(("gas", "density"), "red")
    plot._setup_plots()
    ax = plot.plots[("gas", "density")].axes
    assert_equal(ax.get_facecolor(), (1.0, 0.0, 0.0, 1.0))


def test_set_unit():
    ds = fake_random_ds(32, fields=(("gas", "temperature"),), units=("K",))
    slc = SlicePlot(ds, 2, ("gas", "temperature"))

    orig_array = slc.frb["gas", "temperature"].copy()

    slc.set_unit(("gas", "temperature"), "degF")

    assert str(slc.frb["gas", "temperature"].units) == "°F"
    assert_array_almost_equal(
        np.array(slc.frb["gas", "temperature"]), np.array(orig_array) * 1.8 - 459.67
    )

    # test that a plot modifying function that destroys the frb preserves the
    # new unit
    slc.set_buff_size(1000)

    assert str(slc.frb["gas", "temperature"].units) == "°F"

    slc.set_buff_size(800)

    slc.set_unit(("gas", "temperature"), "K")
    assert str(slc.frb["gas", "temperature"].units) == "K"
    assert_array_almost_equal(slc.frb["gas", "temperature"], orig_array)

    slc.set_unit(("gas", "temperature"), "keV", equivalency="thermal")
    assert str(slc.frb["gas", "temperature"].units) == "keV"
    assert_array_almost_equal(
        slc.frb["gas", "temperature"], (orig_array * kboltz).to("keV")
    )

    # test that a plot modifying function that destroys the frb preserves the
    # new unit with an equivalency
    slc.set_buff_size(1000)

    assert str(slc.frb["gas", "temperature"].units) == "keV"

    # test that destroying the FRB then changing the unit using an equivalency
    # doesn't error out, see issue #1316
    slc = SlicePlot(ds, 2, ("gas", "temperature"))
    slc.set_buff_size(1000)
    slc.set_unit(("gas", "temperature"), "keV", equivalency="thermal")
    assert str(slc.frb["gas", "temperature"].units) == "keV"


WD = "WDMerger_hdf5_chk_1000/WDMerger_hdf5_chk_1000.hdf5"
blast_wave = "amrvac/bw_2d0000.dat"


@requires_file(WD)
@requires_file(blast_wave)
def test_plot_2d():
    # Cartesian
    ds = fake_random_ds((32, 32, 1), fields=("temperature",), units=("K",))
    slc = SlicePlot(
        ds,
        "z",
        [("gas", "temperature")],
        width=(0.2, "unitary"),
        center=[0.4, 0.3, 0.5],
    )
    slc2 = plot_2d(
        ds, ("gas", "temperature"), width=(0.2, "unitary"), center=[0.4, 0.3]
    )
    slc3 = plot_2d(
        ds,
        ("gas", "temperature"),
        width=(0.2, "unitary"),
        center=ds.arr([0.4, 0.3], "cm"),
    )
    assert_array_equal(
        slc.frb[("gas", "temperature")], slc2.frb[("gas", "temperature")]
    )
    assert_array_equal(
        slc.frb[("gas", "temperature")], slc3.frb[("gas", "temperature")]
    )
    # Cylindrical
    ds = data_dir_load(WD)
    slc = SlicePlot(ds, "theta", [("gas", "density")], width=(30000.0, "km"))
    slc2 = plot_2d(ds, ("gas", "density"), width=(30000.0, "km"))
    assert_array_equal(slc.frb[("gas", "density")], slc2.frb[("gas", "density")])

    # Spherical
    ds = data_dir_load(blast_wave)
    slc = SlicePlot(ds, "phi", [("gas", "density")], width=(1, "unitary"))
    slc2 = plot_2d(ds, ("gas", "density"), width=(1, "unitary"))
    assert_array_equal(slc.frb[("gas", "density")], slc2.frb[("gas", "density")])


def test_symlog_colorbar():
    ds = fake_random_ds(16)

    def _thresh_density(field, data):
        wh = data[("gas", "density")] < 0.5
        ret = data[("gas", "density")]
        ret[wh] = 0
        return ret

    def _neg_density(field, data):
        return -data[("gas", "threshold_density")]

    ds.add_field(
        ("gas", "threshold_density"),
        function=_thresh_density,
        units="g/cm**3",
        sampling_type="cell",
    )
    ds.add_field(
        ("gas", "negative_density"),
        function=_neg_density,
        units="g/cm**3",
        sampling_type="cell",
    )

    for field in [
        ("gas", "density"),
        ("gas", "threshold_density"),
        ("gas", "negative_density"),
    ]:
        plot = SlicePlot(ds, 2, field)
        plot.set_log(field, linthresh=0.1)
        with tempfile.NamedTemporaryFile(suffix="png") as f:
            plot.save(f.name)


def test_symlog_min_zero():
    # see https://github.com/yt-project/yt/issues/3791
    shape = (32, 16, 1)
    a = np.linspace(0, 1, 16)
    b = np.ones((32, 16))
    c = np.reshape(a * b, shape)
    data = {("gas", "density"): c}

    ds = load_uniform_grid(
        data,
        shape,
        bbox=np.array([[0.0, 5.0], [0, 1], [-0.1, +0.1]]),
    )

    p = SlicePlot(ds, "z", "density")
    im_arr = p["gas", "density"].image.get_array()

    # check that no data value was mapped to a NaN (log(0))
    assert np.all(~np.isnan(im_arr))
    # 0 should be mapped to itself since we expect a symlog norm
    assert np.min(im_arr) == 0.0


def test_symlog_extremely_small_vals():
    # check that the plot can be constructed without crashing
    # see https://github.com/yt-project/yt/issues/3858
    # and https://github.com/yt-project/yt/issues/3944
    shape = (64, 64, 1)
    arr = np.full(shape, 5.0e-324)
    arr[0, 0] = -1e12
    arr[1, 1] = 200

    arr2 = np.full(shape, 5.0e-324)
    arr2[0, 0] = -1e12

    arr3 = arr.copy()
    arr3[4, 4] = 0.0

    d = {"scalar_spans_0": arr, "tiny_vmax": arr2, "scalar_tiny_with_0": arr3}

    ds = load_uniform_grid(d, shape)
    for field in d:
        p = SlicePlot(ds, "z", field)
        p["stream", field]


def test_symlog_linthresh_gt_vmax():
    # check that some more edge cases do not crash

    # linthresh will end up being larger than vmax here. This is OK.
    shape = (64, 64, 1)
    arr = np.full(shape, -1e30)
    arr[1, 1] = -1e27
    arr[2, 2] = 1e-12
    arr[3, 3] = 1e-10

    arr2 = -1 * arr.copy()  # also check the reverse
    d = {"linthresh_gt_vmax": arr, "linthresh_lt_vmin": arr2}

    ds = load_uniform_grid(d, shape)
    for field in d:
        p = SlicePlot(ds, "z", field)
        p["stream", field]


def test_symlog_symmetric():
    # should run ok when abs(min negative) == abs(pos max)
    shape = (64, 64, 1)
    arr = np.full(shape, -1e30)
    arr[1, 1] = -1e27
    arr[2, 2] = 1e10
    arr[3, 3] = 1e30
    d = {"linthresh_symmetric": arr}

    ds = load_uniform_grid(d, shape)
    p = SlicePlot(ds, "z", "linthresh_symmetric")
    p["stream", "linthresh_symmetric"]


def test_nan_data():
    data = np.random.random((16, 16, 16)) - 0.5
    data[:9, :9, :9] = np.nan

    data = {"density": data}

    ds = load_uniform_grid(data, [16, 16, 16])

    plot = SlicePlot(ds, "z", ("gas", "density"))

    with tempfile.NamedTemporaryFile(suffix="png") as f:
        plot.save(f.name)


def test_sanitize_valid_normal_vector():
    # note: we don't test against non-cartesian geometries
    # because the way normal "vectors" work isn't cleary
    # specified and works more as an implementation detail
    # at the moment
    ds = fake_amr_ds(geometry="cartesian")

    # We allow maximal polymorphism for axis-aligned directions:
    # even if 3-component vector is received, we want to use the
    # AxisAligned* plotting class (as opposed to OffAxis*) because
    # it's much easier to optimize so it's expected to be more
    # performant.
    axis_label_from_inputs = {
        "x": ["x", 0, [1, 0, 0], [0.1, 0.0, 0.0], [-10, 0, 0]],
        "y": ["y", 1, [0, 1, 0], [0.0, 0.1, 0.0], [0, -10, 0]],
        "z": ["z", 2, [0, 0, 1], [0.0, 0.0, 0.1], [0, 0, -10]],
    }
    for expected, user_inputs in axis_label_from_inputs.items():
        for ui in user_inputs:
            assert NormalPlot.sanitize_normal_vector(ds, ui) == expected

    # arbitrary 3-floats sequences are also valid input.
    # They should be returned as np.ndarrays, but the norm and orientation
    # could be altered. What's important is that their direction is preserved.
    for ui in [(1, 1, 1), [0.0, -3, 1e9], np.ones(3, dtype="int8")]:
        res = NormalPlot.sanitize_normal_vector(ds, ui)
        assert isinstance(res, np.ndarray)
        assert res.dtype == np.float64
        assert_array_equal(
            np.cross(ui, res),
            [0, 0, 0],
        )


def test_reject_invalid_normal_vector():
    ds = fake_amr_ds(geometry="cartesian")
    for ui in [0.0, 1.0, 2.0, 3.0]:
        # acceptable scalar numeric values are restricted to integers.
        # Floats might be a sign that something went wrong upstream
        # e.g., rounding errors, parsing error...
        assert_raises(TypeError, NormalPlot.sanitize_normal_vector, ds, ui)
    for ui in [
        "X",
        "xy",
        "not-an-axis",
        (0, 0, 0),
        [0, 0, 0],
        np.zeros(3),
        [1, 0, 0, 0],
        [1, 0],
        [1],
        [0],
        3,
        10,
    ]:
        assert_raises(ValueError, NormalPlot.sanitize_normal_vector, ds, ui)


def test_dispatch_plot_classes():
    ds = fake_random_ds(16)
    p1 = ProjectionPlot(ds, "z", ("gas", "density"))
    p2 = ProjectionPlot(ds, (1, 2, 3), ("gas", "density"))
    s1 = SlicePlot(ds, "z", ("gas", "density"))
    s2 = SlicePlot(ds, (1, 2, 3), ("gas", "density"))
    assert isinstance(p1, AxisAlignedProjectionPlot)
    assert isinstance(p2, OffAxisProjectionPlot)
    assert isinstance(s1, AxisAlignedSlicePlot)
    assert isinstance(s2, OffAxisSlicePlot)


@requires_module("cartopy")
def test_invalid_swap_projection():
    # projections and transforms will not work
    ds = fake_amr_ds(geometry="geographic")
    slc = SlicePlot(ds, "altitude", ds.field_list[0], origin="native")
    slc.set_mpl_projection("Robinson")
    slc.swap_axes()  # should raise mylog.warning and not toggle _swap_axes
    assert slc._has_swapped_axes is False


def test_set_font():
    # simply check that calling the set_font method doesn't raise an error
    # https://github.com/yt-project/yt/issues/4263
    ds = fake_amr_ds()
    slc = SlicePlot(ds, "x", "Density")
    slc.set_font(
        {
            "family": "sans-serif",
            "style": "italic",
            "weight": "bold",
            "size": 24,
            "color": "blue",
        }
    )
