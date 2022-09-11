import os
import shutil
import tempfile
import unittest

from nose.plugins.attrib import attr

import yt
from yt.testing import ANSWER_TEST_TAG, assert_allclose_units, fake_random_ds
from yt.utilities.answer_testing.framework import PhasePlotAttributeTest

ATTR_ARGS = {
    "annotate_text": [
        (((5e-29, 5e7), "Hello YT"), {}),
        (((5e-29, 5e7), "Hello YT"), {"color": "b"}),
    ],
    "set_title": [((("gas", "mass"), "A phase plot."), {})],
    "set_log": [((("gas", "mass"), False), {})],
    "set_unit": [((("gas", "mass"), "Msun"), {})],
    "set_xlim": [((1e-27, 1e-24), {})],
    "set_ylim": [((1e2, 1e6), {})],
}


@attr(ANSWER_TEST_TAG)
def test_phase_plot_attributes():
    """

    This iterates over the all the plot modification functions in
    ATTR_ARGS. Each time, it compares the images produced by
    PhasePlot to the gold standard.


    """

    x_field = ("gas", "density")
    y_field = ("gas", "temperature")
    z_field = ("gas", "mass")
    decimals = 12
    ds = fake_random_ds(16, fields=("density", "temperature"), units=("g/cm**3", "K"))
    for attr_name in ATTR_ARGS.keys():
        for args in ATTR_ARGS[attr_name]:
            test = PhasePlotAttributeTest(
                ds, x_field, y_field, z_field, attr_name, args, decimals
            )
            test.prefix = f"{attr_name}_{args}"
            test.answer_name = "phase_plot_attributes"
            yield test


def test_set_units():
    fields = ("density", "temperature")
    units = (
        "g/cm**3",
        "K",
    )
    ds = fake_random_ds(16, fields=fields, units=units)
    sp = ds.sphere("max", (1.0, "Mpc"))
    p1 = yt.ProfilePlot(sp, ("index", "radius"), ("gas", "density"))
    p2 = yt.PhasePlot(sp, ("gas", "density"), ("gas", "temperature"), ("gas", "mass"))
    # make sure we can set the units using the tuple without erroring out
    p1.set_unit(("gas", "density"), "Msun/kpc**3")
    p2.set_unit(("gas", "temperature"), "R")


def test_set_labels():
    ds = fake_random_ds(16)
    ad = ds.all_data()
    plot = yt.ProfilePlot(
        ad,
        ("index", "radius"),
        [("gas", "velocity_x"), ("gas", "density")],
        weight_field=None,
    )
    # make sure we can set the labels without erroring out
    plot.set_ylabel("all", "test ylabel")
    plot.set_xlabel("test xlabel")


def test_create_from_dataset():
    ds = fake_random_ds(16)
    plot1 = yt.ProfilePlot(
        ds,
        ("index", "radius"),
        [("gas", "velocity_x"), ("gas", "density")],
        weight_field=None,
    )
    plot2 = yt.ProfilePlot(
        ds.all_data(),
        ("index", "radius"),
        [("gas", "velocity_x"), ("gas", "density")],
        weight_field=None,
    )
    assert_allclose_units(
        plot1.profiles[0][("gas", "density")], plot2.profiles[0][("gas", "density")]
    )
    assert_allclose_units(
        plot1.profiles[0]["velocity_x"], plot2.profiles[0]["velocity_x"]
    )

    plot1 = yt.PhasePlot(ds, ("gas", "density"), ("gas", "velocity_x"), ("gas", "mass"))
    plot2 = yt.PhasePlot(
        ds.all_data(), ("gas", "density"), ("gas", "velocity_x"), ("gas", "mass")
    )
    assert_allclose_units(plot1.profile["mass"], plot2.profile["mass"])


class TestAnnotations(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp()
        cls.curdir = os.getcwd()
        os.chdir(cls.tmpdir)

        ds = fake_random_ds(16)
        ad = ds.all_data()
        cls.fields = [
            ("gas", "velocity_x"),
            ("gas", "velocity_y"),
            ("gas", "velocity_z"),
        ]
        cls.plot = yt.ProfilePlot(
            ad, ("index", "radius"), cls.fields, weight_field=None
        )

    @classmethod
    def tearDownClass(cls):
        os.chdir(cls.curdir)
        shutil.rmtree(cls.tmpdir)

    def test_annotations(self):
        # make sure we can annotate without erroring out
        # annotate the plot with only velocity_x
        self.plot.annotate_title("velocity_x plot", self.fields[0])
        self.plot.annotate_text(1e-1, 1e1, "Annotated velocity_x")

        # annotate the plots with velocity_y and velocity_z with
        # the same annotations
        self.plot.annotate_title("Velocity Plots (Y or Z)", self.fields[1:])
        self.plot.annotate_text(1e-1, 1e1, "Annotated vel_y, vel_z", self.fields[1:])
        self.plot.save()

    def test_annotations_wrong_fields(self):
        from yt.utilities.exceptions import YTFieldNotFound

        with self.assertRaises(YTFieldNotFound):
            self.plot.annotate_title("velocity_x plot", "wrong_field_name")

        with self.assertRaises(YTFieldNotFound):
            self.plot.annotate_text(1e-1, 1e1, "Annotated text", "wrong_field_name")


def test_phaseplot_set_log():
    fields = ("density", "temperature")
    units = (
        "g/cm**3",
        "K",
    )
    ds = fake_random_ds(16, fields=fields, units=units)
    sp = ds.sphere("max", (1.0, "Mpc"))
    p1 = yt.ProfilePlot(sp, ("index", "radius"), ("gas", "density"))
    p2 = yt.PhasePlot(sp, ("gas", "density"), ("gas", "temperature"), ("gas", "mass"))
    # make sure we can set the log-scaling using the tuple without erroring out
    p1.set_log(("gas", "density"), False)
    p2.set_log(("gas", "temperature"), False)
    assert not p1.y_log["gas", "density"]
    assert not p2.y_log

    # make sure we can set the log-scaling using a string without erroring out
    p1.set_log(("gas", "density"), True)
    p2.set_log(("gas", "temperature"), True)
    assert p1.y_log["gas", "density"]
    assert p2.y_log

    # make sure we can set the log-scaling using a field object
    p1.set_log(ds.fields.gas.density, False)
    p2.set_log(ds.fields.gas.temperature, False)
    assert not p1.y_log["gas", "density"]
    assert not p2.y_log


def test_phaseplot_showhide_colorbar_axes():
    fields = ("density", "temperature")
    units = (
        "g/cm**3",
        "K",
    )
    ds = fake_random_ds(16, fields=fields, units=units)
    ad = ds.all_data()
    plot = yt.PhasePlot(ad, ("gas", "density"), ("gas", "temperature"), ("gas", "mass"))

    # make sure we can hide colorbar
    plot.hide_colorbar()
    with tempfile.NamedTemporaryFile(suffix="png") as f1:
        plot.save(f1.name)

    # make sure we can show colorbar
    plot.show_colorbar()
    with tempfile.NamedTemporaryFile(suffix="png") as f2:
        plot.save(f2.name)

    # make sure we can hide axes
    plot.hide_axes()
    with tempfile.NamedTemporaryFile(suffix="png") as f3:
        plot.save(f3.name)

    # make sure we can show axes
    plot.show_axes()
    with tempfile.NamedTemporaryFile(suffix="png") as f4:
        plot.save(f4.name)
