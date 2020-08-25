import os
import shutil
import tempfile
import unittest

import pytest

import yt
from yt.data_objects.profiles import create_profile
from yt.testing import assert_allclose_units, assert_array_almost_equal, fake_random_ds
from yt.utilities.answer_testing.answer_tests import generic_image, phase_plot_attribute
from yt.visualization.profile_plotter import PhasePlot, ProfilePlot


def compare(plot):
    def image_from_plot():
        tmpfd, tmpfname = tempfile.mkstemp(suffix=".png")
        os.close(tmpfd)
        return plot.save(tmpfname)

    gi = generic_image(image_from_plot)
    return gi


@pytest.mark.answer_test
@pytest.mark.usefixtures("temp_dir", "hashing")
class TestProfilePlots:
    answer_file = None

    def test_phase_plot_attributes(self, ds_random, attr_name, args):
        """
        This iterates over the all the plot modification functions in
        ATTR_ARGS. Each time, it compares the images produced by
        PhasePlot to the gold standard.
        """
        x_field = "density"
        y_field = "temperature"
        z_field = "cell_mass"
        ppat = phase_plot_attribute(
            ds_random, x_field, y_field, z_field, attr_name, args
        )
        self.hashes.update({"phase_plot_attribute": ppat})

    def test_profile_plot(self, ds_test, region, x_field, y_field):
        profiles = []
        profiles.append(ProfilePlot(region, x_field, y_field))
        profiles.append(
            ProfilePlot(region, x_field, y_field, fractional=True, accumulation=True)
        )
        p1d = create_profile(region, x_field, y_field)
        profiles.append(ProfilePlot.from_profiles(p1d))
        p1 = create_profile(ds_test.all_data(), "density", "temperature")
        p2 = create_profile(ds_test.all_data(), "density", "velocity_x")
        profiles.append(
            ProfilePlot.from_profiles([p1, p2], labels=["temperature", "velocity"])
        )
        profiles[0]._repr_html_()
        for idx, plot in enumerate(profiles):
            gi = compare(plot)
            if "generic_image" not in self.hashes:
                self.hashes.update({"generic_image": {str(idx): gi}})
            else:
                self.hashes["generic_image"].update({str(idx): gi})

    def test_phase_plot(self, ds_test, region, x_field, y_field, z_field):
        phases = []
        # set n_bins to [16, 16] since matplotlib's postscript
        # renderer is slow when it has to write a lot of polygons
        phases.append(
            PhasePlot(region, x_field, y_field, z_field, x_bins=16, y_bins=16)
        )
        phases.append(
            PhasePlot(
                region,
                x_field,
                y_field,
                z_field,
                fractional=True,
                accumulation=True,
                x_bins=16,
                y_bins=16,
            )
        )
        p2d = create_profile(region, [x_field, y_field], z_field, n_bins=[16, 16])
        phases.append(PhasePlot.from_profile(p2d))
        pp = PhasePlot(ds_test.all_data(), "density", "temperature", "cell_mass")
        pp.set_xlim(0.3, 0.8)
        pp.set_ylim(0.4, 0.6)
        pp._setup_plots()
        xlim = pp.plots["cell_mass"].axes.get_xlim()
        ylim = pp.plots["cell_mass"].axes.get_ylim()
        assert_array_almost_equal(xlim, (0.3, 0.8))
        assert_array_almost_equal(ylim, (0.4, 0.6))
        phases.append(pp)
        phases[0]._repr_html_()
        for idx, plot in enumerate(phases):
            gi = compare(plot)
            if "generic_image" not in self.hashes:
                self.hashes.update({"generic_image": {str(idx): gi}})
            else:
                self.hashes["generic_image"].update({str(idx): gi})

    def test_profile_plot_multiple_field_multiple_plot(self, ds_mult):
        sphere = ds_mult.sphere("max", (1.0, "Mpc"))
        profiles = []
        profiles.append(
            yt.create_profile(sphere, ["radius"], fields=["density"], n_bins=32)
        )
        profiles.append(
            yt.create_profile(sphere, ["radius"], fields=["density"], n_bins=64)
        )
        profiles.append(
            yt.create_profile(
                sphere, ["radius"], fields=["dark_matter_density"], n_bins=64
            )
        )
        plot = yt.ProfilePlot.from_profiles(profiles)
        gi = compare(plot)
        self.hashes.update({"generic_image": gi})


def test_set_units():
    fields = ("density", "temperature")
    units = (
        "g/cm**3",
        "K",
    )
    ds = fake_random_ds(16, fields=fields, units=units)
    sp = ds.sphere("max", (1.0, "Mpc"))
    p1 = yt.ProfilePlot(sp, "radius", ("gas", "density"))
    p2 = yt.PhasePlot(sp, ("gas", "density"), ("gas", "temperature"), "cell_mass")
    # make sure we can set the units using the tuple without erroring out
    p1.set_unit(("gas", "density"), "Msun/kpc**3")
    p2.set_unit(("gas", "temperature"), "R")


def test_set_labels():
    ds = fake_random_ds(16)
    ad = ds.all_data()
    plot = yt.ProfilePlot(ad, "radius", ["velocity_x", "density"], weight_field=None)
    # make sure we can set the labels without erroring out
    plot.set_ylabel("all", "test ylabel")
    plot.set_xlabel("test xlabel")


def test_create_from_dataset():
    ds = fake_random_ds(16)
    plot1 = yt.ProfilePlot(ds, "radius", ["velocity_x", "density"], weight_field=None)
    plot2 = yt.ProfilePlot(
        ds.all_data(), "radius", ["velocity_x", "density"], weight_field=None
    )
    assert_allclose_units(plot1.profiles[0]["density"], plot2.profiles[0]["density"])
    assert_allclose_units(
        plot1.profiles[0]["velocity_x"], plot2.profiles[0]["velocity_x"]
    )

    plot1 = yt.PhasePlot(ds, "density", "velocity_x", "cell_mass")
    plot2 = yt.PhasePlot(ds.all_data(), "density", "velocity_x", "cell_mass")
    assert_allclose_units(plot1.profile["cell_mass"], plot2.profile["cell_mass"])


class TestAnnotations(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp()
        cls.curdir = os.getcwd()
        os.chdir(cls.tmpdir)

        ds = fake_random_ds(16)
        ad = ds.all_data()
        cls.fields = ["velocity_x", "velocity_y", "velocity_z"]
        cls.plot = yt.ProfilePlot(ad, "radius", cls.fields, weight_field=None)

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
    p1 = yt.ProfilePlot(sp, "radius", ("gas", "density"))
    p2 = yt.PhasePlot(sp, ("gas", "density"), ("gas", "temperature"), "cell_mass")
    # make sure we can set the log-scaling using the tuple without erroring out
    p1.set_log(("gas", "density"), False)
    p2.set_log(("gas", "temperature"), False)
    assert p1.y_log["gas", "density"] is False
    assert p2.y_log is False

    # make sure we can set the log-scaling using a string without erroring out
    p1.set_log("density", True)
    p2.set_log("temperature", True)
    assert p1.y_log["gas", "density"] is True
    assert p2.y_log is True

    # make sure we can set the log-scaling using a field object
    p1.set_log(ds.fields.gas.density, False)
    p2.set_log(ds.fields.gas.temperature, False)
    assert p1.y_log["gas", "density"] is False
    assert p2.y_log is False


def test_phaseplot_showhide_colorbar_axes():
    fields = ("density", "temperature")
    units = (
        "g/cm**3",
        "K",
    )
    ds = fake_random_ds(16, fields=fields, units=units)
    ad = ds.all_data()
    plot = yt.PhasePlot(ad, ("gas", "density"), ("gas", "temperature"), "cell_mass")

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
