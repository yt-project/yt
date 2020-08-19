"""
Test suite for Particle Plots



"""

# -----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------
import os
import shutil
import tempfile
import unittest

import numpy as np
import pytest

from yt.convenience import load
from yt.data_objects.particle_filters import add_particle_filter
from yt.data_objects.profiles import create_profile
from yt.testing import (
    assert_allclose,
    assert_array_almost_equal,
    assert_fname,
    fake_particle_ds,
    requires_file,
)
from yt.units.yt_array import YTArray
from yt.utilities.answer_testing import utils
from yt.utilities.answer_testing.answer_tests import (
    phase_plot_attribute_test,
    plot_window_attribute_test,
)
from yt.visualization.api import ParticlePhasePlot, ParticlePlot, ParticleProjectionPlot
from yt.visualization.tests.test_plotwindow import WIDTH_SPECS

TEST_FLNMS = [None, "test", "test.png", "test.eps", "test.ps", "test.pdf"]

CENTER_SPECS = (
    "c",
    "C",
    "center",
    "Center",
    [0.5, 0.5, 0.5],
    [[0.2, 0.3, 0.4], "cm"],
    YTArray([0.3, 0.4, 0.7], "cm"),
)

WEIGHT_FIELDS = (
    None,
    "particle_ones",
    ("all", "particle_mass"),
)

PHASE_FIELDS = [
    ("particle_velocity_x", "particle_position_z", "particle_mass"),
    ("particle_position_x", "particle_position_y", "particle_ones"),
    ("particle_velocity_x", "particle_velocity_y", ["particle_mass", "particle_ones"]),
]

g30 = "IsolatedGalaxy/galaxy0030/galaxy0030"


@pytest.mark.answer_test
@pytest.mark.big_data
class TestParticlePlotAnswer:
    @pytest.mark.usefixtures("hashing")
    @utils.requires_ds(g30)
    def test_particle_projection_answers(self, axis, attr_name, attr_args):
        """
        This iterates over the all the plot modification functions in
        PROJ_ATTR_ARGS. Each time, it compares the images produced by
        ParticleProjectionPlot to the gold standard.
        """
        plot_field = "particle_mass"
        ds = utils.data_dir_load(g30)
        pw = plot_window_attribute_test(
            ds, plot_field, axis, attr_name, attr_args, "ParticleProjectionPlot"
        )
        self.hashes.update({"plot_window_attribute": pw})

    @pytest.mark.usefixtures("hashing")
    @utils.requires_ds(g30)
    def test_particle_projection_filter(self, axis, attr_args):
        """
        This tests particle projection plots for filter fields.
        """

        def formed_star(pfilter, data):
            filter = data["all", "creation_time"] > 0
            return filter

        add_particle_filter(
            "formed_star",
            function=formed_star,
            filtered_type="all",
            requires=["creation_time"],
        )
        plot_field = ("formed_star", "particle_mass")
        ds = utils.data_dir_load(g30)
        ds.add_particle_filter("formed_star")
        pw = plot_window_attribute_test(
            ds, plot_field, axis, "set_log", attr_args, "ParticleProjectionPlot"
        )
        self.hashes.update({"plot_window_attribute": pw})

    @pytest.mark.usefixtures("hashing")
    @utils.requires_ds(g30)
    def test_particle_phase_answers(self, attr_name, attr_args):
        """
        This iterates over the all the plot modification functions in
        PHASE_ATTR_ARGS. Each time, it compares the images produced by
        ParticlePhasePlot to the gold standard.
        """
        ds = utils.data_dir_load(g30)
        x_field = "particle_velocity_x"
        y_field = "particle_velocity_y"
        z_field = "particle_mass"
        pp = phase_plot_attribute_test(
            ds, x_field, y_field, z_field, attr_name, attr_args, "ParticlePhasePlot"
        )
        self.hashes.update({"phase_plot_attribute": pp})


class TestParticlePhasePlotSave(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.curdir = os.getcwd()
        os.chdir(self.tmpdir)

    def tearDown(self):
        os.chdir(self.curdir)
        shutil.rmtree(self.tmpdir)

    def test_particle_phase_plot(self):
        test_ds = fake_particle_ds()
        data_sources = [
            test_ds.region([0.5] * 3, [0.4] * 3, [0.6] * 3),
            test_ds.all_data(),
        ]
        particle_phases = []

        for source in data_sources:
            for x_field, y_field, z_fields in PHASE_FIELDS:
                particle_phases.append(
                    ParticlePhasePlot(
                        source, x_field, y_field, z_fields, x_bins=16, y_bins=16
                    )
                )

                particle_phases.append(
                    ParticlePhasePlot(
                        source,
                        x_field,
                        y_field,
                        z_fields,
                        x_bins=16,
                        y_bins=16,
                        deposition="cic",
                    )
                )

                pp = create_profile(
                    source,
                    [x_field, y_field],
                    z_fields,
                    weight_field="particle_ones",
                    n_bins=[16, 16],
                )

                particle_phases.append(ParticlePhasePlot.from_profile(pp))
        particle_phases[0]._repr_html_()
        for p in particle_phases:
            for fname in TEST_FLNMS:
                assert_fname(p.save(fname)[0])


tgal = "TipsyGalaxy/galaxy.00300"


@requires_file(tgal)
def test_particle_phase_plot_semantics():
    ds = load(tgal)
    ad = ds.all_data()
    dens_ex = ad.quantities.extrema(("Gas", "density"))
    temp_ex = ad.quantities.extrema(("Gas", "temperature"))
    plot = ParticlePlot(
        ds, ("Gas", "density"), ("Gas", "temperature"), ("Gas", "particle_mass")
    )
    plot.set_log(("Gas", "density"), True)
    plot.set_log(("Gas", "temperature"), True)
    p = plot.profile

    # bin extrema are field extrema
    assert dens_ex[0] - np.spacing(dens_ex[0]) == p.x_bins[0]
    assert dens_ex[-1] + np.spacing(dens_ex[-1]) == p.x_bins[-1]
    assert temp_ex[0] - np.spacing(temp_ex[0]) == p.y_bins[0]
    assert temp_ex[-1] + np.spacing(temp_ex[-1]) == p.y_bins[-1]

    # bins are evenly spaced in log space
    logxbins = np.log10(p.x_bins)
    dxlogxbins = logxbins[1:] - logxbins[:-1]
    assert_allclose(dxlogxbins, dxlogxbins[0])

    logybins = np.log10(p.y_bins)
    dylogybins = logybins[1:] - logybins[:-1]
    assert_allclose(dylogybins, dylogybins[0])

    plot.set_log(("Gas", "density"), False)
    plot.set_log(("Gas", "temperature"), False)
    p = plot.profile

    # bin extrema are field extrema
    assert dens_ex[0] - np.spacing(dens_ex[0]) == p.x_bins[0]
    assert dens_ex[-1] + np.spacing(dens_ex[-1]) == p.x_bins[-1]
    assert temp_ex[0] - np.spacing(temp_ex[0]) == p.y_bins[0]
    assert temp_ex[-1] + np.spacing(temp_ex[-1]) == p.y_bins[-1]

    # bins are evenly spaced in log space
    dxbins = p.x_bins[1:] - p.x_bins[:-1]
    assert_allclose(dxbins, dxbins[0])

    dybins = p.y_bins[1:] - p.y_bins[:-1]
    assert_allclose(dybins, dybins[0])


@requires_file(tgal)
def test_set_units():
    ds = load(tgal)
    sp = ds.sphere("max", (1.0, "Mpc"))
    pp = ParticlePhasePlot(
        sp, ("Gas", "density"), ("Gas", "temperature"), ("Gas", "particle_mass")
    )
    # make sure we can set the units using the tuple without erroring out
    pp.set_unit(("Gas", "particle_mass"), "Msun")


@requires_file(tgal)
def test_switch_ds():
    """
    Tests the _switch_ds() method for ParticleProjectionPlots that as of
    25th October 2017 requires a specific hack in plot_container.py
    """
    ds = load(tgal)
    ds2 = load(tgal)

    plot = ParticlePlot(
        ds,
        ("Gas", "particle_position_x"),
        ("Gas", "particle_position_y"),
        ("Gas", "density"),
    )

    plot._switch_ds(ds2)

    return


class TestParticleProjectionPlotSave(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.curdir = os.getcwd()
        os.chdir(self.tmpdir)

    def tearDown(self):
        os.chdir(self.curdir)
        shutil.rmtree(self.tmpdir)

    def test_particle_plot(self):
        test_ds = fake_particle_ds()
        for dim in range(3):
            pplot = ParticleProjectionPlot(test_ds, dim, "particle_mass")
            for fname in TEST_FLNMS:
                assert_fname(pplot.save(fname)[0])

    def test_particle_plot_ds(self):
        test_ds = fake_particle_ds()
        ds_region = test_ds.region([0.5] * 3, [0.4] * 3, [0.6] * 3)
        for dim in range(3):
            pplot_ds = ParticleProjectionPlot(
                test_ds, dim, "particle_mass", data_source=ds_region
            )
            pplot_ds.save()

    def test_particle_plot_c(self):
        test_ds = fake_particle_ds()
        for center in CENTER_SPECS:
            for dim in range(3):
                pplot_c = ParticleProjectionPlot(
                    test_ds, dim, "particle_mass", center=center
                )
                pplot_c.save()

    def test_particle_plot_wf(self):
        test_ds = fake_particle_ds()
        for dim in range(3):
            for weight_field in WEIGHT_FIELDS:
                pplot_wf = ParticleProjectionPlot(
                    test_ds, dim, "particle_mass", weight_field=weight_field
                )
                pplot_wf.save()

    def test_creation_with_width(self):
        test_ds = fake_particle_ds()
        for width, (xlim, ylim, pwidth, aun) in WIDTH_SPECS.items():
            plot = ParticleProjectionPlot(test_ds, 0, "particle_mass", width=width)

            xlim = [plot.ds.quan(el[0], el[1]) for el in xlim]
            ylim = [plot.ds.quan(el[0], el[1]) for el in ylim]
            pwidth = [plot.ds.quan(el[0], el[1]) for el in pwidth]

            [assert_array_almost_equal(px, x, 14) for px, x in zip(plot.xlim, xlim)]
            [assert_array_almost_equal(py, y, 14) for py, y in zip(plot.ylim, ylim)]
            [assert_array_almost_equal(pw, w, 14) for pw, w in zip(plot.width, pwidth)]


def test_particle_plot_instance():
    """
    Tests the type of plot instance returned by ParticlePlot.

    If x_field and y_field are any combination of valid particle_position in x,
    y or z axis,then ParticleProjectionPlot instance is expected.


    """
    ds = fake_particle_ds()
    x_field = ("all", "particle_position_x")
    y_field = ("all", "particle_position_y")
    z_field = ("all", "particle_velocity_x")

    plot = ParticlePlot(ds, x_field, y_field)
    assert isinstance(plot, ParticleProjectionPlot)

    plot = ParticlePlot(ds, y_field, x_field)
    assert isinstance(plot, ParticleProjectionPlot)

    plot = ParticlePlot(ds, x_field, z_field)
    assert isinstance(plot, ParticlePhasePlot)
