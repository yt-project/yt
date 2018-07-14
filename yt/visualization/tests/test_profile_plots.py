"""
Testsuite for ProfilePlot and PhasePlot



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
import os
import tempfile
import shutil
import unittest
import yt
from yt.data_objects.profiles import create_profile
from yt.testing import \
    fake_random_ds, \
    assert_array_almost_equal, \
    requires_file, \
    assert_fname, \
    assert_allclose_units
from yt.visualization.profile_plotter import \
    ProfilePlot, PhasePlot
from yt.visualization.tests.test_plotwindow import \
    TEST_FLNMS
from yt.utilities.answer_testing.framework import \
    PhasePlotAttributeTest, \
    requires_ds, \
    data_dir_load

ATTR_ARGS = {"annotate_text": [(((5e-29, 5e7), "Hello YT"), {}), 
                               (((5e-29, 5e7), "Hello YT"), {'color':'b'})],
             
             "set_title": [(('cell_mass', 'A phase plot.'), {})],
             "set_log": [(('cell_mass', False), {})],
             "set_unit": [(('cell_mass', 'Msun'), {})],
             "set_xlim": [((1e-27, 1e-24), {})],
             "set_ylim": [((1e2, 1e6), {})]}


g30 = "IsolatedGalaxy/galaxy0030/galaxy0030"
ETC46 = "enzo_tiny_cosmology/DD0046/DD0046"

@requires_ds(g30, big_data=True)
def test_phase_plot_attributes():
    '''

    This iterates over the all the plot modification functions in 
    ATTR_ARGS. Each time, it compares the images produced by 
    PhasePlot to the gold standard.
    

    '''

    x_field = 'density'
    y_field = 'temperature'
    z_field = 'cell_mass'
    decimals = 12
    ds = data_dir_load(g30)
    for attr_name in ATTR_ARGS.keys():
        for args in ATTR_ARGS[attr_name]:
            test = PhasePlotAttributeTest(ds, x_field, y_field, z_field, 
                                          attr_name, args, decimals)
            test_phase_plot_attributes.__name__ = test.description
            yield test

class TestProfilePlotSave(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.curdir = os.getcwd()
        os.chdir(self.tmpdir)

    def tearDown(self):
        os.chdir(self.curdir)
        shutil.rmtree(self.tmpdir)

    def test_profile_plot(self):
        fields = ('density', 'temperature', 'velocity_x', 'velocity_y',
                  'velocity_z')
        units = ('g/cm**3', 'K', 'cm/s', 'cm/s', 'cm/s')
        test_ds = fake_random_ds(16, fields=fields, units=units)
        regions = [test_ds.region([0.5]*3, [0.4]*3, [0.6]*3), test_ds.all_data()]
        pr_fields = [('density', 'temperature'), ('density', 'velocity_x'),
                     ('temperature', 'cell_mass'), ('density', 'radius'),
                     ('velocity_magnitude', 'cell_mass')]
        profiles = []
        for reg in regions:
            for x_field, y_field in pr_fields:
                profiles.append(ProfilePlot(reg, x_field, y_field))
                profiles.append(ProfilePlot(reg, x_field, y_field,
                                            fractional=True, accumulation=True))
                p1d = create_profile(reg, x_field, y_field)
                profiles.append(ProfilePlot.from_profiles(p1d))
        p1 = create_profile(test_ds.all_data(), 'density', 'temperature')
        p2 = create_profile(test_ds.all_data(), 'density', 'velocity_x')
        profiles.append(ProfilePlot.from_profiles(
            [p1, p2], labels=['temperature', 'velocity']))
        profiles[0]._repr_html_()
        for p in profiles:
            for fname in TEST_FLNMS:
                assert_fname(p.save(fname)[0])

    def test_phase_plot(self):
        fields = ('density', 'temperature', 'velocity_x', 'velocity_y',
                  'velocity_z')
        units = ('g/cm**3', 'K', 'cm/s', 'cm/s', 'cm/s')
        test_ds = fake_random_ds(16, fields=fields, units=units)
        regions = [test_ds.region([0.5]*3, [0.4]*3, [0.6]*3), test_ds.all_data()]
        phases = []
        ph_fields = [('density', 'temperature', 'cell_mass'),
                     ('density', 'velocity_x', 'cell_mass'),
                     ('radius', 'temperature', 'velocity_magnitude')]
        for reg in regions:
            for x_field, y_field, z_field in ph_fields:
                # set n_bins to [16, 16] since matplotlib's postscript
                # renderer is slow when it has to write a lot of polygons
                phases.append(PhasePlot(reg, x_field, y_field, z_field,
                                        x_bins=16, y_bins=16))
                phases.append(PhasePlot(reg, x_field, y_field, z_field,
                                        fractional=True, accumulation=True,
                                        x_bins=16, y_bins=16))
                p2d = create_profile(reg, [x_field, y_field], z_field,
                                     n_bins=[16, 16])
                phases.append(PhasePlot.from_profile(p2d))
        pp = PhasePlot(test_ds.all_data(), 'density', 'temperature', 'cell_mass')
        pp.set_xlim(0.3, 0.8)
        pp.set_ylim(0.4, 0.6)
        pp._setup_plots()
        xlim = pp.plots['cell_mass'].axes.get_xlim()
        ylim = pp.plots['cell_mass'].axes.get_ylim()
        assert_array_almost_equal(xlim, (0.3, 0.8))
        assert_array_almost_equal(ylim, (0.4, 0.6))
        phases.append(pp)
        phases[0]._repr_html_()
        for p in phases:
            for fname in TEST_FLNMS:
                assert_fname(p.save(fname)[0])

    @requires_file(ETC46)
    def test_profile_plot_multiple_field_multiple_plot(self):
        ds = yt.load(ETC46)
        sphere = ds.sphere("max", (1.0, "Mpc"))
        profiles = []
        profiles.append(yt.create_profile(
            sphere, ["radius"],
            fields=["density"], n_bins=32))
        profiles.append(yt.create_profile(
            sphere, ["radius"],
            fields=["density"], n_bins=64))
        profiles.append(yt.create_profile(
            sphere, ["radius"],
            fields=["dark_matter_density"], n_bins=64))

        plot = yt.ProfilePlot.from_profiles(profiles)
        assert_fname(plot.save()[0])

@requires_file(ETC46)
def test_set_units():
    ds = yt.load(ETC46)
    sp = ds.sphere("max", (1.0, "Mpc"))
    p1 = yt.ProfilePlot(sp, "radius", ("enzo", "Density"))
    p2 = yt.PhasePlot(sp, ("enzo", "Density"), ("enzo", "Temperature"), "cell_mass")
    # make sure we can set the units using the tuple without erroring out
    p1.set_unit(("enzo", "Density"), "Msun/kpc**3")
    p2.set_unit(("enzo", "Temperature"), "R")

def test_set_labels():
    ds = fake_random_ds(16)
    ad = ds.all_data()
    plot = yt.ProfilePlot(ad, "radius", ["velocity_x", "density"], weight_field=None)
    # make sure we can set the labels without erroring out
    plot.set_ylabel("all", "test ylabel")
    plot.set_xlabel("test xlabel")

def test_create_from_dataset():
    ds = fake_random_ds(16)
    plot1 = yt.ProfilePlot(ds, "radius", ["velocity_x", "density"],
                           weight_field=None)
    plot2 = yt.ProfilePlot(ds.all_data(), "radius", ["velocity_x", "density"],
                           weight_field=None)
    assert_allclose_units(
        plot1.profiles[0]['density'], plot2.profiles[0]['density'])
    assert_allclose_units(
        plot1.profiles[0]['velocity_x'], plot2.profiles[0]['velocity_x'])

    plot1 = yt.PhasePlot(ds, 'density', 'velocity_x', 'cell_mass')
    plot2 = yt.PhasePlot(ds.all_data(), 'density', 'velocity_x', 'cell_mass')
    assert_allclose_units(
        plot1.profile['cell_mass'], plot2.profile['cell_mass'])

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
            self.plot.annotate_title("velocity_x plot",  "wrong_field_name")

        with self.assertRaises(YTFieldNotFound):
            self.plot.annotate_text(1e-1, 1e1, "Annotated text", "wrong_field_name")

@requires_file(ETC46)
def test_phaseplot_set_log():
    ds = yt.load(ETC46)
    sp = ds.sphere("max", (1.0, "Mpc"))
    p1 = yt.ProfilePlot(sp, "radius", ("enzo", "Density"))
    p2 = yt.PhasePlot(sp, ("enzo", "Density"), ("enzo", "Temperature"), "cell_mass")
    # make sure we can set the log-scaling using the tuple without erroring out
    p1.set_log(("enzo", "Density"), False)
    p2.set_log(("enzo", "Temperature"), False)
    assert p1.y_log["enzo", "Density"] is False
    assert p2.y_log is False

    # make sure we can set the log-scaling using a string without erroring out
    p1.set_log("Density", True)
    p2.set_log("Temperature", True)
    assert p1.y_log["enzo", "Density"] is True
    assert p2.y_log is True

    # make sure we can set the log-scaling using a field object
    p1.set_log(ds.fields.enzo.Density, False)
    p2.set_log(ds.fields.enzo.Temperature, False)
    assert p1.y_log["enzo", "Density"] is False
    assert p2.y_log is False
