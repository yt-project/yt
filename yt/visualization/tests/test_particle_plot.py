"""
Test suite for Particle Plots



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

import numpy as np
from nose.plugins.attrib import attr

from yt.data_objects.profiles import create_profile
from yt.visualization.tests.test_plotwindow import \
    WIDTH_SPECS, ATTR_ARGS
from yt.data_objects.particle_filters import add_particle_filter
from yt.testing import \
    fake_particle_ds, \
    assert_array_almost_equal, \
    assert_allclose, \
    fake_random_ds, ANSWER_TEST_TAG
from yt.utilities.answer_testing.framework import \
    PlotWindowAttributeTest, \
    PhasePlotAttributeTest, GenericImageTest
from yt.visualization.api import \
    ParticlePlot, \
    ParticleProjectionPlot, \
    ParticlePhasePlot
from yt.units.yt_array import YTArray

def setup():
    """Test specific setup."""
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"

#  override some of the plotwindow ATTR_ARGS
PROJ_ATTR_ARGS = ATTR_ARGS.copy()
PROJ_ATTR_ARGS["set_cmap"] = [(('particle_mass', 'RdBu'), {}),
                                  (('particle_mass', 'kamae'), {})]
PROJ_ATTR_ARGS["set_log"] = [(('particle_mass', False), {})]
PROJ_ATTR_ARGS["set_zlim"] = [(('particle_mass', 1e-25, 1e-23), {}),
                                  (('particle_mass', 1e-25, None),
                                   {'dynamic_range': 4})]
PROJ_ATTR_ARGS["set_buff_size"] = [((1600, ), {}), ((600, ), {})]
PHASE_ATTR_ARGS = {"annotate_text": [(((5e-29, 5e7), "Hello YT"), {}),
                               (((5e-29, 5e7), "Hello YT"), {'color':'b'})],
                   "set_title": [(('particle_mass', 'A phase plot.'), {})],
                   "set_log": [(('particle_mass', False), {})],
                   "set_unit": [(('particle_mass', 'Msun'), {})],
                   "set_xlim": [((-4e7, 4e7), {})],
                   "set_ylim": [((-4e7, 4e7), {})]}

CENTER_SPECS = (
    "c",
    "C",
    "center",
    "Center",
    [0.5, 0.5, 0.5],
    [[0.2, 0.3, 0.4], "cm"],
    YTArray([0.3, 0.4, 0.7], "cm")
)

WEIGHT_FIELDS = (
    None,
    'particle_ones',
    ('all', 'particle_mass'),
)

PHASE_FIELDS = [('particle_velocity_x', 'particle_position_z', 'particle_mass'),
                ('particle_position_x', 'particle_position_y', 'particle_ones'),
                ('particle_velocity_x', 'particle_velocity_y',
                 ['particle_mass', 'particle_ones'])]

def compare(ds, plot, test_prefix, test_name, decimals=12):
    def image_from_plot(filename_prefix):
        return plot.save(filename_prefix)

    image_from_plot.__name__ = "particle_{}".format(test_prefix)
    test = GenericImageTest(ds, image_from_plot, decimals)
    test.prefix = test_prefix
    test.answer_name = test_name
    return test

@attr(ANSWER_TEST_TAG)
def test_particle_projection_answers():
    '''

    This iterates over the all the plot modification functions in
    PROJ_ATTR_ARGS. Each time, it compares the images produced by
    ParticleProjectionPlot to the gold standard.


    '''

    plot_field = 'particle_mass'
    decimals = 12
    ds = fake_particle_ds()
    for ax in 'xyz':
        for attr_name in PROJ_ATTR_ARGS.keys():
            for args in PROJ_ATTR_ARGS[attr_name]:
                test = PlotWindowAttributeTest(ds, plot_field, ax,
                                               attr_name,
                                               args, decimals,
                                               'ParticleProjectionPlot')
                test.prefix = "%s_%s_%s" % (ax, attr_name, plot_field)
                test.answer_name = "particle_projection_attributes"
                yield test


@attr(ANSWER_TEST_TAG)
def test_particle_projection_filter():
    '''

    This tests particle projection plots for filter fields.


    '''

    def star(pfilter, data):
        filter = data["all", "particle_velocity_x"] > data["all", "particle_velocity_y"]
        return filter

    add_particle_filter("star", function=star, filtered_type='all',
                        requires=["particle_velocity_x", "particle_velocity_y"])

    plot_field = ('star', 'particle_mass')

    decimals = 12
    ds = fake_particle_ds()
    ds.add_particle_filter('star')
    for ax in 'xyz':
        attr_name = "set_log"
        for args in PROJ_ATTR_ARGS[attr_name]:
            test = PlotWindowAttributeTest(ds, plot_field, ax,
                                           attr_name,
                                           args, decimals,
                                           'ParticleProjectionPlot')
            test.prefix = "%s_%s_%s" % (ax, attr_name, plot_field)
            test.answer_name = "particle_projection_filter"
            yield test


@attr(ANSWER_TEST_TAG)
def test_particle_phase_answers():
    '''

    This iterates over the all the plot modification functions in
    PHASE_ATTR_ARGS. Each time, it compares the images produced by
    ParticlePhasePlot to the gold standard.

    '''

    decimals = 12
    ds = fake_particle_ds()

    x_field = 'particle_velocity_x'
    y_field = 'particle_velocity_y'
    z_field = 'particle_mass'
    for attr_name in PHASE_ATTR_ARGS.keys():
        for args in PHASE_ATTR_ARGS[attr_name]:
            test = PhasePlotAttributeTest(ds, x_field, y_field, z_field,
                                          attr_name, args, decimals,
                                          'ParticlePhasePlot')
            test.prefix = "%s_%s" % (attr_name, attr_name)
            test.answer_name = "particle_phase_attributes"
            yield test

@attr(ANSWER_TEST_TAG)
def test_particle_phase_plot():
    test_ds = fake_particle_ds()
    data_sources = [test_ds.region([0.5] * 3, [0.4] * 3, [0.6] * 3),
                    test_ds.all_data()]
    particle_phases = []

    for source in data_sources:
        for x_field, y_field, z_fields in PHASE_FIELDS:
            particle_phases.append(ParticlePhasePlot(source,
                                                     x_field,
                                                     y_field,
                                                     z_fields,
                                                     x_bins=16,
                                                     y_bins=16))

            particle_phases.append(ParticlePhasePlot(source,
                                                     x_field,
                                                     y_field,
                                                     z_fields,
                                                     x_bins=16,
                                                     y_bins=16,
                                                     deposition='cic'))

            pp = create_profile(source, [x_field, y_field], z_fields,
                                weight_field='particle_ones',
                                n_bins=[16, 16])

            particle_phases.append(ParticlePhasePlot.from_profile(pp))
    particle_phases[0]._repr_html_()
    for idx, plot in enumerate(particle_phases):
        test_prefix = "%s_%s" % (plot.plots.keys(), idx)
        yield compare(test_ds, plot, test_prefix=test_prefix,
                      test_name="particle_phase_plots")

def test_particle_phase_plot_semantics():
    ds = fake_particle_ds()
    ad = ds.all_data()
    dens_ex = ad.quantities.extrema(('gas', 'cell_volume'))
    temp_ex = ad.quantities.extrema(('gas', 'x'))
    plot = ParticlePlot(ds,
                        ('gas', 'cell_volume'),
                        ('gas', 'x'),
                        ('gas', 'path_element_x'))
    plot.set_log(('gas', 'cell_volume'), True)
    plot.set_log(('gas', 'x'), True)
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

    plot.set_log(('gas', 'cell_volume'), False)
    plot.set_log(('gas', 'x'), False)
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

def test_set_units():
    ds = fake_random_ds(16, nprocs=4, particles=16 ** 3,
                        fields=("density", "temperature"),
                        units=('g/cm**3', 'K'))
    sp = ds.sphere("max", (1.0, "Mpc"))
    pp = ParticlePhasePlot(sp, ("gas", "density"), ("gas", "temperature"),
                           ("gas", "cell_mass"))
    # make sure we can set the units using the tuple without erroring out
    pp.set_unit(("gas", "cell_mass"), "Msun")

def test_switch_ds():
    """
    Tests the _switch_ds() method for ParticleProjectionPlots that as of
    25th October 2017 requires a specific hack in plot_container.py
    """
    ds = fake_particle_ds()
    ds2 = fake_particle_ds()

    plot = ParticlePlot(
        ds,
        ("gas", "path_element_x"),
        ("gas", "path_element_y"),
        ("gas", "cell_volume"),
    )

    plot._switch_ds(ds2)

    return

@attr(ANSWER_TEST_TAG)
def test_particle_plot():
    test_ds = fake_particle_ds()
    for dim in range(3):
        pplot = ParticleProjectionPlot(test_ds, dim, "particle_mass")
        test_prefix = "particle_mass_%s" % dim
        yield compare(test_ds, pplot, test_prefix=test_prefix,
                      test_name="particle_plot")

@attr(ANSWER_TEST_TAG)
def test_particle_plot_ds():
    test_ds = fake_particle_ds()
    ds_region = test_ds.region([0.5] * 3, [0.4] * 3, [0.6] * 3)
    for dim in range(3):
        pplot_ds = ParticleProjectionPlot(test_ds, dim, "particle_mass",
                                          data_source=ds_region)
        test_prefix = "particle_mass_%s" % dim
        yield compare(test_ds, pplot_ds, test_prefix=test_prefix,
                      test_name="particle_plot_ds")

@attr(ANSWER_TEST_TAG)
def test_particle_plot_c():
    test_ds = fake_particle_ds()
    for center in CENTER_SPECS:
        for dim in range(3):
            pplot_c = ParticleProjectionPlot(
                test_ds, dim, "particle_mass", center=center)
            test_prefix = "particle_mass_%s_%s" % (center, dim)
            yield compare(test_ds, pplot_c, test_prefix=test_prefix,
                          test_name="particle_plot_c")

@attr(ANSWER_TEST_TAG)
def test_particle_plot_wf():
    test_ds = fake_particle_ds()
    for dim in range(3):
        for weight_field in WEIGHT_FIELDS:
            pplot_wf = ParticleProjectionPlot(test_ds, dim, "particle_mass",
                                              weight_field=weight_field)
            test_prefix = "particle_mass_%s_%s" % (dim, weight_field)
            yield compare(test_ds, pplot_wf, test_prefix=test_prefix,
                          test_name="particle_plot_wf")

class TestParticleProjectionPlotSave(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.curdir = os.getcwd()
        os.chdir(self.tmpdir)

    def tearDown(self):
        os.chdir(self.curdir)
        shutil.rmtree(self.tmpdir)

    def test_creation_with_width(self):
        test_ds = fake_particle_ds()
        for width, (xlim, ylim, pwidth, aun) in WIDTH_SPECS.items():
            plot = ParticleProjectionPlot(
                test_ds, 0, 'particle_mass', width=width)

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
    x_field = ('all', 'particle_position_x')
    y_field = ('all', 'particle_position_y')
    z_field = ('all', 'particle_velocity_x')

    plot = ParticlePlot(ds, x_field, y_field)
    assert isinstance(plot, ParticleProjectionPlot)

    plot = ParticlePlot(ds, y_field, x_field)
    assert isinstance(plot, ParticleProjectionPlot)

    plot = ParticlePlot(ds, x_field, z_field)
    assert isinstance(plot, ParticlePhasePlot)
