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

from yt.data_objects.profiles import create_profile
from yt.visualization.tests.test_plotwindow import \
    WIDTH_SPECS, ATTR_ARGS
from yt.convenience import load
from yt.data_objects.particle_filters import add_particle_filter
from yt.frontends.stream.data_structures import load_particles
from yt.testing import \
    fake_particle_ds, \
    assert_array_almost_equal, \
    requires_file, \
    assert_allclose, \
    assert_fname, \
    assert_array_equal
from yt.utilities.answer_testing.framework import \
    requires_ds, \
    data_dir_load, \
    PlotWindowAttributeTest, \
    PhasePlotAttributeTest
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

PHASE_ATTR_ARGS = {"annotate_text": [(((5e-29, 5e7), "Hello YT"), {}), 
                               (((5e-29, 5e7), "Hello YT"), {'color':'b'})],
                   "set_title": [(('particle_mass', 'A phase plot.'), {})],
                   "set_log": [(('particle_mass', False), {})],
                   "set_unit": [(('particle_mass', 'Msun'), {})],
                   "set_xlim": [((-4e7, 4e7), {})],
                   "set_ylim": [((-4e7, 4e7), {})]}

TEST_FLNMS = [None, 'test', 'test.png', 'test.eps',
              'test.ps', 'test.pdf']

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


g30 = "IsolatedGalaxy/galaxy0030/galaxy0030"

@requires_ds(g30, big_data=True)
def test_particle_projection_answers():
    '''

    This iterates over the all the plot modification functions in 
    PROJ_ATTR_ARGS. Each time, it compares the images produced by 
    ParticleProjectionPlot to the gold standard.
    

    '''

    plot_field = 'particle_mass'
    decimals = 12
    ds = data_dir_load(g30)
    for ax in 'xyz':
        for attr_name in PROJ_ATTR_ARGS.keys():
            for args in PROJ_ATTR_ARGS[attr_name]:
                test = PlotWindowAttributeTest(ds, plot_field, ax, 
                                               attr_name,
                                               args, decimals, 
                                               'ParticleProjectionPlot')
                test_particle_projection_answers.__name__ = test.description
                yield test


@requires_ds(g30, big_data=True)
def test_particle_projection_filter():
    '''

    This tests particle projection plots for filter fields.
    

    '''

    def formed_star(pfilter, data):
        filter = data["all", "creation_time"] > 0
        return filter

    add_particle_filter("formed_star", function=formed_star, filtered_type='all',
                        requires=["creation_time"])

    plot_field = ('formed_star', 'particle_mass')

    decimals = 12
    ds = data_dir_load(g30)
    ds.add_particle_filter('formed_star')
    for ax in 'xyz':
        attr_name = "set_log"
        for args in PROJ_ATTR_ARGS[attr_name]:
            test = PlotWindowAttributeTest(ds, plot_field, ax,
                                           attr_name,
                                           args, decimals,
                                           'ParticleProjectionPlot')
            test_particle_projection_filter.__name__ = test.description
            yield test


@requires_ds(g30, big_data=True)
def test_particle_phase_answers():
    '''

    This iterates over the all the plot modification functions in 
    PHASE_ATTR_ARGS. Each time, it compares the images produced by 
    ParticlePhasePlot to the gold standard.

    '''

    decimals = 12
    ds = data_dir_load(g30)

    x_field = 'particle_velocity_x'
    y_field = 'particle_velocity_y'
    z_field = 'particle_mass'
    for attr_name in PHASE_ATTR_ARGS.keys():
        for args in PHASE_ATTR_ARGS[attr_name]:
            test = PhasePlotAttributeTest(ds, x_field, y_field, z_field,
                                          attr_name, args, decimals,
                                          'ParticlePhasePlot')
                
            test_particle_phase_answers.__name__ = test.description
            yield test

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
        for p in particle_phases:
            for fname in TEST_FLNMS:
                assert assert_fname(p.save(fname)[0])

tgal = 'TipsyGalaxy/galaxy.00300'
@requires_file(tgal)
def test_particle_phase_plot_semantics():
    ds = load(tgal)
    ad = ds.all_data()
    dens_ex = ad.quantities.extrema(('Gas', 'density'))
    temp_ex = ad.quantities.extrema(('Gas', 'temperature'))
    plot = ParticlePlot(ds,
                        ('Gas', 'density'),
                        ('Gas', 'temperature'),
                        ('Gas', 'particle_mass'))
    plot.set_log('density', True)
    plot.set_log('temperature', True)
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

    plot.set_log('density', False)
    plot.set_log('temperature', False)
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
    pp = ParticlePhasePlot(sp, ("Gas", "density"), ("Gas", "temperature"), ("Gas", "particle_mass"))
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
                assert assert_fname(pplot.save(fname)[0])

    def test_particle_plot_ds(self):
        test_ds = fake_particle_ds()
        ds_region = test_ds.region([0.5] * 3, [0.4] * 3, [0.6] * 3)
        for dim in range(3):
            pplot_ds = ParticleProjectionPlot(
                test_ds, dim, "particle_mass", data_source=ds_region)
            pplot_ds.save()

    def test_particle_plot_c(self):
        test_ds = fake_particle_ds()
        for center in CENTER_SPECS:
            for dim in range(3):
                pplot_c = ParticleProjectionPlot(
                    test_ds, dim, "particle_mass", center=center)
                pplot_c.save()

    def test_particle_plot_wf(self):
        test_ds = fake_particle_ds()
        for dim in range(3):
            for weight_field in WEIGHT_FIELDS:
                pplot_wf = ParticleProjectionPlot(
                    test_ds, dim, "particle_mass", weight_field=weight_field)
                pplot_wf.save()

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

def test_projection_orientation():

    npart = 7

    # one particle at the origin, one particle along x-axis, two along y,
    # three along z
    data = {
        'particle_position_x': (
            np.array([0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]), 'cm'),
        'particle_position_y': (
            np.array([0.0, 0.0, 1.0, 2.0, 0.0, 0.0, 0.0]), 'cm'),
        'particle_position_z': (
            np.array([0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0]), 'cm'),
        'particle_mass': (np.ones(npart), 'g'),
        'particle_velocity_x': (np.zeros(npart), 'cm/s'),
        'particle_velocity_y': (np.zeros(npart), 'cm/s'),
        'particle_velocity_z': (np.zeros(npart), 'cm/s'),
    }

    bbox = np.array([[-4, 4], [-4, 4], [-4, 4]])

    ds = load_particles(data=data, length_unit=1.0, bbox=bbox)

    ANSWERS = {
        # 1 at origin, 3 particles along z and 2 particles along y
        ('x', (100, 100)): ([50, 50, 50, 62, 75, 87], [50, 62, 75, 50, 50, 50]),
        ('x', (100, 50 )): ([25, 25, 25, 31, 37, 43], [50, 62, 75, 50, 50, 50]),
        ('x', (50 , 100)): ([50, 50, 50, 62, 75, 87], [25, 31, 37, 25, 25, 25]),

        # 1 at origin, 3 particles long z and 1 particle along x
        ('y', (100, 100)): ([50, 50, 50, 50, 62], [50, 62, 75, 87, 50]),
        ('y', (100, 50 )): ([25, 25, 25, 25, 31], [50, 62, 75, 87, 50]),
        ('y', (50 , 100)): ([50, 50, 50, 50, 62], [25, 31, 37, 43, 25]),

        # 1 at origin, 2 particles along y and 1 particle along x
        ('z', (100, 100)): ([50, 50, 62, 75], [50, 62, 50, 50]),
        ('z', (100, 50 )): ([25, 25, 31, 37], [50, 62, 50, 50]),
        ('z', (50 , 100)): ([50, 50, 62, 75], [25, 31, 25, 25]),
    }

    for shape in [(100, 100), (100, 50), (50, 100)]:
        for ax in 'xyz':
            p = ParticleProjectionPlot(ds, ax, 'particle_mass', origin='native')
            p.set_buff_size(shape)
            p.set_log('particle_mass', False)
            p.set_zlim('particle_mass', 0, 4)
            image = p.frb['io', 'particle_mass']
            assert_array_equal(ANSWERS[ax, shape], np.where(np.isfinite(image)))
