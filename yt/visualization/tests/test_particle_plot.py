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
import itertools
import os
import tempfile
import shutil
import unittest
from yt.data_objects.profiles import create_profile
from yt.extern.parameterized import parameterized, param
from yt.visualization.tests.test_plotwindow import \
    assert_fname, WIDTH_SPECS
from yt.testing import \
    fake_particle_ds, assert_array_almost_equal
from yt.visualization.api import \
    ParticleProjectionPlot, ParticlePhasePlot
from yt.units.yt_array import YTArray


def setup():
    """Test specific setup."""
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"


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


class TestParticlePhasePlotSave(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
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

        cls.particle_phases = particle_phases
        cls.ds = test_ds

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.curdir = os.getcwd()
        os.chdir(self.tmpdir)

    def tearDown(self):
        os.chdir(self.curdir)
        shutil.rmtree(self.tmpdir)

    @parameterized.expand(param.explicit((fname, )) for fname in TEST_FLNMS)
    def test_particle_phase_plot(self, fname):
        for p in self.particle_phases:
            assert assert_fname(p.save(fname)[0])

    def test_ipython_repr(self):
        self.particle_phases[0]._repr_html_()


class TestParticleProjectionPlotSave(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        test_ds = fake_particle_ds()
        ds_region = test_ds.region([0.5] * 3, [0.4] * 3, [0.6] * 3)
        pplots = []
        pplots_ds = []
        pplots_c = []
        pplots_wf = []
        pplots_w = {}
        pplots_m = []
        for dim in range(3):
            pplots.append(ParticleProjectionPlot(test_ds, dim,
                                                 "particle_mass"))

            pplots_ds.append(ParticleProjectionPlot(test_ds, dim,
                                                    "particle_mass",
                                                    data_source=ds_region))

        for center in CENTER_SPECS:
            pplots_c.append(ParticleProjectionPlot(test_ds, dim,
                                                   "particle_mass",
                                                    center=center))

        for width in WIDTH_SPECS:
            pplots_w[width] = ParticleProjectionPlot(test_ds, dim,
                                                     'particle_mass',
                                                     width=width)
        for wf in WEIGHT_FIELDS:
            pplots_wf.append(ParticleProjectionPlot(test_ds, dim,
                                                    "particle_mass",
                                                    weight_field=wf))
        cls.pplots = pplots
        cls.pplots_ds = pplots_ds
        cls.pplots_c = pplots_c
        cls.pplots_wf = pplots_wf
        cls.pplots_w = pplots_w
        cls.pplots_m = pplots_m

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
    def test_particle_plot(self, dim, fname):
        assert assert_fname(self.pplots[dim].save(fname)[0])

    @parameterized.expand([(0, ), (1, ), (2, )])
    def test_particle_plot_ds(self, dim):
        self.pplots_ds[dim].save()

    @parameterized.expand([(i, ) for i in range(len(CENTER_SPECS))])
    def test_particle_plot_c(self, dim):
        self.pplots_c[dim].save()

    @parameterized.expand([(i, ) for i in range(len(WEIGHT_FIELDS))])
    def test_particle_plot_wf(self, dim):
        self.pplots_wf[dim].save()

    @parameterized.expand(
        param.explicit((width, ))
        for width in WIDTH_SPECS)
    def test_creation_with_width(self, width):
        xlim, ylim, pwidth, aun = WIDTH_SPECS[width]
        plot = self.pplots_w[width]

        xlim = [plot.ds.quan(el[0], el[1]) for el in xlim]
        ylim = [plot.ds.quan(el[0], el[1]) for el in ylim]
        pwidth = [plot.ds.quan(el[0], el[1]) for el in pwidth]

        [assert_array_almost_equal(px, x, 14) for px, x in zip(plot.xlim, xlim)]
        [assert_array_almost_equal(py, y, 14) for py, y in zip(plot.ylim, ylim)]
        [assert_array_almost_equal(pw, w, 14) for pw, w in zip(plot.width, pwidth)]
