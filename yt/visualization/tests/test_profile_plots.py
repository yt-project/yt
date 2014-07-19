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
import itertools
import os
import tempfile
import shutil
import unittest
from yt.data_objects.profiles import create_profile
from yt.extern.parameterized import\
    parameterized, param
from yt.testing import fake_random_ds
from yt.visualization.profile_plotter import \
    ProfilePlot, PhasePlot
from yt.visualization.tests.test_plotwindow import \
    assert_fname, TEST_FLNMS

class TestProfilePlotSave(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        fields = ('density', 'temperature', 'velocity_x', 'velocity_y',
                  'velocity_z')
        units = ('g/cm**3', 'K', 'cm/s', 'cm/s', 'cm/s')
        test_ds = fake_random_ds(64, fields=fields, units=units)
        regions = [test_ds.region([0.5]*3, [0.4]*3, [0.6]*3), test_ds.all_data()]
        profiles = []
        phases = []
        pr_fields = [('density', 'temperature'), ('density', 'velocity_x'),
                     ('temperature', 'cell_mass'), ('density', 'radius'),
                     ('velocity_magnitude', 'cell_mass')]
        ph_fields = [('density', 'temperature', 'cell_mass'),
                     ('density', 'velocity_x', 'cell_mass'),
                     ('radius', 'temperature', 'velocity_magnitude')]
        for reg in regions:
            for x_field, y_field in pr_fields:
                profiles.append(ProfilePlot(reg, x_field, y_field))
                profiles.append(ProfilePlot(reg, x_field, y_field,
                                            fractional=True, accumulation=True))
                p1d = create_profile(reg, x_field, y_field)
                profiles.append(ProfilePlot.from_profiles(p1d))
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
        cls.profiles = profiles
        cls.phases = phases
        cls.ds = test_ds

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.curdir = os.getcwd()
        os.chdir(self.tmpdir)

    def tearDown(self):
        os.chdir(self.curdir)
        shutil.rmtree(self.tmpdir)

    @parameterized.expand(param.explicit((fname, )) for fname in TEST_FLNMS)
    def test_profile_plot(self, fname):
        for p in self.profiles:
            yield assert_fname(p.save(fname)[0])

    @parameterized.expand(param.explicit((fname, )) for fname in TEST_FLNMS)
    def test_phase_plot(self, fname):
        for p in self.phases:
            assert assert_fname(p.save(fname)[0])
