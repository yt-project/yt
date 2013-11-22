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
import itertools
import os
import tempfile
import shutil
import unittest
from yt.extern.parameterized import parameterized, param
from yt.testing import \
    fake_random_pf, assert_equal, assert_rel_equal
from yt.utilities.answer_testing.framework import \
    requires_pf, data_dir_load, PlotWindowAttributeTest
from yt.visualization.api import \
    SlicePlot, ProjectionPlot, OffAxisSlicePlot, OffAxisProjectionPlot


def setup():
    """Test specific setup."""
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"


def assert_fname(fname):
    """Function that checks file type using libmagic"""
    if fname is None:
        return

    with open(fname, 'rb') as fimg:
        data = fimg.read()
    data = str(data)
    image_type = ''

    # see http://www.w3.org/TR/PNG/#5PNG-file-signature
    if data.startswith('\211PNG\r\n\032\n'):
        image_type = '.png'
    # see http://www.mathguide.de/info/tools/media-types/image/jpeg
    elif data.startswith('\377\330'):
        image_type = '.jpeg'
    elif data.startswith('%!PS-Adobe'):
        if 'EPSF' in data[:data.index('\n')]:
            image_type = '.eps'
        else:
            image_type = '.ps'
    elif data.startswith('%PDF'):
        image_type = '.pdf'

    return image_type == os.path.splitext(fname)[1]


TEST_FLNMS = [None, 'test.png', 'test.eps',
              'test.ps', 'test.pdf']
M7 = "DD0010/moving7_0010"
WT = "WindTunnel/windtunnel_4lev_hdf5_plt_cnt_0030"

ATTR_ARGS = {"pan": [(((0.1, 0.1), ), {})],
             "pan_rel": [(((0.1, 0.1), ), {})],
             "set_axes_unit": [(("kpc", ), {}),
                               (("Mpc", ), {}),
                               ((("kpc", "kpc"),), {}),
                               ((("kpc", "Mpc"),), {})],
             "set_buff_size": [((1600, ), {}),
                               (((600, 800), ), {})],
             "set_center": [(((0.4, 0.3), ), {})],
             "set_cmap": [(('Density', 'RdBu'), {}),
                          (('Density', 'kamae'), {})],
             "set_font": [(({'family': 'sans-serif', 'style': 'italic',
                             'weight': 'bold', 'size': 24}, ), {})],
             "set_log": [(('Density', False), {})],
             "set_window_size": [((7.0, ), {})],
             "set_zlim": [(('Density', 1e-25, 1e-23), {}),
                          (('Density', 1e-25, None), {'dynamic_range': 4})],
             "zoom": [((10, ), {})]}


@requires_pf(M7)
def test_attributes():
    """Test plot member functions that aren't callbacks"""
    plot_field = 'Density'
    decimals = 3

    pf = data_dir_load(M7)
    for ax in 'xyz':
        for attr_name in ATTR_ARGS.keys():
            for args in ATTR_ARGS[attr_name]:
                test = PlotWindowAttributeTest(pf, plot_field, ax, attr_name,
                                               args, decimals)
                test_attributes.__name__ = test.description
                yield test


@requires_pf(WT)
def test_attributes_wt():
    plot_field = 'Density'
    decimals = 3

    pf = data_dir_load(WT)
    ax = 'z'
    for attr_name in ATTR_ARGS.keys():
        for args in ATTR_ARGS[attr_name]:
            yield PlotWindowAttributeTest(pf, plot_field, ax, attr_name,
                                          args, decimals)


class TestSetWidth(unittest.TestCase):

    pf = None

    def setUp(self):
        if self.pf is None:
            self.pf = fake_random_pf(64)
            self.slc = SlicePlot(self.pf, 0, 'Density')

    def _assert_15kpc(self):
        assert_rel_equal([self.slc.xlim, self.slc.ylim, self.slc.width],
                         [(-7.5 / self.pf['kpc'], 7.5 / self.pf['kpc']),
                          (-7.5 / self.pf['kpc'], 7.5 / self.pf['kpc']),
                          (15.0 / self.pf['kpc'], 15. / self.pf['kpc'])], 15)

    def _assert_15_10kpc(self):
        assert_rel_equal([self.slc.xlim, self.slc.ylim, self.slc.width],
                         [(-7.5 / self.pf['kpc'], 7.5 / self.pf['kpc']),
                          (-5.0 / self.pf['kpc'], 5.0 / self.pf['kpc']),
                          (15.0 / self.pf['kpc'], 10. / self.pf['kpc'])], 15)

    def test_set_width_one(self):
        assert_equal([self.slc.xlim, self.slc.ylim, self.slc.width],
                     [(0.0, 1.0), (0.0, 1.0), (1.0, 1.0)])

    def test_set_width_nonequal(self):
        self.slc.set_width((0.5, 0.8))
        assert_rel_equal([self.slc.xlim, self.slc.ylim, self.slc.width],
                         [(0.25, 0.75), (0.1, 0.9), (0.5, 0.8)], 15)

    def test_twoargs_eq(self):
        self.slc.set_width(15, 'kpc')
        self._assert_15kpc()

    def test_tuple_eq(self):
        self.slc.set_width((15, 'kpc'))
        self._assert_15kpc()

    def test_tuple_of_tuples_neq(self):
        self.slc.set_width(((15, 'kpc'), (10, 'kpc')))
        self._assert_15_10kpc()

    def test_tuple_of_tuples_neq2(self):
        self.slc.set_width(((15, 'kpc'), (10000, 'pc')))
        self._assert_15_10kpc()

    def test_pair_of_tuples_neq(self):
        self.slc.set_width((15, 'kpc'), (10000, 'pc'))
        self._assert_15_10kpc()


class TestPlotWindowSave(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        test_pf = fake_random_pf(64)
        normal = [1, 1, 1]
        ds_region = test_pf.h.region([0.5] * 3, [0.4] * 3, [0.6] * 3)
        projections = []
        projections_ds = []
        for dim in range(3):
            projections.append(ProjectionPlot(test_pf, dim, 'Density'))
            projections_ds.append(ProjectionPlot(test_pf, dim, 'Density',
                                                 data_source=ds_region))

        cls.slices = [SlicePlot(test_pf, dim, 'Density') for dim in range(3)]
        cls.projections = projections
        cls.projections_ds = projections_ds
        cls.offaxis_slice = OffAxisSlicePlot(test_pf, normal, 'Density')
        cls.offaxis_proj = OffAxisProjectionPlot(test_pf, normal, 'Density')

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
    def test_slice_plot(self, dim, fname):
        assert assert_fname(self.slices[dim].save(fname)[0])

    @parameterized.expand(
        param.explicit(item)
        for item in itertools.product(range(3), TEST_FLNMS))
    def test_projection_plot(self, dim, fname):
        assert assert_fname(self.projections[dim].save(fname)[0])

    @parameterized.expand([(0, ), (1, ), (2, )])
    def test_projection_plot_ds(self, dim):
        self.projections_ds[dim].save()

    @parameterized.expand(
        param.explicit((fname, ))
        for fname in TEST_FLNMS)
    def test_offaxis_slice_plot(self, fname):
        assert assert_fname(self.offaxis_slice.save(fname)[0])

    @parameterized.expand(
        param.explicit((fname, ))
        for fname in TEST_FLNMS)
    def test_offaxis_projection_plot(self, fname):
        assert assert_fname(self.offaxis_proj.save(fname)[0])
