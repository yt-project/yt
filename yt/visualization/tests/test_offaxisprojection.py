"""
Test for off_axis_projection and write_projection



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
import os
import os.path
import shutil
import tempfile
import unittest

from yt.mods import write_projection
from yt.testing import \
    fake_random_ds, assert_equal, expand_keywords
from yt.visualization.volume_rendering.api import off_axis_projection


class TestOffAxisProjection(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp()
        cls.curdir = os.getcwd()
        os.chdir(cls.tmpdir)

    @classmethod
    def tearDownClass(cls):
        os.chdir(cls.curdir)
        shutil.rmtree(cls.tmpdir)

    def test_oap(self):
        """Tests functionality of off_axis_projection and write_projection."""

        # args for off_axis_projection
        test_ds = fake_random_ds(64)
        c = test_ds.domain_center
        norm = [0.5, 0.5, 0.5]
        W = test_ds.arr([0.5,0.5,1.0], 'unitary')
        N = 256
        field = ("gas","density")
        oap_args = [test_ds, c, norm, W, N, field]

        # kwargs for off_axis_projection
        oap_kwargs = {}
        oap_kwargs['weight'] = (None, 'cell_mass')
        oap_kwargs['no_ghost'] = (True, False)
        oap_kwargs['interpolated'] = (False,)
        oap_kwargs['north_vector'] = ((1, 0, 0), (0, 0.5, 1.0))
        oap_kwargs_list = expand_keywords(oap_kwargs)

        # args or write_projection
        fn = "test_%d.png"

        # kwargs for write_projection
        wp_kwargs = {}
        wp_kwargs['colorbar'] = (True, False)
        wp_kwargs['colorbar_label'] = ('test')
        wp_kwargs['title'] = ('test')
        wp_kwargs['limits'] = (None, (1e3, 1e5))
        wp_kwargs['take_log'] = (True, False)
        wp_kwargs['figsize'] = ((8,6), [1,1])
        wp_kwargs['dpi'] = (100, 50)
        wp_kwargs['cmap_name'] = ('arbre', 'kelp')
        wp_kwargs_list = expand_keywords(wp_kwargs)

        # test all off_axis_projection kwargs and write_projection kwargs
        # make sure they are able to be projected, then remove and try next
        # iteration
        for i, oap_kwargs in enumerate(oap_kwargs_list):
            image = off_axis_projection(*oap_args, **oap_kwargs)
            for wp_kwargs in wp_kwargs_list:
                write_projection(image, fn % i, **wp_kwargs)
                assert_equal(os.path.exists(fn % i), True)

        # Test remaining parameters of write_projection
        write_projection(image, "test_2", xlabel="x-axis", ylabel="y-axis")
        write_projection(image, "test_3.pdf", xlabel="x-axis", ylabel="y-axis")
        write_projection(image, "test_4.eps", xlabel="x-axis", ylabel="y-axis")
