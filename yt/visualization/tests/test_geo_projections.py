"""
Tests for making unstructured mesh slices

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
from nose.plugins.attrib import attr

import yt
import unittest
import cartopy.crs
from yt.testing import ANSWER_TEST_TAG, fake_amr_ds
from yt.utilities.answer_testing.framework import GenericImageTest
from yt.visualization.geo_plot_utils import transform_list, get_mpl_transform

def setup():
    """Test specific setup."""
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"

def compare(ds, field, idir, test_prefix, test_name, transform,
            decimals=12, annotate=False):
    def slice_image(filename_prefix):
        sl = yt.SlicePlot(ds, idir, field)
        if annotate:
            sl.annotate_mesh_lines()
        sl.set_log('all', False)
        image_file = sl.save(filename_prefix)
        return image_file

    slice_image.__name__ = "slice_{}".format(test_prefix)
    test = GenericImageTest(ds, slice_image, decimals)
    test.prefix = test_prefix
    test.answer_name = test_name
    return test

@attr(ANSWER_TEST_TAG)
def test_geo_slices_amr():
    ds = fake_amr_ds(geometry="geographic")
    for transform in transform_list:
        for field in ds.field_list:
            prefix = "%s_%s_%s" % (field[0], field[1], transform)
            yield compare(ds, field, 'altitude', test_prefix=prefix,
                          test_name="geo_slices_amr", transform=transform)

class TestGeoProjections(unittest.TestCase):

    def setUp(self):
        self.ds = fake_amr_ds(geometry="geographic")
        self.slc = yt.SlicePlot(self.ds, "altitude", "Density")

    def tearDown(self):
        del self.ds
        del self.slc

    def test_projection_setup(self):

        assert isinstance(self.slc._projection, cartopy.crs.PlateCarree)
        assert isinstance(self.slc._transform, cartopy.crs.PlateCarree)
        assert isinstance(self.slc._projection,
                          type(self.slc.plots['Density'].axes.projection))

    def test_projection_transform(self):
        # remove UTM as a transform for testing
        # since it has a required arg
        transform_list.remove('UTM')
        for transform in transform_list:
            self.slc.set_mpl_projection(transform)
            proj_type = type(get_mpl_transform(transform))

            assert isinstance(self.slc._projection, proj_type)
            assert isinstance(self.slc._transform, cartopy.crs.PlateCarree)
            assert isinstance(self.slc.plots['Density'].axes.projection,
                              proj_type)
