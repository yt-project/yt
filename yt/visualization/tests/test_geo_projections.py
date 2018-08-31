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
from yt.testing import ANSWER_TEST_TAG, fake_amr_ds, requires_module
from yt.utilities.answer_testing.framework import GenericImageTest
from yt.visualization.geo_plot_utils import transform_list, get_mpl_transform

transform_list.remove('UTM')

def setup():
    """Test specific setup."""
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"

def compare(ds, field, idir, test_prefix, test_name, transform,
            decimals=12, annotate=False):
    def slice_image(filename_prefix):
        sl = yt.SlicePlot(ds, idir, field)
        sl.set_mpl_projection(transform)
        if annotate:
            sl._setup_plots()
            sl.annotate_mesh_lines()
        sl.set_log('all', False)
        image_file = sl.save(filename_prefix)
        return image_file

    slice_image.__name__ = "slice_{}".format(test_prefix)
    test = GenericImageTest(ds, slice_image, decimals)
    test.prefix = test_prefix
    test.answer_name = test_name
    return test

@requires_module("cartopy")
@attr(ANSWER_TEST_TAG)
def test_geo_slices_amr():
    ds = fake_amr_ds(geometry="geographic")
    for transform in transform_list:
        for field in ds.field_list:
            prefix = "%s_%s_%s" % (field[0], field[1], transform)
            yield compare(ds, field, 'altitude', test_prefix=prefix,
                          test_name="geo_slices_amr", transform=transform)

@requires_module("cartopy")
class TestGeoProjections(unittest.TestCase):

    def setUp(self):
        self.ds = fake_amr_ds(geometry="geographic")

    def tearDown(self):
        del self.ds
        del self.slc

    def test_projection_setup(self):

        from yt.utilities.on_demand_imports import _cartopy as cartopy
        self.slc = yt.SlicePlot(self.ds, "altitude", "Density")

        assert isinstance(self.slc._projection, cartopy.crs.PlateCarree)
        assert isinstance(self.slc._transform, cartopy.crs.PlateCarree)
        assert isinstance(self.slc._projection,
                          type(self.slc.plots['Density'].axes.projection))

    def test_projection_setup_modified(self):
        from yt.utilities.on_demand_imports import _cartopy as cartopy

        for transform in transform_list:
            self.slc = yt.SlicePlot(self.ds, "altitude", "Density",
                                    geo_projection=transform)
            proj_type = type(get_mpl_transform(transform))

            assert isinstance(self.slc._projection, proj_type)
            assert isinstance(self.slc._transform, cartopy.crs.PlateCarree)
            assert isinstance(self.slc.plots['Density'].axes.projection,
                              proj_type)

    def test_projection_transform(self):
        from yt.utilities.on_demand_imports import _cartopy as cartopy
        self.slc = yt.SlicePlot(self.ds, "altitude", "Density")

        for transform in transform_list:
            self.slc.set_mpl_projection(transform)
            proj_type = type(get_mpl_transform(transform))

            assert isinstance(self.slc._projection, proj_type)
            assert isinstance(self.slc._transform, cartopy.crs.PlateCarree)
            assert isinstance(self.slc.plots['Density'].axes.projection,
                              proj_type)
