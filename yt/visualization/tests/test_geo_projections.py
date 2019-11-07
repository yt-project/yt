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
import os
import tempfile

import pytest

import yt
import unittest
from yt.testing import \
    fake_amr_ds, \
    requires_module
import yt.utilities.answer_testing.framework as fw
from yt.visualization.geo_plot_utils import transform_list, get_mpl_transform


@pytest.mark.answer_test
@pytest.mark.usefixtures('temp_dir', 'answer_file', 'hashing')
class TestGeoSlicesAMR(fw.AnswerTest):
    @requires_module("cartopy")
    def test_geo_slices_amr(self, transform, field, ds):
        if transform not in ('UTM', 'OSNI'):
            tmpfd, tmpfname = tempfile.mkstemp(suffix='.png')
            os.close(tmpfd)
            sl = yt.SlicePlot(ds, 'altitude', field)
            sl.set_mpl_projection(transform)
            sl.set_log('all', False)
            sl.save(tmpfname)
            gi_hd = self.generic_image_test(tmpfname)
            self.hashes.update({'generic_image' : gi_hd})

@requires_module("cartopy")
class TestGeoProjections(unittest.TestCase):

    def setUp(self):
        self.ds = fake_amr_ds(geometry="geographic")

    def tearDown(self):
        del self.ds
        del self.slc

    def test_geo_projection_setup(self):

        from yt.utilities.on_demand_imports import _cartopy as cartopy
        axis = "altitude"
        self.slc = yt.SlicePlot(self.ds, axis, "Density")

        assert isinstance(self.slc._projection, cartopy.crs.Mollweide)
        assert isinstance(self.slc._transform, cartopy.crs.PlateCarree)
        assert self.ds.coordinates.data_projection[axis] == "Mollweide"
        assert self.ds.coordinates.data_transform[axis] == "PlateCarree"
        assert isinstance(self.slc._projection,
                          type(self.slc.plots['Density'].axes.projection))

    def test_geo_projections(self):
        from yt.utilities.on_demand_imports import _cartopy as cartopy
        self.slc = yt.SlicePlot(self.ds, "altitude", "Density")

        for transform in transform_list:
            if transform == 'UTM':
                # this requires special arguments so let's skip it
                continue
            if transform == 'OSNI':
                # avoid crashes, see https://github.com/SciTools/cartopy/issues/1177
                continue
            self.slc.set_mpl_projection(transform)
            proj_type = type(get_mpl_transform(transform))

            assert isinstance(self.slc._projection, proj_type)
            assert isinstance(self.slc._transform, cartopy.crs.PlateCarree)
            assert isinstance(self.slc.plots['Density'].axes.projection,
                              proj_type)

    def test_projection_object(self):
        from yt.utilities.on_demand_imports import _cartopy as cartopy
        shortlist = ['Orthographic', 'PlateCarree', 'Mollweide']

        for transform in shortlist:
            projection = get_mpl_transform(transform)
            proj_type = type(projection)
            self.slc = yt.SlicePlot(self.ds, "altitude", "Density")
            self.slc.set_mpl_projection(projection)

            assert isinstance(self.slc._projection, proj_type)
            assert isinstance(self.slc._transform, cartopy.crs.PlateCarree)
            assert isinstance(self.slc.plots['Density'].axes.projection,
                              proj_type)

    def test_nondefault_transform(self):
        from yt.utilities.on_demand_imports import _cartopy as cartopy
        axis = "altitude"
        self.ds.coordinates.data_transform[axis] = "Miller"
        self.slc = yt.SlicePlot(self.ds, axis, "Density")

        shortlist = ['Orthographic', 'PlateCarree', 'Mollweide']

        for transform in shortlist:

            self.slc.set_mpl_projection(transform)
            proj_type = type(get_mpl_transform(transform))

            assert isinstance(self.slc._projection, proj_type)
            assert isinstance(self.slc._transform, cartopy.crs.Miller)
            assert self.ds.coordinates.data_projection[axis] == "Mollweide"
            assert self.ds.coordinates.data_transform[axis] == "Miller"
            assert isinstance(self.slc.plots['Density'].axes.projection,
                              proj_type)

class TestNonGeoProjections(unittest.TestCase):

    def setUp(self):
        self.ds = fake_amr_ds()

    def tearDown(self):
        del self.ds
        del self.slc

    def test_projection_setup(self):
        axis = "x"
        self.slc = yt.SlicePlot(self.ds, axis, "Density")

        assert self.ds.coordinates.data_projection[axis] is None
        assert self.ds.coordinates.data_transform[axis] is None
        assert self.slc._projection is None
        assert self.slc._transform is None

