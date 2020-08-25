import os
import tempfile
import unittest

import pytest

import yt
from yt.testing import fake_amr_ds, requires_module
from yt.utilities.answer_testing.answer_tests import generic_image
from yt.visualization.geo_plot_utils import get_mpl_transform, transform_list


def compare(ds, field, idir, projection, annotate=False):
    def slice_image():
        tmpfd, tmpfname = tempfile.mkstemp(suffix=".png")
        os.close(tmpfd)
        sl = yt.SlicePlot(ds, idir, field, origin="native")
        sl.set_mpl_projection(projection)
        if annotate:
            sl._setup_plots()
            sl.annotate_mesh_lines()
        sl.set_log("all", False)
        image_file = sl.save(tmpfname)
        return image_file

    gi = generic_image(slice_image)
    return gi


@pytest.mark.answer_test
@pytest.mark.usefixtures("temp_dir", "hashing")
class TestGeoSlicesAMR:
    answer_file = None
    saved_hashes = None

    @requires_module("cartopy")
    def test_geo_slices_amr(self, transform, field, ds):
        if transform not in ("UTM", "OSNI"):
            gi = compare(ds, field, "altitude")
            self.hashes.update({"generic_image": gi})


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
        assert isinstance(
            self.slc._projection, type(self.slc.plots["Density"].axes.projection)
        )

    def test_geo_projections(self):
        from yt.utilities.on_demand_imports import _cartopy as cartopy

        self.slc = yt.SlicePlot(self.ds, "altitude", "Density")

        for transform in transform_list:
            if transform == "UTM":
                # this requires special arguments so let's skip it
                continue
            if transform == "OSNI":
                # avoid crashes, see https://github.com/SciTools/cartopy/issues/1177
                continue
            self.slc.set_mpl_projection(transform)
            proj_type = type(get_mpl_transform(transform))

            assert isinstance(self.slc._projection, proj_type)
            assert isinstance(self.slc._transform, cartopy.crs.PlateCarree)
            assert isinstance(self.slc.plots["Density"].axes.projection, proj_type)

    def test_projection_object(self):
        from yt.utilities.on_demand_imports import _cartopy as cartopy

        shortlist = ["Orthographic", "PlateCarree", "Mollweide"]

        for transform in shortlist:
            projection = get_mpl_transform(transform)
            proj_type = type(projection)
            self.slc = yt.SlicePlot(self.ds, "altitude", "Density")
            self.slc.set_mpl_projection(projection)

            assert isinstance(self.slc._projection, proj_type)
            assert isinstance(self.slc._transform, cartopy.crs.PlateCarree)
            assert isinstance(self.slc.plots["Density"].axes.projection, proj_type)

    def test_nondefault_transform(self):
        from yt.utilities.on_demand_imports import _cartopy as cartopy

        axis = "altitude"
        self.ds.coordinates.data_transform[axis] = "Miller"
        self.slc = yt.SlicePlot(self.ds, axis, "Density")

        shortlist = ["Orthographic", "PlateCarree", "Mollweide"]

        for transform in shortlist:

            self.slc.set_mpl_projection(transform)
            proj_type = type(get_mpl_transform(transform))

            assert isinstance(self.slc._projection, proj_type)
            assert isinstance(self.slc._transform, cartopy.crs.Miller)
            assert self.ds.coordinates.data_projection[axis] == "Mollweide"
            assert self.ds.coordinates.data_transform[axis] == "Miller"
            assert isinstance(self.slc.plots["Density"].axes.projection, proj_type)


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
