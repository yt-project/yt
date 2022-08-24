import unittest

import numpy as np
from nose.plugins.attrib import attr

import yt
from yt.testing import ANSWER_TEST_TAG, fake_amr_ds, requires_module
from yt.utilities.answer_testing.framework import GenericImageTest
from yt.visualization.geo_plot_utils import get_mpl_transform, transform_list


def setup():
    """Test specific setup."""
    from yt.config import ytcfg

    ytcfg["yt", "internals", "within_testing"] = True


def compare(
    ds, field, idir, test_prefix, test_name, projection, decimals=12, annotate=False
):
    def slice_image(filename_prefix):
        sl = yt.SlicePlot(ds, idir, field, origin="native")
        sl.set_mpl_projection(projection)
        if annotate:
            sl._setup_plots()
            sl.annotate_mesh_lines()
        sl.set_log("all", False)
        image_file = sl.save(filename_prefix)
        return image_file

    slice_image.__name__ = f"slice_{test_prefix}"
    test = GenericImageTest(ds, slice_image, decimals)
    test.prefix = test_prefix
    test.answer_name = test_name
    return test


@requires_module("cartopy")
@attr(ANSWER_TEST_TAG)
def test_geo_slices_amr():
    ds = fake_amr_ds(geometry="geographic")
    for transform in transform_list:
        if transform == "UTM":
            # requires additional argument so we skip
            continue
        if transform == "OSNI":
            # avoid crashes, see https://github.com/SciTools/cartopy/issues/1177
            continue
        for field in ds.field_list:
            prefix = f"{field[0]}_{field[1]}_{transform}"
            yield compare(
                ds,
                field,
                "altitude",
                test_prefix=prefix,
                test_name="geo_slices_amr",
                projection=transform,
            )


class TestGeoProjections(unittest.TestCase):
    @requires_module("cartopy")
    def setUp(self):
        self.ds = fake_amr_ds(geometry="geographic")
        # switch off the log plot to avoid some unrelated matplotlib issues
        f = self.ds._get_field_info("stream", "Density")
        f.take_log = False

    @requires_module("cartopy")
    def tearDown(self):
        del self.ds

    @requires_module("cartopy")
    def test_geo_projection_setup(self):

        from yt.utilities.on_demand_imports import _cartopy as cartopy

        axis = "altitude"
        self.slc = yt.SlicePlot(self.ds, axis, ("stream", "Density"), origin="native")

        assert isinstance(self.slc._projection, cartopy.crs.Mollweide)
        assert isinstance(self.slc._transform, cartopy.crs.PlateCarree)
        assert self.ds.coordinates.data_projection[axis] == "Mollweide"
        assert self.ds.coordinates.data_transform[axis] == "PlateCarree"
        assert isinstance(
            self.slc._projection,
            type(self.slc.plots[("stream", "Density")].axes.projection),
        )

    @requires_module("cartopy")
    def test_geo_projections(self):
        from yt.utilities.on_demand_imports import _cartopy as cartopy

        self.slc = yt.SlicePlot(
            self.ds, "altitude", ("stream", "Density"), origin="native"
        )

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
            assert isinstance(
                self.slc.plots[("stream", "Density")].axes.projection, proj_type
            )

    @requires_module("cartopy")
    def test_projection_object(self):
        from yt.utilities.on_demand_imports import _cartopy as cartopy

        shortlist = ["Orthographic", "PlateCarree", "Mollweide"]

        for transform in shortlist:
            projection = get_mpl_transform(transform)
            proj_type = type(projection)
            self.slc = yt.SlicePlot(
                self.ds, "altitude", ("stream", "Density"), origin="native"
            )
            self.slc.set_mpl_projection(projection)

            assert isinstance(self.slc._projection, proj_type)
            assert isinstance(self.slc._transform, cartopy.crs.PlateCarree)
            assert isinstance(
                self.slc.plots[("stream", "Density")].axes.projection, proj_type
            )

    @requires_module("cartopy")
    def test_nondefault_transform(self):
        from yt.utilities.on_demand_imports import _cartopy as cartopy

        axis = "altitude"
        # Note: The Miller transform has an extent of approx. +/- 180 in x,
        # +/-132 in y (in Miller, x is longitude, y is a factor of latitude).
        # So by changing the projection in this way, the dataset goes from
        # covering the whole globe (+/- 90 latitude), to covering part of the
        # globe (+/-72 latitude). Totally fine for testing that the code runs,
        # but good to be aware of!
        self.ds.coordinates.data_transform[axis] = "Miller"
        self.slc = yt.SlicePlot(self.ds, axis, ("stream", "Density"), origin="native")

        shortlist = ["Orthographic", "PlateCarree", "Mollweide"]

        for transform in shortlist:

            self.slc.set_mpl_projection(transform)
            proj_type = type(get_mpl_transform(transform))

            assert isinstance(self.slc._projection, proj_type)
            assert isinstance(self.slc._transform, cartopy.crs.Miller)
            assert self.ds.coordinates.data_projection[axis] == "Mollweide"
            assert self.ds.coordinates.data_transform[axis] == "Miller"
            assert isinstance(
                self.slc.plots[("stream", "Density")].axes.projection, proj_type
            )

    @requires_module("cartopy")
    def test_extent(self):
        # checks that the axis extent is narrowed when doing a subselection
        axis = "altitude"
        slc = yt.SlicePlot(self.ds, axis, ("stream", "Density"), origin="native")
        ax = slc.plots[("stream", "Density")].axes
        full_extent = np.abs(ax.get_extent())

        slc = yt.SlicePlot(
            self.ds, axis, ("stream", "Density"), origin="native", width=(80.0, 50.0)
        )
        ax = slc.plots[("stream", "Density")].axes
        extent = np.abs(ax.get_extent())
        assert np.all(extent < full_extent)


class TestNonGeoProjections(unittest.TestCase):
    def setUp(self):
        self.ds = fake_amr_ds()

    def tearDown(self):
        del self.ds
        del self.slc

    def test_projection_setup(self):
        axis = "x"
        self.slc = yt.SlicePlot(self.ds, axis, ("stream", "Density"), origin="native")

        assert self.ds.coordinates.data_projection[axis] is None
        assert self.ds.coordinates.data_transform[axis] is None
        assert self.slc._projection is None
        assert self.slc._transform is None
