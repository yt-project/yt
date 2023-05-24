import pytest

import yt
from yt.testing import fake_amr_ds, requires_module_pytest as requires_module


class TestGeoTransform:
    # Cartopy require pykdtree *or* scipy to enable the projections
    # we test here. We require scipy for simplicity because it offers
    # better support for various platforms via PyPI at the time of writing.abs

    # the following projections are skipped (reason)
    # - UTM (requires additional arguments)
    # - OSNI (avoid crashes, see https://github.com/SciTools/cartopy/issues/1177)

    @classmethod
    def setup_class(cls):
        cls.ds = fake_amr_ds(geometry="geographic")

    @requires_module("cartopy", "scipy")
    @pytest.mark.parametrize(
        "transform",
        [
            "PlateCarree",
            "LambertConformal",
            "LambertCylindrical",
            "Mercator",
            "Miller",
            "Mollweide",
            "Orthographic",
            "Robinson",
            "Stereographic",
            "TransverseMercator",
            "InterruptedGoodeHomolosine",
            "RotatedPole",
            "OSGB",
            "EuroPP",
            "Geostationary",
            "Gnomonic",
            "NorthPolarStereo",
            "SouthPolarStereo",
            "AlbersEqualArea",
            "AzimuthalEquidistant",
            "Sinusoidal",
            "NearsidePerspective",
            "LambertAzimuthalEqualArea",
        ],
    )
    @pytest.mark.mpl_image_compare
    def test_geo_tranform(self, transform):
        field = ("stream", "Density")
        sl = yt.SlicePlot(self.ds, "altitude", field, origin="native")
        sl.set_mpl_projection(transform)
        sl.set_log(field, False)
        sl.render()
        return sl.plots[field].figure
