import pytest

import yt
from yt.testing import fake_amr_ds, requires_module


class TestGeoTransform:
    # the following projections are skipped (reason)
    # - UTM (requires additional arguments)
    # - OSNI (avoid crashes, see https://github.com/SciTools/cartopy/issues/1177)

    @classmethod
    def setup_class(cls):
        cls.ds = fake_amr_ds(geometry="geographic")

    @requires_module("cartopy")
    @pytest.mark.mpl_image_compare(filename="geo_tranform_PlateCarree.png")
    def test_geo_tranform_PlateCarree(self):
        field = ("stream", "Density")
        sl = yt.SlicePlot(self.ds, "altitude", field, origin="native")
        sl.set_mpl_projection("PlateCarree")
        sl.set_log(field, False)
        return sl.plots[field].figure

    @requires_module("cartopy")
    @pytest.mark.mpl_image_compare(filename="geo_tranform_LambertConformal.png")
    def test_geo_tranform_LambertConformal(self):
        field = ("stream", "Density")
        sl = yt.SlicePlot(self.ds, "altitude", field, origin="native")
        sl.set_mpl_projection("LambertConformal")
        sl.set_log(field, False)
        return sl.plots[field].figure

    @requires_module("cartopy")
    @pytest.mark.mpl_image_compare(filename="geo_tranform_LambertCylindrical.png")
    def test_geo_tranform_LambertCylindrical(self):
        field = ("stream", "Density")
        sl = yt.SlicePlot(self.ds, "altitude", field, origin="native")
        sl.set_mpl_projection("LambertCylindrical")
        sl.set_log(field, False)
        return sl.plots[field].figure

    @requires_module("cartopy")
    @pytest.mark.mpl_image_compare(filename="geo_tranform_Mercator.png")
    def test_geo_tranform_Mercator(self):
        field = ("stream", "Density")
        sl = yt.SlicePlot(self.ds, "altitude", field, origin="native")
        sl.set_mpl_projection("Mercator")
        sl.set_log(field, False)
        return sl.plots[field].figure

    @requires_module("cartopy")
    @pytest.mark.mpl_image_compare(filename="geo_tranform_Miller.png")
    def test_geo_tranform_Miller(self):
        field = ("stream", "Density")
        sl = yt.SlicePlot(self.ds, "altitude", field, origin="native")
        sl.set_mpl_projection("Miller")
        sl.set_log(field, False)
        return sl.plots[field].figure

    @requires_module("cartopy")
    @pytest.mark.mpl_image_compare(filename="geo_tranform_Mollweide.png")
    def test_geo_tranform_Mollweide(self):
        field = ("stream", "Density")
        sl = yt.SlicePlot(self.ds, "altitude", field, origin="native")
        sl.set_mpl_projection("Mollweide")
        sl.set_log(field, False)
        return sl.plots[field].figure

    @requires_module("cartopy")
    @pytest.mark.mpl_image_compare(filename="geo_tranform_Orthographic.png")
    def test_geo_tranform_Orthographic(self):
        field = ("stream", "Density")
        sl = yt.SlicePlot(self.ds, "altitude", field, origin="native")
        sl.set_mpl_projection("Orthographic")
        sl.set_log(field, False)
        return sl.plots[field].figure

    @requires_module("cartopy")
    @pytest.mark.mpl_image_compare(filename="geo_tranform_Robinson.png")
    def test_geo_tranform_Robinson(self):
        field = ("stream", "Density")
        sl = yt.SlicePlot(self.ds, "altitude", field, origin="native")
        sl.set_mpl_projection("Robinson")
        sl.set_log(field, False)
        return sl.plots[field].figure

    @requires_module("cartopy")
    @pytest.mark.mpl_image_compare(filename="geo_tranform_Stereographic.png")
    def test_geo_tranform_Stereographic(self):
        field = ("stream", "Density")
        sl = yt.SlicePlot(self.ds, "altitude", field, origin="native")
        sl.set_mpl_projection("Stereographic")
        sl.set_log(field, False)
        return sl.plots[field].figure

    @requires_module("cartopy")
    @pytest.mark.mpl_image_compare(filename="geo_tranform_TransverseMercator.png")
    def test_geo_tranform_TransverseMercator(self):
        field = ("stream", "Density")
        sl = yt.SlicePlot(self.ds, "altitude", field, origin="native")
        sl.set_mpl_projection("TransverseMercator")
        sl.set_log(field, False)
        return sl.plots[field].figure

    @requires_module("cartopy")
    @pytest.mark.mpl_image_compare(
        filename="geo_tranform_InterruptedGoodeHomolosine.png"
    )
    def test_geo_tranform_InterruptedGoodeHomolosine(self):
        field = ("stream", "Density")
        sl = yt.SlicePlot(self.ds, "altitude", field, origin="native")
        sl.set_mpl_projection("InterruptedGoodeHomolosine")
        sl.set_log(field, False)
        return sl.plots[field].figure

    @requires_module("cartopy")
    @pytest.mark.mpl_image_compare(filename="geo_tranform_RotatedPole.png")
    def test_geo_tranform_RotatedPole(self):
        field = ("stream", "Density")
        sl = yt.SlicePlot(self.ds, "altitude", field, origin="native")
        sl.set_mpl_projection("RotatedPole")
        sl.set_log(field, False)
        return sl.plots[field].figure

    @requires_module("cartopy")
    @pytest.mark.mpl_image_compare(filename="geo_tranform_OSGB.png")
    def test_geo_tranform_OSGB(self):
        field = ("stream", "Density")
        sl = yt.SlicePlot(self.ds, "altitude", field, origin="native")
        sl.set_mpl_projection("OSGB")
        sl.set_log(field, False)
        return sl.plots[field].figure

    @requires_module("cartopy")
    @pytest.mark.mpl_image_compare(filename="geo_tranform_EuroPP.png")
    def test_geo_tranform_EuroPP(self):
        field = ("stream", "Density")
        sl = yt.SlicePlot(self.ds, "altitude", field, origin="native")
        sl.set_mpl_projection("EuroPP")
        sl.set_log(field, False)
        return sl.plots[field].figure

    @requires_module("cartopy")
    @pytest.mark.mpl_image_compare(filename="geo_tranform_Geostationary.png")
    def test_geo_tranform_Geostationary(self):
        field = ("stream", "Density")
        sl = yt.SlicePlot(self.ds, "altitude", field, origin="native")
        sl.set_mpl_projection("Geostationary")
        sl.set_log(field, False)
        return sl.plots[field].figure

    @requires_module("cartopy")
    @pytest.mark.mpl_image_compare(filename="geo_tranform_Gnomonic.png")
    def test_geo_tranform_Gnomonic(self):
        field = ("stream", "Density")
        sl = yt.SlicePlot(self.ds, "altitude", field, origin="native")
        sl.set_mpl_projection("Gnomonic")
        sl.set_log(field, False)
        return sl.plots[field].figure

    @requires_module("cartopy")
    @pytest.mark.mpl_image_compare(filename="geo_tranform_NorthPolarStereo.png")
    def test_geo_tranform_NorthPolarStereo(self):
        field = ("stream", "Density")
        sl = yt.SlicePlot(self.ds, "altitude", field, origin="native")
        sl.set_mpl_projection("NorthPolarStereo")
        sl.set_log(field, False)
        return sl.plots[field].figure

    @requires_module("cartopy")
    @pytest.mark.mpl_image_compare(filename="geo_tranform_SouthPolarStereo.png")
    def test_geo_tranform_SouthPolarStereo(self):
        field = ("stream", "Density")
        sl = yt.SlicePlot(self.ds, "altitude", field, origin="native")
        sl.set_mpl_projection("SouthPolarStereo")
        sl.set_log(field, False)
        return sl.plots[field].figure

    @requires_module("cartopy")
    @pytest.mark.mpl_image_compare(filename="geo_tranform_AlbersEqualArea.png")
    def test_geo_tranform_AlbersEqualArea(self):
        field = ("stream", "Density")
        sl = yt.SlicePlot(self.ds, "altitude", field, origin="native")
        sl.set_mpl_projection("AlbersEqualArea")
        sl.set_log(field, False)
        return sl.plots[field].figure

    @requires_module("cartopy")
    @pytest.mark.mpl_image_compare(filename="geo_tranform_AzimuthalEquidistant.png")
    def test_geo_tranform_AzimuthalEquidistant(self):
        field = ("stream", "Density")
        sl = yt.SlicePlot(self.ds, "altitude", field, origin="native")
        sl.set_mpl_projection("AzimuthalEquidistant")
        sl.set_log(field, False)
        return sl.plots[field].figure

    @requires_module("cartopy")
    @pytest.mark.mpl_image_compare(filename="geo_tranform_Sinusoidal.png")
    def test_geo_tranform_Sinusoidal(self):
        field = ("stream", "Density")
        sl = yt.SlicePlot(self.ds, "altitude", field, origin="native")
        sl.set_mpl_projection("Sinusoidal")
        sl.set_log(field, False)
        return sl.plots[field].figure

    @requires_module("cartopy")
    @pytest.mark.mpl_image_compare(filename="geo_tranform_NearsidePerspective.png")
    def test_geo_tranform_NearsidePerspective(self):
        field = ("stream", "Density")
        sl = yt.SlicePlot(self.ds, "altitude", field, origin="native")
        sl.set_mpl_projection("NearsidePerspective")
        sl.set_log(field, False)
        return sl.plots[field].figure

    @requires_module("cartopy")
    @pytest.mark.mpl_image_compare(
        filename="geo_tranform_LambertAzimuthalEqualArea.png"
    )
    def test_geo_tranform_LambertAzimuthalEqualArea(self):
        field = ("stream", "Density")
        sl = yt.SlicePlot(self.ds, "altitude", field, origin="native")
        sl.set_mpl_projection("LambertAzimuthalEqualArea")
        sl.set_log(field, False)
        return sl.plots[field].figure
