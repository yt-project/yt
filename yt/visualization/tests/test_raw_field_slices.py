import pytest

import yt
from yt.utilities.answer_testing.framework import data_dir_load_v2

try:
    ds = data_dir_load_v2("Laser/plt00015")
except Exception:

    pytest.skip(
        reason="dataset 'Laser/plt00015' isn't installed", allow_module_level=True
    )


class TestRawFieldSlice:
    @classmethod
    def setup_class(cls):
        cls.ds = ds
        fields = [
            ("raw", "Bx"),
            ("raw", "By"),
            ("raw", "Bz"),
            ("raw", "Ex"),
            ("raw", "Ey"),
            ("raw", "Ez"),
            ("raw", "jx"),
            ("raw", "jy"),
            ("raw", "jz"),
        ]
        cls.sl = yt.SlicePlot(cls.ds, "z", fields)
        cls.sl.set_log("all", False)
        cls.sl._setup_plots()

    @pytest.mark.mpl_image_compare(filename="raw_fields_slice_Bx.png")
    def test_raw_field_slice_Bx(self):
        return self.sl.plots["raw", "Bx"].figure

    @pytest.mark.mpl_image_compare(filename="raw_fields_slice_By.png")
    def test_raw_field_slice_By(self):
        return self.sl.plots["raw", "By"].figure

    @pytest.mark.mpl_image_compare(filename="raw_fields_slice_Bz.png")
    def test_raw_field_slice_Bz(self):
        return self.sl.plots["raw", "Bz"].figure

    @pytest.mark.mpl_image_compare(filename="raw_fields_slice_Ex.png")
    def test_raw_field_slice_Ex(self):
        return self.sl.plots["raw", "Ex"].figure

    @pytest.mark.mpl_image_compare(filename="raw_fields_slice_Ey.png")
    def test_raw_field_slice_Ey(self):
        return self.sl.plots["raw", "Ey"].figure

    @pytest.mark.mpl_image_compare(filename="raw_fields_slice_Ez.png")
    def test_raw_field_slice_Ez(self):
        return self.sl.plots["raw", "Ez"].figure

    @pytest.mark.mpl_image_compare(filename="raw_fields_slice_jx.png")
    def test_raw_field_slice_jx(self):
        return self.sl.plots["raw", "jx"].figure

    @pytest.mark.mpl_image_compare(filename="raw_fields_slice_jy.png")
    def test_raw_field_slice_jy(self):
        return self.sl.plots["raw", "jy"].figure

    @pytest.mark.mpl_image_compare(filename="raw_fields_slice_jz.png")
    def test_raw_field_slice_jz(self):
        return self.sl.plots["raw", "jz"].figure
