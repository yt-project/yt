import pytest

import yt
from yt.testing import requires_file
from yt.utilities.answer_testing.framework import data_dir_load_v2


class TestRawFieldSlice:
    @requires_file("Laser/plt00015")
    @classmethod
    def setup_class(cls):
        cls.ds = data_dir_load_v2("Laser/plt00015")
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
        cls.sl.render()

    @pytest.mark.parametrize(
        "fname", ["Bx", "By", "Bz", "Ex", "Ey", "Ez", "jx", "jy", "jz"]
    )
    @pytest.mark.mpl_image_compare
    def test_raw_field_slice(self, fname):
        return self.sl.plots["raw", "Bx"].figure
