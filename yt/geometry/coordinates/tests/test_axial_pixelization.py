import pytest

from yt.testing import fake_amr_ds
from yt.utilities.answer_testing.answer_tests import axial_pixelization


@pytest.mark.answer_test
class TestAxialPixelization:
    @pytest.mark.usefixtures("hashing")
    def test_axial_pixelization(self, geom, axis):
        ds = fake_amr_ds(geometry=geom)
        ap = axial_pixelization(ds)
        if axis == "x_axis":
            ap = ap[0]
        elif axis == "y_axis":
            ap = ap[1]
        self.hashes.update({"axial_pixelization": ap})
